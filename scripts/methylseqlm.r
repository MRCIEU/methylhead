library(parallel)

#' Performs a modified version of seqlm
#' (https://github.com/raivokolde/seqlm/)
#' DMR analysis to methyl-seq data but in a slightly
#' more flexible and efficient way.
#' Flexibility is permitted by allowing coveriates
#' in the analysis, unlike the original seqlm.
#' Computationl efficiency is improved by
#' meta-analysing model coefficients
#' rather than refitting models for multi-CpG DMR candidates.
#' Memory requirements are reduced by
#' performing the dynamic programming in a way
#' that slides the saved statistics matrix across
#' the genome dropping statistics that won't influence
#' current or later DMR decisions (based
#' on a maximum DMR size). 

fit_linear_models = function(y,x) {
  ## fit a linear model for each outcome variable in y
  ## y: outcomes, rows=outcomes, cols=samples
  ## x: design matrix, rows=samples, cols=covariates
  ## returns model fit coefficients and standard errors
  ##    (rows=outcomes in y,cols=covariates in x)
  stopifnot(ncol(y) == nrow(x))
  stats = lm.fit(x,t(y))

  ## least squares coefficients for covariate
  coefs = stats$coefficients
  rownames(coefs) = colnames(x)
  colnames(coefs) = rownames(y)

  ## standard errors of coefficients for each covariate
  rss = colSums(stats$residuals^2)
  df = nrow(x)-ncol(x)
  sigma2 = rss/df
  xx_inv = diag(solve(crossprod(x)))
  ses = sqrt(matrix(xx_inv,ncol=1) %*% matrix(sigma2,nrow=1))
  dimnames(ses) = dimnames(coefs)

  list(coef=t(coefs), se=t(ses))
}

calc_residuals = function(y,x,coefs) {
  ## calculate differences between each outcome in y
  ## and values predicted by the model
  ##
  ## y: outcomes, rows=outcomes, cols=samples
  ## x: model design matrix, rows=samples, cols=covariates
  ## coefs: model coefficients, one for each model covariate
  ## returns matrix of residuals (rows=outcomes,cols=samples)
  stopifnot(ncol(x) == length(coefs))
  stopifnot(nrow(x) == ncol(y))
  preds = as.vector(x %*% matrix(coefs,ncol=1))
  t(apply(y,1,function(out) out-preds))
}

calc_sigma = function(rss,n,df,n0=0,sigma0=0) {
  ## calculate variance of residual sum of squares
  ## rss: residual sum of squares
  ## n: number of samples
  ## df: residual degrees of freedom
  ## sigma0: prior for sigma
  ## n0: weight for prior
  sigma = sqrt(rss/df)
  if (sigma0 > 0)
    sigma = sqrt((n0*sigma0^2 + n*sigma^2)/(n0+n-2))
  sigma
}

calc_sigma_prior = function(y,x) {
  ## calculate prior for sigma for a small random selection of outcomes
  ## y: outcomes, rows=outcomes, cols=samples
  ## x: model design matrix, rows=samples, cols=covariates
  stopifnot(ncol(y) == nrow(x))
  if (nrow(y) > 100)
    y = y[sample(1:nrow(y), 100),]
  n = ncol(y)
  rank = ncol(x)
  df = n - rank
  fits = fit_linear_models(y,x)
  mean(sapply(1:nrow(y), function(i) {
    res = calc_residuals(y[i,,drop=F], x, fits$coef[i,])
    calc_sigma(sum(res^2),n,df)
  }))
}

calc_loglike = function(y,x,fits,n0,sigma0) {
  ## calculate log-likelihood of all outcomes in y given x
  ## y: outcomes, rows=outcomes, cols=samples
  ## x: model design matrix, row=samples, cols=covariates
  ## fits: fit_linear_models(y,x) output
  stopifnot(ncol(y) == nrow(x))
  stopifnot(nrow(fits$coef) == nrow(y))
  stopifnot(ncol(fits$coef) == ncol(x))
  ## meta-analyse associations between x and each outcome in y
  wt = 1/fits$se^2
  coefs = colSums(fits$coef*wt)/colSums(wt)
  ## calculate log-likelihood of the meta-analysed model fit
  res = calc_residuals(y,x,coefs)
  n = ncol(y)*nrow(y)
  rank = ncol(x)
  df = n - rank
  rss = sum(res^2)
  sigma = calc_sigma(rss,n,df,n0,sigma0)
  if (!is.na(sigma) && sigma > 0)
    sum(log2(dnorm(x=res, mean=0, sd=sigma)))
  else
    0
}

segment = function(meth,pos,x,max_dmr,n0=1,m0=10,alpha=2) {
  ## segment the methylome to minimize description length in terms of x
  ## meth: methylation matrix (rows=CpG sites, cols=samples)
  ## pos: genomic position of each CpG site
  ##      (pos should be sorted and match rows of meth)
  ## x: design matrix (rows=samples, cols=covariates)
  ## max_dmr: length in bases of the maximum DMR
  ## n0: weight for sigma prior
  ## m0,alpha: parameters penalizing model complexity
  if (length(pos) == 1)
    return(cbind(start=1,end=1))
  stopifnot(length(pos) > 1)
  stopifnot(length(pos) == nrow(meth))
  stopifnot(ncol(meth) == nrow(x))
  stopifnot(!is.unsorted(pos))
  n = ncol(meth)
  p = nrow(meth)
  rank = ncol(x)
  ## calculate priors
  sigma0 = calc_sigma_prior(meth,x)
  n0 = n*n0
  m0 = n*m0
  ## set penalty for each candidate segment
  penalty = alpha*(rank*log2(m0)/2 + log2(max_dmr))
  ## fit model for each outcome in y
  fits = fit_linear_models(meth,x)

  ## prepare variables for dynamic programming
  best_cost = rep(Inf,p) ## for each genomic position, minimum description length (MDL)
  best_start = rep(NA,p) ## for each genomic position, start of final segment
  ## that yeilds MDL from methylome start to that position
  min_start = 1          ## for current genomic position, how far to look back
  ## to satisfy max_dmr

  ## perform dynamic programming starting at the beginning of the methylome
  ## i.e. progressively fill in best_cost and best_start starting at beginning
  for (end in 1:p) { ## for each CpG site
    ## determine how far to look back given max_dmr
    starts = min_start:end
    min_start = starts[which(pos[starts] >= pos[end]-max_dmr)[1]]
    ## calculate contribution of each segment ending at the current position to MDL
    costs = sapply(min_start:end, function(start) {
      ## for the segment from start..end
      ## just consider associations and methylation from position start..end
      fits = list(coef=fits$coef[start:end,,drop=F], se=fits$se[start:end,,drop=F])
      meth = meth[start:end,,drop=F]
      ## calculate contribution of segment to MDL
      cost = -calc_loglike(meth,x,fits,n0,sigma0) + penalty
      ## if start is not the beginning
      ## then MDL for this scheme includes MDL from 1..(start-1)
      if (start > 1)
        cost = cost + best_cost[start-1]
      cost
    })
    ## given all schemes end at the current position, choose the minimum
    idx = which.min(costs)
    best_cost[end] = costs[idx]
    best_start[end] = min_start + idx - 1
  }

  ## mdl for the whole methylome corresponds to last one calculated
  ## identify all segments by working backward to the beginning
  is_best = rep(F,length(best_start))
  idx = length(best_start)
  while (idx > 1) {
    idx = best_start[idx]
    is_best[idx] = T
    idx = idx - 1
  }
  starts = which(is_best)
  ends = c(starts[-1] - 1, length(is_best))
  cbind(start=starts,end=ends)
}

is.genome.ordered = function(chr,pos) {
  stopifnot(length(chr) == length(pos))
  chrs = unique(chr)
  ((sum(tail(chr,-1) != head(chr,-1)) == length(chrs) - 1)
    && all(tail(pos,-1) > head(pos,-1) | tail(chr,-1) != head(chr,-1)))
}

methylseqlm.segment = function(meth,chr,pos,x,max_dmr,n0=1,m0=10,alpha=2) {
  stopifnot(is.matrix(meth))
  stopifnot(is.matrix(x))
  stopifnot(!any(is.na(x)))
  stopifnot(!any(is.na(meth)))
  stopifnot(length(pos) == length(chr))
  stopifnot(length(pos) == nrow(meth))
  stopifnot(ncol(meth) == nrow(x))
  stopifnot(is.genome.ordered(chr,pos))

  ## ensure that x^tx is invertible before continuing
  xx_inv = solve(crossprod(x))  

  ## segment the genome
  segments = mclapply(unique(chr), function(cur) {
    idx = which(chr==cur)
    segments = segment(meth[idx,,drop=F],pos[idx],x,max_dmr,n0,m0,alpha)
    segments[,"start"] = idx[segments[,"start"]]
    segments[,"end"] = idx[segments[,"end"]]
    data.frame(
      chr=cur,
      start=pos[segments[,"start"]],
      end=pos[segments[,"end"]],
      idx = segments)
  })
  segments = do.call(rbind, segments)
  rownames(segments) = with(segments, paste0(chr,":",start,"-",end))
  segments
}

methylseqlm.merge = function(meth,segments,coverage=NULL) {
  stopifnot(is.matrix(meth))
  stopifnot(is.data.frame(segments))
  stopifnot(all(c("idx.start","idx.end") %in% colnames(segments)))
  stopifnot(all(segments$idx.start <= segments$idx.end))
  stopifnot(all(1 <= segments$idx.start))
  stopifnot(all(segments$idx.end <= nrow(meth)))
  if (is.null(coverage))
    coverage = matrix(1,ncol=ncol(meth),nrow=nrow(meth))
  stopifnot(is.matrix(coverage))
  stopifnot(nrow(coverage) == nrow(meth) && ncol(coverage) == ncol(meth))
  meth.merged = t(sapply(1:nrow(segments), function(i) {
    idx = segments[i,"idx.start"]:segments[i,"idx.end"]
    mcounts = colSums(meth[idx,,drop=F]*coverage[idx,,drop=F],na.rm=T)
    totals = colSums(coverage[idx,,drop=F],na.rm=T)
    mcounts/totals
  }))
  colnames(meth.merged) = colnames(meth)
  rownames(meth.merged) = rownames(segments)
  meth.merged
}


