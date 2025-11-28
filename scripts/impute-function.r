impute.fun = function(x,margin,FUN,...) {
  idx = which(is.na(x), arr.ind=T)
  if (nrow(idx) > 0) {
    vs = apply(x,margin,FUN,...)
    x[idx] = vs[idx[,margin]]
  }
  x
}
