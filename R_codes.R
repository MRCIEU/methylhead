---
title: "DNA Methylation: Bisulfite Sequencing Report"
author: "Onur"
date: now
format: 
  html:
    code-fold: true
    toc: true
    number-sections: true
    css: styles.css
---

```{r}
#| echo: false
#| label: libraries
#| warning: false
#| results: hide
library(methylKit)
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(GenomicRanges)
library(genomation)
library(meffonym)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(data.table)
```

```{r}
#| echo: false
#| label: upload files from bismark results
#| warning: false

setwd("/user/work/ag24712/new_data/nf/next_f/results/results/")
file.vector <- list.files(pattern = "bismark\\.cov\\.gz", full.names = FALSE)
sample.ids <- gsub("_bismark_bt2_pe.deduplicated.bismark.cov.gz", "", basename(file.vector))
sample.ids<-as.list(sample.ids)
sample_info<-fread("Sample_info.txt")
file.list<-as.list(file.vector)
myobj <- methRead(file.list,
           sample.id=sample.ids,
           pipeline = "bismarkCoverage",                
           assembly="hg38",
           treatment=sample_info$treatment,
           mincov = 10)
## fig-cap: "Checking number of samples"
myobj
```


```{r}
#| echo: false
#| label: Get a histogram of the read coverage per sample
#| warning: false

## "coverage stats"
for (i in 1:length(myobj)) {
  coverage_stats <- getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE)
  cat("Sample", i, "coverage stats:\n")
  print(coverage_stats)
}
```

```{r}
#| echo: false
#| label: Filtering Normalization and get percent methylation matrix
#| warning: false
myobj.filt <- filterByCoverage(myobj,
                      lo.count=10,
                      lo.perc=NULL,
                      hi.count=NULL,
                      hi.perc=99.9)
## "Filtered Results"
myobj.filt 
myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")
meth <- unite(myobj.filt.norm, destrand=FALSE)
## "Normalizaion"
head(meth)
pm=percMethylation(meth)
pm=pm/100
```

```{r}
#| echo: false
#| label: Merge data and annotation
#| warning: false
meth_df<-data.frame(meth)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann4<-data.frame(ann450k)
meth_df<-cbind(meth_df,(pm))
merged_data <- left_join(meth_df, ann4, by = c("start" = "pos"))
data1<-merged_data[,c(1,2,3,17,18,19,20,23,5,8,11,14)]
final_data<-data1
final_data2<-na.omit(final_data)
```

```{r}
#| echo: false
#| label: Preparing Data
#| warning: false   
methylation<-final_data2[,c(4,5,6,7)]
rownames(methylation)<-(final_data2$Name)
meth_matrix<-data.matrix(methylation)
```

```{r}
#| echo: false
#| label: DNA methylation indices of exposure and phenotype (meffonym)
#| warning: false
ret <- meffonym.score(m2, "hillary")
age<-c(60,58,39,24)
cor(ret$score, age)
```

```{r}
#| echo: false
#| label: Correlation Plots
#| warning: false
cor_data<-cbind(age,ret$score)
colnames(cor_data)<-c("Age","Score")
p.mat <- cor_pmat(cor_data)
corr <- round(cor(cor_data), 1)
ggcorrplot(corr)
ggcorrplot(corr, hc.order = TRUE, outline.col = "white")
ggcorrplot(corr, hc.order = TRUE, type = "lower",outline.col = "white")
ggcorrplot(corr, hc.order = TRUE, type = "lower",lab = TRUE)
ggplot(cor_data, aes(x=Age, y=Score)) + geom_point()
## Add the regression line
ggplot(cor_data, aes(x=Age, y=Score)) +  geom_point()+geom_smooth(method=lm)
## Remove the confidence interval
ggplot(cor_data, aes(x=Age, y=Score)) + geom_point()+geom_smooth(method=lm, se=FALSE)
## Loess method
ggplot(cor_data, aes(x=Age, y=Score)) +  geom_point()+geom_smooth()
```
