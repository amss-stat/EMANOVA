# EMANOVA


## Overview
EMANOVA provides tools for testing the association between microbiome and biological factors.

## Key features
E-MANOVA (Ensemble multivariate analysis of variance using distance matrices) method addressing the limitations of Traditional PERMANOVA. Traditional PERMANOVA is not robust to distances and association signals, which can lead to power reduction in certain scenarios. Using the idea of ensemble learning, we take similarity matrix to the $r$-th power to construct base test and then combine multiple tests to construct ensemble test. Our test statistic demonstrates high power and robustness compared to other existing methods. We also use direct moment approximation and Pearson type III distribution to approximate the permutation null distribution, completely avoiding the computationally intensive permutation procedure. Finally, we utilize the Cauchy combination method to aggregate p-values from multiple distances, eliminating the need to pre-specify distance measure before analysis. 

## Installation
```r
remotes::install_github("amss-stat/E-MANOVA")
```

## Quick start
```r
library(EMANOVA)
library(GUniFrac)
data("throat.otu.tab")
data("throat.tree")
data("throat.meta")

####transfer the outcome variable into binary.
smoke_stat<-rep(0,60)
for (i in 1:60) {
  if(throat.meta$SmokingStatus[i]==throat.meta$SmokingStatus[1]){
    smoke_stat[i]<-0
  }
  else{
    smoke_stat[i]<-1
  }
}

####Calculate relative abundance.
row.sum<-apply(throat.otu.tab,1,sum)
std_out<-throat.otu.tab/row.sum

####calculated the distance, D.05 is weighted UniFrac distance,
####Du unweighted UniFrac distance, D.BC is Bray-Curtis distance.
unifracs = GUniFrac::GUniFrac(std_out, throat.tree, alpha = c(0,0.5,1))$unifracs
Du  = unifracs[,,"d_UW"]
D.05= unifracs[,,"d_0.5"]
D.BC= as.matrix(vegan::vegdist(std_out, method="bray"))

combine_distance<-abind::abind(D.05, Du, D.BC, along = 3)
EMANOVA(combine_distance, smoke_stat, confonding_stat = FALSE, r_vec = c(0.125,0.25,0.5,1,2))
```






