# multivariate probit Spike-and-Slab LASSO

A working repository developing a multivariate probit Spike-and-Slab LASSO for binary data. 

It generally follows multivariate Spike-and-Slab LASSO by [Deshpande et al. 2018](https://arxiv.org/abs/1708.08911), with an EM algorithm account for the probit structure follows [Gessner et al. 2019](arxiv.org/abs/1910.09328)'s sampling methods.

It is not on CRAN and likely won't in the near future, to install, use:

```r
devtools::install_github("YunyiShen/mpSSL")
```

Note we do not have documentation yet. 