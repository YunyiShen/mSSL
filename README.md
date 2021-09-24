# collection of multivariate Spike-and-Slab LASSO

A working repository collect multivariate Spike-and-Slab model for different type of responses centered on Gaussian. With effective posterior explorations. 

- For continuous data (unit tested)
    - multivariate SSL by [Deshpande et al. 2018](https://arxiv.org/abs/1708.08911), this uses a mean covariance parameterization of Gaussian model
    - chain graphical SSL by Shen, Solis-Lemus and Deshpande 2021+ (not published), uses a Gaussian chain graph parameterization
    - chain graphical VAR(1) SSL described in the same Shen, Solis-Lemus and Deshpande 2021+, uses a Gaussian chain graphical VAR(1) model for time series

- For binary data (unit tested)
    - multivariate probit SSL follows multivariate Spike-and-Slab LASSO by [Deshpande et al. 2018](https://arxiv.org/abs/1708.08911), with an EM algorithm account for the probit structure follows [Gessner et al. 2019](arxiv.org/abs/1910.09328)'s sampling methods. (tested)
    - the chain graphical version of the above

- For counting
    - multivariate STAR SSL (currently only allow fixed links) follows multivariate Spike-and-Slab LASSO by [Deshpande et al. 2018](https://arxiv.org/abs/1708.08911), with EM algorithm described in [Kowal and Canale 2020](arxiv.org/abs/1906.11653)'s STAR model. (unit tested)
    - cg version of the above



It is not yet on CRAN and likely won't in the near future, to install, use:

```r
devtools::install_github("YunyiShen/mSSL")
```

Note we do not have documentation yet. 