% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simu_data.R
\name{rCVARm}
\alias{rCVARm}
\title{Sampling from Conditional Vector Autoregression (1) model, another simple Gaussian chain graph model}
\usage{
rCVARm(n, B, Sigma)
}
\arguments{
\item{n}{number of samples}

\item{B}{regression coefficients}

\item{Sigma}{Covariance matrix}
}
\value{
a matrix of n rows
}
\description{
Take random samples following VCAR(1) model, i.e. Y_t+1~N(Y_tBSigma,Sigma)
}
