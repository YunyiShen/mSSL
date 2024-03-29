% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpSSL.R
\name{mpSSL}
\alias{mpSSL}
\title{Multivariate probit Spike-and-Slab LASSO}
\usage{
mpSSL(
  Y,
  X,
  condexp = FALSE,
  cg = FALSE,
  lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
  xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
  theta_hyper_params = c(1, ncol(X) * ncol(Y)),
  eta_hyper_params = c(1, ncol(Y)),
  diag_penalty = FALSE,
  max_iter = 500,
  eps = 0.001,
  s_max_condition = 10 * nrow(X),
  obj_counter_max = 5,
  verbose = FALSE,
  nrep = 200,
  nskp = 1
)
}
\arguments{
\item{Y}{response matrix}

\item{X}{design matrix}

\item{condexp}{bool, whether to do the fast conditional posterior exploration (dcpe), default is FALSE for doing the dynamic posterior exploration (dpe)}

\item{cg}{bool, whether use chain graphical parameterization in latent normal, default FALSE, meaning using mean-cov parameterization, be careful about the interpretation of each.}

\item{lambdas}{hyperparameters to be explored by the algorithm, penalty on B}

\item{xis}{hyperparameters to be explored by the algorithm, penalty on Omega}

\item{theta_hyper_params}{hyperparameter to be set, prior on spike weight of B}

\item{eta_hyper_params}{hyperparameter to be set, prior on spike weight of Omega}

\item{diag_penalty}{bool, whether to penalize the diagonal, default no}

\item{max_iter}{maximum iterations for the EM algorithm}

\item{eps}{tolerance for convergence}

\item{s_max_condition}{only used in dpe, maximum tolerance for the condition number}

\item{obj_counter_max}{only used in dpe, maximum number of couting the objective function}

\item{verbose}{bool, whether to print intermidate notes}

\item{nrep}{number of sample to take during E step of the latent normal}

\item{nskp}{one point taken per nskp samples during E step of the latent normal}
}
\value{
A list with dynamic exploration result, point estimates are in `$Omega` and `$B`.
}
\description{
Main posterior exploration algorithm for multivariate probit spike-and-slab LASSO.
}
