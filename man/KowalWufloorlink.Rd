% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/starSSL.R
\name{KowalWufloorlink}
\alias{KowalWufloorlink}
\title{Empirical link b Kowal & Wu 2021 and a floor rounding
calculate the lower and upper bound by under Kowal & Wu 2021's transformation deterimination procedure}
\usage{
KowalWufloorlink(Y, distribution = "np", overall = TRUE, grid = NULL)
}
\arguments{
\item{Y}{the raw counting data}

\item{distribution}{which distribution to use?  must be one of
\itemize{
\item "np" (empirical CDF)
\item "pois" (moment-matched marginal Poisson CDF)
\item "neg-bin" (moment-matched marginal Negative Binomial CDF)
}}

\item{overall}{whether to find the transformation using a collected CDF or find it for each coordinate}

\item{grid}{the grid to calculate the inverse transformation}
}
\value{
A list of lower and upper bounds, as well as the empirical paramters used in that transformation, to be used for prediction.
}
\description{
Empirical link b Kowal & Wu 2021 and a floor rounding
calculate the lower and upper bound by under Kowal & Wu 2021's transformation deterimination procedure
}
