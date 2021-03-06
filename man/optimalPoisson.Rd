% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcpois.R
\name{optimalPoisson}
\alias{optimalPoisson}
\title{Kelly-fraction under geometric compensated Poisson dynamics}
\usage{
optimalPoisson(tt, a, b, lambda, rate = 0, N = 1000)
}
\arguments{
\item{tt}{maturity to simulate until}

\item{a}{jump size}

\item{b}{compensator size}

\item{lambda}{mean-rate of jumps}

\item{rate}{risk-free rate of return}

\item{N}{number of time-subintervals to take}
}
\value{
data.frame of time, stock, and the portfolio values.
}
\description{
{Generate a sample path of \eqn{\log X_t} of the Kelly portfolio
under a stock driven by a geometric compensated Poisson process.}
}
\details{
{The optimal fraction is analytically known through the equation
\eqn{\lambda/(r+b)-1/(e^a-1)} where the log dynamics follow
\eqn{X_t=aN_t-bt}, a scaled and compensated Poisson process. A basic
Euler-Maruyama scheme is then used to generate sample paths of both the stock
and the Kelly-portfolio.}
}
