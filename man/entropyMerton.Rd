% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merton.R
\name{entropyMerton}
\alias{entropyMerton}
\title{Kelly-criterion entropy under Merton's jump diffusion}
\usage{
entropyMerton(
  param,
  rate = 0,
  iterations = 500,
  subdivisions = 500,
  tol = 10^-6
)
}
\arguments{
\item{param}{vector of parameters defining the jump-diffusion dynamics. See details}

\item{rate}{the risk-free rate, or money-market account interest rate}

\item{iterations}{number of iterations to use in the Newton-Raphson method for finding the optimal control as a root}

\item{subdivisions}{number of subdivisions to use in numerical integrations}

\item{tol}{tolerance for the zero to be found.}
}
\value{
list
}
\description{
{Compute the entropy, i.e. growth rate of the Kelly allocation
under Merton's jump diffusion. }
}
