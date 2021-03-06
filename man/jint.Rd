% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merton.R
\name{jint}
\alias{jint}
\title{Jumps integral in Kelly criterion for Merton dynamics}
\usage{
jint(a, g, alpha, beta, subdivisions = 200)
}
\arguments{
\item{a}{the control value}

\item{g}{the function to compute expectation of \eqn{J = e^Y-1}}

\item{alpha}{the mean jump size}

\item{beta}{the standard deviation of the jump size}

\item{subdivisions}{number of subdivisions in the numerical integration}
}
\value{
numeric
}
\description{
{Compute various integrals against a density function for the
log jump size.}
}
