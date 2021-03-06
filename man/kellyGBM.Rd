% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometricBrownianMotion.R
\name{kellyGBM}
\alias{kellyGBM}
\title{The Kelly-fraction for a stock following GBM dynamics}
\usage{
kellyGBM(drift, volat, rate = 0, restraint = NULL)
}
\arguments{
\item{drift}{the mean drift rate of the GBM}

\item{volat}{the volatility of the GBM}

\item{rate}{the risk-free rate of the bond}

\item{restraint}{a number between 0 and 1, null for (possibly) leveraged Kelly}
}
\value{
numeric
}
\description{
{The famous Kelly-fraction, the amount to invest
in a risky stock following a geometric Brownian motion.}
}
