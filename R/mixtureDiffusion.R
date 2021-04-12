#' Kelly-criterion for mixture diffusions
#'
#' @param t current time
#' @param s current price
#' @param spot initial price
#' @param rate the discounting rate/money-market account return
#' @param param matrix with rows probs, mus, sigmas for the components
#' @param FUN boolean whether to return numeric or function
#'
#' @description {Straightforward generalization of the classic Kelly-fraction for
#' functions of time and space, specialized to mixture coefficients. A wrapper to \code{kellyItoProcess}.}
#' @return numeric
#' @export kellyMixtureDiffusion
kellyMixtureDiffusion <- function(t, s, spot, rate, param, FUN = FALSE)
{
  mixtureParameterCheck(param)
  probs <- param[1, ]
  mus <- param[2, ]
  sigmas <- param[3, ]

  # Drift and volatility coefficients for mixture diffusion
  dynamics <- list(drift = function(t, s) findistr::mixture_drift(t, log(s/spot), param),
                   diffusion = function(t, s) findistr::mixture_vol(t, log(s/spot), param)
  )
  z <- kellyItoProcess(t, s, dynamics, rate, FUN)
  return(z)
}

#' Compute the log-growth rate under optimally controlled mixture diffusions
#'
#' @param t time horizon
#' @param s current price
#' @param spot initial price
#' @param rate risk-free rate
#' @param param parameters of mixture
#' @param resolution time and space resolution for the PDE solver
#'
#' @description {A special case wrapper to \code{entropyItoProcess}. Compute the entropy rate as a function of the fraction. This
#' involves solving a Feynman-Kac PDE with running cost the sum of the risk-free rate
#' and half the square of the market price of risk.}
#' @return numeric
#' @export entropyMixtureDiffusion
entropyMixtureDiffusion <- function(t, s, spot, rate, param, resolution = c(200, 200))
{
  mixtureParameterCheck(param)
  probs <- param[1, ]
  mus <- param[2, ]
  sigmas <- param[3, ]
  dynamics <- list(drift = function(t, s) findistr::mixture_drift(t, log(s/spot), param),
                   diffusion = function(t, s) findistr::mixture_vol(t, log(s/spot), param)
  )
  z <- entropyItoProcess(t, s, dynamics, rate, resolution)
  return(z)
}


#' Simulate log-optimal strategy on mixture diffusions
#'
#' @param bankroll initial bankroll to invest
#' @param t time horizon to trade over
#' @param spot the initial stock price
#' @param rate the return of the bond
#' @param param matrix of parameters of mixture, see details
#'
#' @description {Simulate log-optimal strategy for mixture diffusion prices. A wrapper
#' to \code{optimalItoProcess}.}
#' @details {The matrix must contain a row of probabilities, a row of means, and
#' a row of standard deviations.}
#' @return data.frame of solution
#' @export optimalMixtureDiffusion
optimalMixtureDiffusion <- function(bankroll, t, spot, rate, param)
{
  mixtureParameterCheck(param)
  probs <- param[1, ]
  mus <- param[2, ]
  sigmas <- param[3, ]
  dynamics <- list(drift = function(t, s) findistr::mixture_drift(t, log(s/spot), param),
                   diffusion = function(t, s) findistr::mixture_vol(t, log(s/spot), param)
  )
  z <- optimalItoProcess(bankroll, t, spot, rate, dynamics)
  return(z)
}

#' Input error handling for mixture model
#'
#' @param param three rows containing probabilities, drifts, volatilities
#'
#' @return error or null
mixtureParameterCheck <- function(param)
{
  if(nrow(param) != 3)
  {
    stop("param must be three rows of probabilities, drifts, and volatilities")
  }
}

