#' Optimal fraction under CEV diffusions
#'
#' @param spot current spot price
#' @param rate risk-free rate
#' @param parameters vector of drift, volatility and cev term.
#'
#' @description {Generalization of Kelly-fraction to CEV model.}
#' @return numeric
#' @export kellyCEV
kellyCEV <- function(spot, rate, parameters)
{
  cevParameterCheck(parameters)
  mu <- parameters[1]
  volat <- parameters[2]
  cev <- parameters[3]
  return((mu-rate)/(volat*spot^(cev-1))^2)
}


#' Optimal growth-rate under CEV diffusions
#'
#' @param t time-horizon
#' @param spot current spot price
#' @param rate risk-free rate
#' @param parameters vector of drift, volatility and cev term.
#'
#' @description {Log return under Kelly-fraction in the CEV model.}
#' @return numeric
#' @export entropyCEV
entropyCEV <- function(t, spot, rate, parameters)
{
  cevParameterCheck(parameters)
  mu <- parameters[1]
  volat <- parameters[2]
  cev <- parameters[3]
  dynamics <- list(function(t, s) mu,
                   function(t, s) volat*s^(cev-1))
  z <- entropyItoProcess(t, spot, dynamics, rate)
  return(z)
}

#' Simulate log-optimal strategy on CEV diffusions
#'
#' @param bankroll initial bankroll to invest
#' @param t time horizon to trade over
#' @param spot the initial stock price
#' @param rate the return of the bond
#' @param parameters vector of cev parameters
#'
#' @description {Simulate log-optimal strategy for cev diffusion prices. A wrapper
#' to \code{optimalItoProcess}.}
#' @details {parameters must be a vector of drift, volatility and cev.}
#' @return data.frame of solution
#' @export optimalCEV
optimalCEV <- function(bankroll, t, spot, rate, parameters)
{

  cevParameterCheck(parameters)
  mu <- parameters[1]
  volat <- parameters[2]
  cev <- parameters[3]
  dynamics <- list(function(t, s) mu,
                   function(t, s) volat*s^(cev-1))
  z <- optimalItoProcess(bankroll, t, spot, rate, dynamics)
  return(z)
}

#' Input error handling for cev model
#'
#' @param parameters vector of three numbers, drift, volatility and cev
#'
#' @return error or null
cevParameterCheck <- function(parameters)
{
  if(length(parameters) != 3)
  {
    stop("parameters must be a vector of drift, volatility and constant of elasticity of variance")
  }
  if(parameters[3] > 1)
  {
    stop("cev (parameters[3]) must be less than one")
  }
  if(parameters[2] <= 0)
  {
    stop("volatility must be positive")
  }
}
