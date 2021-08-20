#' The optimal chance of reaching a wealth goal by a terminal date
#'
#' @param t current time
#' @param x initial wealth level
#' @param b target wealth level
#' @param maturity terminal date to hit the goal by
#' @param mu the drift of the stock
#' @param volat the volatility of the stock
#' @param rate the risk-free rate from a cash-account, bond, etc
#'
#' @description {The optimal value from maximizing the probability of the terminal wealth
#' of a GBM-cash portfolio being greater than a fixed level, conditional on the initial wealth.}
#' @return numeric, matrix, vector
#' @export wealthGoalChance
wealthGoalChance <- function(t, x, b, maturity, mu, volat, rate)
{
  n <- length(t)
  m <- length(x)
  v <- matrix(0, n, m)
  lambda <- (mu-rate)/volat
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      if(x[j] > b*exp(-rate*(maturity-t[i])))
      {
        stop("initial wealth must be less than discounted goal")
      }
      z <- sqrt(2)*sqrt(log(b/x[j])-rate*(maturity-t[i]))-lambda*sqrt(maturity-t[i])
      v[i, j] <- 1-stats::pnorm(z)
    }
  }
  return(v)
}

#' The optimal control for reaching a wealth goal by a terminal date
#'
#' @param t current time
#' @param x initial wealth level
#' @param b target wealth level
#' @param maturity terminal date to hit the goal by
#' @param mu the drift of the stock
#' @param volat the volatility of the stock
#' @param rate the risk-free rate from a cash-account, bond, etc
#'
#' @description {The optimal control for maximizing the probability of the terminal wealth
#' of a GBM-cash portfolio being greater than a fixed level, conditional on the initial wealth.}
#' @return numeric, matrix, vector
#' @export wealthGoalControl
wealthGoalControl <- function(t, x, b, maturity, mu, volat, rate)
{

  n <- length(t)
  m <- length(x)
  v <- matrix(0, n, m)
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      if(x[j] > b*exp(-rate*(maturity-t[i])))
      {
        stop("initial wealth must be less than discounted goal")
      }
      v[i, j] <- sign(mu-rate)*(sqrt(2)/(volat*sqrt(maturity-t[i])))*sqrt(log(b/x[j])-rate*(maturity-t[i]))
    }
  }
  return(v)
}

#' The sub-optimal log-growth obtained from reaching a wealth goal by a terminal date
#'
#' @param t current time
#' @param x initial wealth level
#' @param b target wealth level
#' @param maturity terminal date to hit the goal by
#' @param mu the drift of the stock
#' @param volat the volatility of the stock
#' @param rate the risk-free rate from a cash-account, bond, etc
#'
#' @description {The sub-optimal log-growth from maximizing the probability of the terminal wealth
#' of a GBM-cash portfolio being greater than a fixed level, conditional on the initial wealth. This will
#' always be less than the entropy from the pure Kelly case.}
#' @return numeric, matrix, vector
#' @export wealthGoalEntropy
wealthGoalEntropy <- function(t, x, b, maturity, mu, volat, rate)
{
  lambda <- (mu-rate)/volat
  n <- length(t)
  m <- length(x)
  v <- matrix(0, n, m)
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      if(x[j] > b*exp(-rate*(maturity-t[i])))
      {
        stop("initial wealth must be less than discounted goal")
      }
      v[i,j] <- 2*rate-log(b/x[j])/(maturity-t[i])+sqrt(2)*lambda*sqrt(log(b/x[j])/(maturity-t[i])-rate)
    }
  }
  return(v)

}
