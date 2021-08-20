#' Kelly-fraction for European options under Black-Scholes dynamics
#'
#' @param mu real drift of the underlying stock price
#' @param volat volatility of the stock price and volatility used for the option pricing
#' @param rate the risk-neutral rate
#' @param spot the underlying spot price
#' @param premium the premium paid for the option
#' @param delta the delta of the option, rate of change with respect to spot-price
#'
#' @description {A Kelly-fraction generalization for trading a portfolio consiting of
#' a European option priced under Black-Scholes dynamics with an assumed real-volatility,
#' and a cash-account or bond at the risk-neutral rate. TWo functions are provided in this package,
#' one for computing it specified by the exact values, and one computed via specifying the
#' option via strike, expiry, etc and computing the option value and delta via a pricing method.}
#' @return numeric
#' @export kellyEuroBS1
kellyEuroBS1 <- function(mu, volat, rate, spot, premium, delta)
{
  excessReturn <- mu-rate
  dollarDelta <- delta*spot
  numerator <- excessReturn*premium
  denominator <- (volat^2)*dollarDelta
  kellyRatio <- numerator/denominator
  return(kellyRatio)
}

#' Kelly-fraction for European options under Black-Scholes dynamics
#'
#' @param strike strike price of the option
#' @param expiry time to maturity of option in trading years
#' @param mu real drift of the underlying stock price
#' @param volat volatility of the stock price and volatility used for the option pricing
#' @param rate the risk-neutral rate
#' @param spot the underlying spot price
#' @param type type of option, "call" or "put"
#' @param N number of sub-intervals in time-grid
#' @param M number of sub-intervals in space-grid
#'
#' @description {A Kelly-fraction generalization for trading a portfolio consiting of
#' a European option priced under Black-Scholes dynamics with an assumed real-volatility,
#' and a cash-account or bond at the risk-neutral rate. TWo functions are provided in this package,
#' one for computing it specified by the exact values, and one computed via specifying the
#' option via strike, expiry, etc and computing the option value and delta via a pricing method.}
#' @return numeric
#' @export kellyEuroBS2
kellyEuroBS2 <- function(strike, expiry, mu, volat, rate, spot, type, N = 100, M = 100)
{
  w <- pricing::bs_greeks(strike, expiry, spot, volat, rate, type, N, M, american = FALSE)
  premium <- w[2]
  delta <- w[3]
  return(kellyEuroBS1(mu, volat, rate, spot, premium, delta))
}
