#' Optimized log-growth for Heston dynamics
#'
#' @param v0 current volatility state
#' @param tt time-horizon
#' @param param the vector \code{c(kappa, theta, xi, mu)}
#' @param rate the risk-free rate
#' @param resolution the time and space interval resolution
#'
#' @description {Compute the log-growth under the Kelly-criterion
#' for the stochastic Heston volatility model.}
#' @return numeric
#' @export entropyHeston
entropyHeston <- function(v0, tt, param, rate = 0, resolution = c(100, 100))
{
  kappa <- param[1]
  theta <- param[2]
  xi <- param[3]
  mu <- param[4]
  drift <- function(t, x) kappa*(theta-x)
  diffusion <- function(t, x) xi*sqrt(x)
  dynamics <- list(drift = drift, diffusion = diffusion)
  problem <- list(discount = function(t, x) 0,
                  running_cost = function(t, x) (mu-rate)^2/x+rate,
                  terminal_cost = function(x) 0
  )

  region <- c(0, 2*v0)
  w <- pdes::feynman_kac(region, tt, dynamics, problem, FALSE, resolution)
  return(w$u[resolution[1]+1, resolution[2]/2+1])
}


#' Simulate a sample-path of the Kelly-portfolio under Heston dynamics
#'
#' @param x0 initial wealth to invest
#' @param s0 initial spot price of the stock
#' @param v0 initial (hidden) stochastic volatility level
#' @param tt time horizon to simulate a path over
#' @param param vector of \code{kappa, theta, xi, mu, rho}
#' @param rate risk-free rate of cash-account, money-market account, or bond
#' @param n number of time sub-intervals in sample-path
#' @param plotG boolean for plotting the portfolio in the function body
#'
#' @description {A basic Euler-Maruyama scheme for simulating a sample path of
#' a stochastic volatility price model and its corresponding Kelly-strategy.}
#' @return data.frame of time, volatility, stock-price, and portfolio-value
#' @export optimalHeston
optimalHeston <- function(x0, s0, v0, tt, param, rate = 0, n = 1000, plotG = TRUE)
{
  h <- tt/n
  x <- matrix(0, nrow = n+1)
  s <- matrix(0, nrow = n+1)
  v <- matrix(0, nrow = n+1)
  x[1] <- x0
  s[1] <- s0
  v[1] <- v0
  kappa <- param[1]
  theta <- param[2]
  xi <- param[3]
  mu <- param[4]
  rho <- param[5]
  z <- findistr::rcornorm(n, rho = rho)
  for(i in 1:n)
  {
    v[i+1] <- v[i] + kappa*(theta-v[i])*h + xi*sqrt(v[i])*sqrt(h)*z[i, 1]
    s[i+1] <- s[i] + mu*s[i]*h + sqrt(v[i])*s[i]*sqrt(h)*z[i, 2]
    x[i+1] <- x[i] + (rate + ((mu-rate)^2) / v[i])*x[i]*h + ((mu-rate)/sqrt(v[i]))*x[i]*sqrt(h)*z[i, 2]
  }
  output <- data.frame(t = (0:n)*h, v = v, s = s, x = x)
  if(plotG)
  {
    plot(output$t, output$x, type = "l", main = "Heston Kelly-portfolio")
  }
  return(output)
}
