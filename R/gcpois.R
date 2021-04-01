#' Kelly-fraction under geometric compensated Poisson dynamics
#'
#' @param a jump size
#' @param b compensator size
#' @param lambda mean-rate of jumps
#' @param rate risk-free rate of return
#'
#' @description {Computes the optimal allocation for a risky stock
#' driven by a compensated Poisson process in terms of its log dynamics.}
#' @details {The optimal fraction is analytically known through the equation
#' \eqn{\lambda/(r+b)-1/(e^a-1)} where the log dynamics follow
#' \eqn{X_t=aN_t-bt}, a scaled and compensated Poisson process.}
#' @return numeric real number
#' @export kellyPoisson
kellyPoisson <- function(a, b, lambda, rate = 0)
{
  alpha <- lambda/(rate+b)-1/(exp(a)-1)
  return(alpha)
}

#' Kelly-fraction under geometric compensated Poisson dynamics
#'
#' @param tt maturity to simulate until
#' @param a jump size
#' @param b compensator size
#' @param lambda mean-rate of jumps
#' @param rate risk-free rate of return
#' @param N number of time-subintervals to take
#'
#' @description {Generate a sample path of \eqn{\log X_t} of the Kelly portfolio
#' under a stock driven by a geometric compensated Poisson process.}
#' @details {The optimal fraction is analytically known through the equation
#' \eqn{\lambda/(r+b)-1/(e^a-1)} where the log dynamics follow
#' \eqn{X_t=aN_t-bt}, a scaled and compensated Poisson process. A basic
#' Euler-Maruyama scheme is then used to generate sample paths of both the stock
#' and the Kelly-portfolio.}
#' @return data.frame of time, stock, and the portfolio values.
#' @export optimalPoisson
optimalPoisson <- function(tt, a, b, lambda, rate = 0, N = 1000)
{
  alpha <- kellyPoisson(a, b, lambda, rate)
  k <- tt/N
  x <- seq(0, tt, length.out = N+1)
  y <- matrix(0, nrow = N+1)
  y[1] <- 0
  s <- matrix(0, nrow = N+1)
  s[1] <- 0
  for(i in 2:(N+1))
  {
    n <- stats::rpois(1, lambda = lambda*k)
    # The stock process
    s[i] <- s[i-1] + a*n-b*k
    # The controlled process (portfolio)
    y[i] <- y[i-1]+log(1+alpha*(exp(a)-1))*n+(rate-alpha*(rate+b))*k


  }
  graphics::par(mfrow = c(1, 1))
  plot(x, s, type = "l", ylim = c(min(y, s), max(y, s)))
  graphics::lines(x, y, col = "red")
  graphics::legend("topleft", legend = c("stock", "portfolio"),
         col = c("black", "red"), lty = 1, cex = 0.5)
  return(data.frame(t = tt, s = s, x = y))
}
