#' Kelly-criterion for mixture diffusions
#'
#' @param t current time
#' @param s current price
#' @param rate the discounting rate/money-market account return
#' @param parameters matrix with rows probs, mus, sigmas for the components
#' @param spot initial price
#'
#' @description {Straightforward generalization of the classic Kelly-fraction for
#' functions of time and space.}
#' @return numeric
#' @export kellyMixtureDiffusion
kellyMixtureDiffusion <- function(t, s, rate, parameters, spot)
{
  probs <- parameters[1, ]
  mus <- parameters[2, ]
  sigmas <- parameters[3, ]

  # Drift and volatility coefficients for mixture diffusion
  mu <- function(t, s) findistr::drift_lvm(s, t, probs, mus, sigmas, spot)
  v <- function(t, s) findistr::volat_lvm(s, t, probs, mus, sigmas, spot)
  z <- (mu(t, s)-rate)/(v(t, s)^2)
  return(z)
}

#' Compute the log-growth rate under optimally controlled mixture diffusions
#'
#' @param t time horizon
#' @param s current price
#' @param spot initial price
#' @param rate risk-free rate
#' @param parameters parameters of mixture
#'
#' @description {Compute the entropy rate as a function of the fraction. This
#' involves solving a Feynman-Kac PDE with running cost the sum of the risk-free rate
#' and half the square of the market price of risk.}
#' @return numeric
#' @export entropyMixtureDiffusion
entropyMixtureDiffusion <- function(t, s, spot, rate, parameters)
{
  probs <- parameters[1, ]
  mus <- parameters[2, ]
  sigmas <- parameters[3, ]
  mu <- function(t, x) findistr::drift_lvm(x, t, probs, mus, sigmas, spot)
  volat <- function(t, x) findistr::volat_lvm(x, t, probs, mus, sigmas, spot)
  sharpe <- function(t, x) (mu(t, x)- rate)/volat(t, x)

  dynamics <- list(function(t, x) mu(t, x)*x,
                   function(t, x) volat(t, x)*x
  )
  problem <- list(function(t, x) 0,
                  function(t, x) rate+0.5*sharpe(t, x)^2,
                  function(x) 0
  )
  region <- c(t, 0, 2*spot)
  control <- list(N = 200, M = 200, variational = FALSE, output = "price", engine = "c++")
  v <- fkpde::solvePDE(dynamics, problem, region, control)
  return(v)
}


#' Simulate log-optimal strategy on mixture diffusions
#'
#' @param bankroll initial bankroll to invest
#' @param t time horizon to trade over
#' @param spot the initial stock price
#' @param rate the return of the bond
#' @param parameters matrix of parameters of mixture, see details
#'
#' @description {Simulate log-optimal strategy for mixture diffusion prices}
#' @details {The matrix must contain a row of probabilities, a row of means, and
#' a row of standard deviations.}
#' @return data.frame of solution
#' @export simulateMixtureDiffusion
simulateMixtureDiffusion <- function(bankroll, t, spot, rate, parameters)
{
  probs <- parameters[1, ]
  mus <- parameters[2, ]
  sigmas <- parameters[3, ]

  # Drift and volatility coefficients for mixture diffusion
  mu <- function(t, s) findistr::drift_lvm(s, t, probs, mus, sigmas, spot)
  volat <- function(t, s) findistr::volat_lvm(s, t, probs, mus, sigmas, spot)
  IC <- list(s = spot, x = bankroll)
  control <- function(t, s) kellyMixtureDiffusion(t, s, rate, parameters, spot)

  f <- list(function(t, s, x) mu(t, s)*s,
            function(t, s, x) (rate+(mu(t, s)-rate)*control(t, s))*x
  )
  g <- list(function(t, s, x) volat(t, s)*s,
            function(t, s, x) control(t, s)*volat(t, s)*x
  )
  z <- sdes::samplePathSystem(f, g, IC, NULL, t0 = 0, tn = t, n = 1000)
  sol <- z
  sol$s <- log(sol$s)-log(IC$s)
  sol$x <- log(sol$x)-log(IC$x)
  graphics::par(mfrow = c(1, 2))
  odeSolveR::plot_trajectories(sol, legend_names = c("Stock", "Portfolio"),
                               legend_size = 0.4, legend_loc = "topleft")
  odeSolveR::plot_phase_portrait(sol)
  return(z)
}


