#' Kelly-fraction for Ito process with coefficient functions of time and space.
#'
#' @param t current time, numeric
#' @param s current price, numeric greater than zero
#' @param dynamics list of drift and volatility coefficient functions of \code{(t, x)}
#' @param rate risk-free rate of return
#' @param FUN boolean default to false to return  numeric, true to return function
#'
#' @description {The generalized Kelly-fraction for space-time dependent coefficient functions
#' in Ito processes.}
#' @return numeric or function
#' @export kellyItoProcess
kellyItoProcess <- function(t, s, dynamics, rate = 0, FUN = FALSE)
{
  if(s <= 0)
  {
    stop("Price 's' must be positive")
  }
  if(t < 0)
  {
    stop("Time 't' must be non-negative")
  }
  if(dynamics[[2]](t, s) == 0)
  {
    stop("volatility coefficient is zero at (t,s); Kelly-fraction is undefined")
  }
  f <- (dynamics[[1]](t, s)-rate)/(dynamics[[2]](t, s)^2)
  if(FUN)
  {
    return(function(tt, ss) (dynamics[[1]](tt, ss)-rate)/(dynamics[[2]](tt, ss)^2))
  } else
  {
    return(f)
  }

}

#' Log-growth rate under optimal fraction invested in Ito process
#'
#' @param t length of time
#' @param s current spot price
#' @param dynamics list of drift and volatility function
#' @param rate risk-free rate
#' @param resolution the time and space resolution for the PDE solver
#'
#' @description {Compute the entropy rate as a function of the fraction. This
#' involves solving a Feynman-Kac PDE with running cost the sum of the risk-free rate
#' and half the square of the market price of risk.}
#' @return numeric
#' @export entropyItoProcess
entropyItoProcess <- function(t, s, dynamics, rate = 0, resolution = c(200, 200))
{
  # setupPDE handles the time-reversion, so no need to here.
  mu <- function(tt, ss) dynamics$drift(tt, ss)
  volat <- function(tt, ss) dynamics$diffusion(tt, ss)
  sharpe <- function(tt, ss) (mu(tt, ss)- rate)/volat(tt, ss)

  priceDynamics <- list(drift = function(tt, ss) mu(tt, ss)*ss,
                   diffusion = function(tt, ss) volat(tt, ss)*ss
  )
  problem <- list(discount = function(tt, ss) 0,
                  running_cost = function(tt, ss) rate+0.5*sharpe(tt, ss)^2,
                  terminal_cost = function(x) 0
  )
  region <- c(0, 2*s)
  v <- pdes::feynman_kac(region, t, priceDynamics, problem, FALSE, resolution)
  return(v$u[resolution[1]+1, resolution[2]/2+1])

}

#' Simulate log-optimal strategy for a given Ito process
#'
#' @param bankroll initial bankroll to invest
#' @param t time horizon to trade over
#' @param spot the initial stock price
#' @param rate the return of the bond
#' @param dynamics list of drift and volatility function of \code{(t, s)}
#' @param n number of time-steps to use in sample-path simulation
#'
#' @description {Simulate log-optimal strategy for prices driven by an Ito process}
#' @details {The argument \code{dynamics} must be a named list of functions \code{drift(t,x)} and
#' \code{diffusion(t, x)}.}
#' @return data.frame of solution
#' @export optimalItoProcess
optimalItoProcess <- function(bankroll, t, spot, rate, dynamics, n = 1000)
{
  # Drift and volatility coefficients for mixture diffusion
  mu <- function(tt, ss) dynamics$drift(tt, ss)
  volat <- function(tt, ss) dynamics$diffusion(tt, ss)
  IC <- list(s = spot, x = bankroll)
  control <- kellyItoProcess(t, spot, dynamics, rate, TRUE)

  f <- list(function(t, s, x) mu(t, s)*s,
            function(t, s, x) (rate+(mu(t, s)-rate)*control(t, s))*x
  )
  g <- list(function(t, s, x) volat(t, s)*s,
            function(t, s, x) control(t, s)*volat(t, s)*x
  )
  z <- sdes::sde_system(IC, t0 = 0, tn = t, f, g, n = n)
  sol <- z
  sol$s <- log(sol$s)-log(IC$s)
  sol$x <- log(sol$x)-log(IC$x)
  graphics::par(mfrow = c(1, 2))
  odeSolveR::plot_trajectories(sol, legend_names = c("Stock", "Portfolio"),
                               legend_size = 0.4, legend_loc = "topleft")
  odeSolveR::plot_phase_portrait(sol)
  return(z)
}
