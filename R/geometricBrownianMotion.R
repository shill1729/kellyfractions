#' The Kelly-fraction for a stock following GBM dynamics
#'
#' @param drift the mean drift rate of the GBM
#' @param volat the volatility of the GBM
#' @param rate the risk-free rate of the bond
#' @param restraint a number between 0 and 1, null for (possibly) leveraged Kelly
#'
#' @description {The famous Kelly-fraction, the amount to invest
#' in a risky stock following a geometric Brownian motion.}
#' @return numeric
#' @export kellyGBM
kellyGBM <- function(drift, volat, rate = 0, restraint = NULL)
{
  f <- 0
  if(is.null(restraint))
  {
    f <- (drift-rate)/(volat^2)

  } else if (abs(drift-rate) <= restraint*volat^2)
  {
    f <- (drift-rate)/(volat^2)


  } else{
    f <- restraint
  }
  names(f) <- "optimalFraction"
  return(f)

}

#' The optimal log-growth rate under univariate GBM
#'
#' @param drift the mean drift rate of the GBM
#' @param volat the volatility of the GBM
#' @param rate the risk-free rate of the bond
#'
#' @description {The optimal log growth rate under the Kelly-criterion for
#' geometric Brownian motion.}
#' @return numeric
#' @export entropyGBM
entropyGBM <- function(drift, volat, rate = 0)
{
  lambda <- (drift-rate)/volat
  g <- rate+0.5*lambda^2
  return(g)
}


#' Simulate log-optimal strategy on mixture diffusions
#'
#' @param bankroll initial bankroll to invest
#' @param t time horizon to trade over
#' @param spot the initial stock price
#' @param rate the return of the bond
#' @param parameters drift and volatility
#'
#' @description {Simulate log-optimal strategy for prices under geometric Brownian motion. A wrapper
#' to \code{optimalItoProcess}.}
#' @details {The parameters must be a vector of drift and volatility.}
#' @return data.frame of solution
#' @export optimalGBM
optimalGBM <- function(bankroll, t, spot, rate, parameters)
{

  mu <- parameters[1]
  volat <- parameters[2]
  if(volat <= 0)
  {
    stop("volatility must be positive")
  }
  dynamics <- list(function(t, s) mu,
                   function(t, s) volat
  )
  z <- optimalItoProcess(bankroll, t, spot, rate, dynamics)
  return(z)
}

#' Kelly portfolio in continuous time for GBM market
#'
#' @param drift the drift vector of the GBMs
#' @param Sigma the covariance matrix of the GBMs
#' @param rate the rate earned on cash in a money market account, etc
#' @param restraint percentage of wealth to bankroll
#' @param direction the direction of the bet: long or short
#'
#' @description {A wrapper to \code{solve.QP} for computing optimal portfolios
#' under a market of correlated geometric Brownian motions.}
#' @details {The drift of the log-returns and covariance must be passed.}
#' @return list of growth rate and optimal portfolio
#'
#' @importFrom quadprog solve.QP
#' @export kellyPortfolioGBM
kellyPortfolioGBM <- function(drift, Sigma, rate = 0, restraint = 1, direction = "long")
{
  # Get number of assets of portfolio
  N <- length(drift)

  # TODO: Add a semi-positive definitive check here!...
  # Constraint vector and matrix, want <= 1 for budget constraint, and >= 0 for no-short selling
  if(direction == "long")
  {
    # -restraint for budget constraint (to get u 1^T <= restraint), >0, -u>-1 for 0<u<1
    bvec <- c(-restraint, rep(0, N), rep(-1, N))
    # for row is 1s for budget equality constraint, then diag matrices
    Amat <- cbind(matrix(rep(-1, N), nrow = N), diag(x = 1, N, N), diag(x = -1, N, N))
  } else if(direction == "short")
  {
    bvec <- c(-restraint, rep(0, N), rep(-1, N))
    Amat <- cbind(matrix(rep(1, N), nrow = N), diag(x = -1, N, N), diag(x = 1, N, N))
  }

  # Quadratic programming routine
  qp <- quadprog::solve.QP(Dmat = Sigma, dvec = drift-rate, Amat = Amat, bvec)
  optimalWeight <- as.matrix(round(qp$solution, 8))

  growth <- rate+t(drift-rate)%*%(optimalWeight)-0.5*t(optimalWeight)%*%Sigma%*%optimalWeight
  # Append cash weight
  bet <- as.matrix(c(optimalWeight, 1-sum(optimalWeight)))
  rownames(bet) <- c(colnames(Sigma), "cash")

  return(bet)
}

#' Maximum log-growth rate for continuous time GBM portfolios
#'
#' @param drift the drift vector of the GBMs
#' @param Sigma the covariance matrix of the GBMs
#' @param rate the rate earned on cash in a money market account, etc
#' @param restraint percentage of wealth to bankroll
#' @param direction the direction of the bet: long or short
#'
#' @description {A wrapper to \code{solve.QP} for computing optimal portfolios
#' under a market of correlated geometric Brownian motions.}
#' @details {The drift of the log-returns and covariance must be passed.}
#'
#' @return list of growth rate and optimal portfolio
#' @importFrom quadprog solve.QP
#' @export entropyPortfolioGBM
entropyPortfolioGBM <- function(drift, Sigma, rate = 0, restraint = 1, direction = "long")
{
  optimalWeight <- kellyPortfolioGBM(drift, Sigma, rate, restraint, direction)
  # Remove cash holding
  optimalWeight <- optimalWeight[-nrow(optimalWeight),]
  # Compute maximum growth rate
  growth <- rate+t(drift-rate)%*%(optimalWeight)-0.5*t(optimalWeight)%*%Sigma%*%optimalWeight
  return(growth)
}
