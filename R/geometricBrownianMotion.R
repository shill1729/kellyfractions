#' The Kelly-fraction for a stock following GBM dynamics
#'
#' @param drift the mean drift rate of the GBM
#' @param rate the risk-free rate of the bond
#' @param volat the volatility of the GBM
#'
#' @description {The famous Kelly-fraction, the amount to invest
#' in a risky stock following a geometric Brownian motion.}
#' @return numeric
#' @export kellyGBM
kellyGBM <- function(drift, rate, volat)
{
  f <- (drift-rate)/volat^2
  lambda <- (drift-rate)/volat
  g <- rate+0.5*lambda^2
  return(c(fraction = f, growth = g))
}

#' Kelly portfolio in continuous time for GBM market
#'
#' @param log_ret the data-set of (daily) log-returns. Assumes it is a matrix with colnames of tickers
#' @param rate the rate earned on cash in a money market account, etc
#' @param restraint percentage of wealth to bankroll
#' @param direction the direction of the bet: long or short
#' @param sample_mean the mean vector, optional
#' @param sample_cov the covariance matrix, optional
#'
#' @description {For continuous time market of GBMs, the log optimal portfolio is given by a quadratic
#' programming problem.}
#' @details {The function estimates the entire sample set for mean returns and covariance and then runs the optimization routine by calling solve.qp from quadprog package.}
#' @return list of growth rate and optimal portfolio
#' @import quadprog
#' @export portfolioGBM
portfolioGBM <- function(log_ret, rate = 0, restraint = 1, direction = "long", sample_mean = NULL, sample_cov = NULL)
{
  # Get number of assets of portfolio
  N <- dim(log_ret)[2]
  # Compute sample mean and covariance matrix
  if(is.null(sample_mean) || is.null(sample_cov))
  {
    sample_cov <- stats::cov(log_ret)*252
    sample_mean <- apply(log_ret, 2, mean)*252+0.5*diag(sample_cov)

  }
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
  qp <- quadprog::solve.QP(Dmat = sample_cov, dvec = sample_mean-rate, Amat = Amat, bvec)
  optimalWeight <- as.matrix(round(qp$solution, 8))
  growth <- rate+t(sample_mean-rate)%*%(optimalWeight)-0.5*t(optimalWeight)%*%sample_cov%*%optimalWeight
  # append cash weight
  bet <- as.matrix(c(optimalWeight, 1-sum(optimalWeight)))
  rownames(bet) <- c(colnames(log_ret), "cash")

  return(list(growth = growth, bet = bet))
}
