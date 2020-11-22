#' Kelly-criterion for discrete-time financial market models
#'
#' @param distr distribution name for the (daily) arithmetic returns
#' @param param the parameters of the distribution
#' @param rate the risk-neutral rate of the money-market account
#' @param h step size to avoid numeric blow-ups
#'
#' @description {Approximates the integral in the log utility problem and
#' then calls uniroot or Newton-Raphson solver}
#' @return numeric
#' @importFrom findistr dgmm
#' @importFrom findistr dstable
#' @export kellyDTFM
kellyDTFM <- function(distr, param, rate = 0, h = 10^-4)
{
  ddistr <- paste("d", distr, sep = "")
  ddistr <- get(ddistr)

  if(distr == "unif")
  {
    return(kellyUniform(param, rate))
  }
  if(distr == "norm")
  {

    if(param[1] < rate)
    {
      return(-0.5)
    }
    if(param[1]-rate > param[2]^2)
    {
      # warning("Optimal point is leveraged > 1, using probabilistic bounds for range")
      mu <- param[1]
      v <- param[2]
      K <- stats::qnorm(1-h/2, mean = 0, sd = 1)
      LB <- mu-K*v
      UB <- mu+K*v
      region <- c(LB, UB)
      admissible <- c((1+rate)/(rate-UB)+h, (1+rate)/(rate-LB)-h)
    } else
    {
      # print("Using (-1, Inf) for returns support")
      region <- c(-1+h, Inf)
      admissible <- c(h, 1-h)
    }
  }
  if(distr == "gmm")
  {
    mu <- param[1, ]%*%param[2, ]
    v <- sqrt(param[1, ]%*%param[3, ]^2+param[1, ]%*%param[2, ]^2-(param[1, ]%*%param[2, ])^2)
    K <- stats::qnorm(1-h/2, mean = 0, sd = 1)
    LB <- mu-K*v
    UB <- mu+K*v
    region <- c(LB, UB)
    admissible <- c((1+rate)/(rate-UB)+h, (1+rate)/(rate-LB)-h)
    # region <- c(-1+h, Inf)
    # admissible <- c(h, 1-h)
  }
  if(distr == "stable")
  {
    UB <- stats::uniroot(function(x) libstableR::stable_cdf(x, param)-libstableR::stable_cdf(-x, param)-0.99, interval = c(-1+h, 1))
    UB <- UB$root
    LB <- -UB
    region <- c(LB, UB)
    admissible <- c((1+rate)/(rate-UB)+h, (1+rate)/(rate-LB)-h)
  }


  integrand <- function(y, x)
  {
    if(distr == "unif" || distr == "norm")
    {
      args <- c(list(y), as.list(param))
    } else if(distr == "gmm")
    {
      args <- list(c(y), param[1, ], param[2, ], param[3, ])
    } else if(distr == "stable")
    {
      args <- list(c(y), param)
    }
    ((y-rate)/(1+rate+x*(y-rate)))*do.call(ddistr, args)
  }
  gp <- function(x)
  {
    stats::integrate(f = integrand, lower = region[1], upper = region[2], x = x)$value
  }
  optimal <- stats::uniroot(gp, interval = admissible)$root
  return(optimal)
}


#' Log-growth n for discrete-time financial market models
#'
#' @param distr distribution name for the (daily) arithmetic returns
#' @param param the parameters of the distribution
#' @param rate the risk-neutral rate of the money-market account
#' @param h step size for numeric blow-ups
#'
#' @description {Approximates the integral in the log utility problem. This
#' is a new interface for \code{entropy_integral}, which will be deprecated and removed.}
#'
#' @return numeric
#' @export entropyDTFM
entropyDTFM <- function(distr, param, rate = 0, h = 10^-4)
{
  ddistr <- paste("d", distr, sep = "")
  ddistr <- get(ddistr)
  xopt <- kellyDTFM(distr, param, rate, h)
  if(distr == "unif")
  {
    region <- param
  }
  if(distr == "norm")
  {

    if(param[1] < rate)
    {
      return(0.01)
    }
    if(param[1]-rate > param[2]^2)
    {
      # warning("Optimal point is leveraged > 1, using probabilistic bounds for range")
      mu <- param[1]
      v <- param[2]
      K <- stats::qnorm(1-h/2, mean = 0, sd = 1)
      LB <- mu-K*v
      UB <- mu+K*v
      region <- c(LB, UB)
    } else
    {
      # print("Using (-1, Inf) for returns support")
      region <- c(-1+h, Inf)
    }

  }
  if(distr == "gmm")
  {
    mu <- param[1, ]%*%param[2, ]
    v <- sqrt(param[1, ]%*%param[3, ]^2+param[1, ]%*%param[2, ]^2-(param[1, ]%*%param[2, ])^2)
    K <- stats::qnorm(1-h/2, mean = 0, sd = 1)
    LB <- mu-K*v
    UB <- mu+K*v
    region <- c(LB, UB)
    # region <- c(-1+h, Inf)
    # admissible <- c(h, 1-h)
  }
  if(distr == "stable")
  {
    BD <- libstableR::stable_q(0.99, param)
    LB <- -BD
    region <- c(-BD, BD)
  }
  integrand <- function(y, x)
  {
    if(distr == "unif" || distr == "norm")
    {
      args <- c(list(y), as.list(param))
    } else if(distr == "gmm")
    {
      args <- list(c(y), param[1, ], param[2, ], param[3, ])
    } else if(distr == "stable")
    {
      args <- list(c(y), param)
    }

    (log(1+rate+x*(y-rate)))*do.call(ddistr, args)
  }
  g <- stats::integrate(integrand, lower = region[1], upper = region[2], x = xopt)$value
  return(g)

}
