#' Jumps integral in Kelly criterion for Merton dynamics
#'
#' @param a the control value
#' @param g the function to compute expectation of \eqn{J = e^Y-1}
#' @param alpha the mean jump size
#' @param beta the standard deviation of the jump size
#' @param subdivisions number of subdivisions in the numerical integration
#'
#' @description {Compute various integrals against a density function for the
#' log jump size.}
#' @return numeric
#' @export jint
jint <- function(a, g, alpha, beta, subdivisions = 200)
{
  d <- matrix(0, nrow = length(a))
  # Truncate interval
  lb <- stats::qnorm(0.99, alpha, beta, FALSE)
  ub <- stats::qnorm(0.99, alpha, beta)

  # DIY Composite Trapezoid scheme
  y <- seq(lb, ub, length.out = subdivisions+1)
  h <- diff(y)[1]

  # Function for R's base integrate function
  # integrand <- function(y, aa)
  # {
  #   g(exp(y)-1, aa)*stats::dnorm(y, alpha, beta)
  # }
  for(i in 1:length(a))
  {
    # R's base integrate function
    # d[i] <- stats::integrate(integrand, lb, ub, aa = a[i], subdivisions = subdivisions)$value

    # Trapezoid scheme requires vectors of the function at the nodes
    integrand <- g(exp(y)-1, a[i])*stats::dnorm(y, alpha, beta)
    d[i] <- (integrand%*%c(1, rep(2, length(y)-2), 1))*h/2
  }
  return(d)
}

#' Kelly-criterion under Merton's jump diffusion
#'
#' @param param vector of parameters defining the jump-diffusion dynamics. See details
#' @param rate the risk-free rate, or money-market account interest rate
#' @param iterations number of iterations to use in the Newton-Raphson method for finding the optimal control as a root
#' @param subdivisions number of subdivisions to use in numerical integrations
#' @param tol tolerance for the zero to be found.
#'
#' @description {Compute the optimal allocation fraction under Merton's jump diffusion dynamics by
#' root finding via Newton-Raphson's method. }
#' @return list
#' @export kellyMerton
kellyMerton <- function(param, rate = 0, iterations = 500, subdivisions = 500, tol = 10^-6)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  alpha <- param[4]
  beta <- param[5]

  g1 <- function(x, a) (a*x^2)/(1+a*x)
  g <- function(x, a) (x/(1+a*x))^2
  a <- matrix(0, nrow = iterations)
  a[1] <- (mu-rate)/(volat^2)
  for(i in 1:(iterations-1))
  {
    Li <- mu-rate-a[i]*volat^2-lambda*jint(a[i], g1, alpha, beta, subdivisions)
    Lip <- -volat-lambda*jint(a[i], g, alpha, beta, subdivisions)
    a[i+1] <- a[i]-Li/Lip
  }
  rootCheck <- mu-rate-a[iterations]*volat^2-lambda*jint(a[iterations], g1, alpha, beta, subdivisions)
  msg <- ""
  if(abs(rootCheck) < tol)
  {
    msg <- "converged"
  } else
  {
    msg <- "not within tolerance"
  }
  output <- list(root = a[iterations], value = rootCheck, msg = msg)
  return(output)
}

#' Kelly-criterion entropy under Merton's jump diffusion
#'
#' @param param vector of parameters defining the jump-diffusion dynamics. See details
#' @param rate the risk-free rate, or money-market account interest rate
#' @param iterations number of iterations to use in the Newton-Raphson method for finding the optimal control as a root
#' @param subdivisions number of subdivisions to use in numerical integrations
#' @param tol tolerance for the zero to be found.
#'
#' @description {Compute the entropy, i.e. growth rate of the Kelly allocation
#' under Merton's jump diffusion. }
#' @return list
#' @export entropyMerton
entropyMerton <- function(param, rate = 0, iterations = 500, subdivisions = 500, tol = 10^-6)
{
  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  alpha <- param[4]
  beta <- param[5]
  a <- kellyMerton(param, rate, iterations, subdivisions, tol)
  a <- a$root
  g <- function(x, aa) log(1+aa*x)-aa*x
  v <- rate+(mu-rate)*a-0.5*(volat*a)^2+lambda*jint(a, g, alpha, beta, subdivisions)
  return(v)
}

#' Euler-Maruyama scheme for Kelly portfolio under Merton's jump diffusion
#'
#' @param t maturity to simulate under
#' @param param vector of parameters defining the jump-diffusion dynamics. See details
#' @param rate the risk-free rate, or money-market account interest rate
#' @param n number of time-subintervals to use in the Euler-Maruyama scheme
#'
#' @description {Simulate a sample path of \eqn{\log X_t} of the Kelly-portfolio
#' under Merton's jump diffusion dynamics, via the Euler-Maruyama scheme}
#' @return list
#' @export optimalMertonPath
optimalMertonPath <- function(t, param, rate = 0, n = 1000)
{

  mu <- param[1]
  volat <- param[2]
  lambda <- param[3]
  alpha <- param[4]
  beta <- param[5]
  eta <- exp(alpha+0.5*beta^2)-1
  a <- kellyfractions::kellyMerton(param, rate)$root
  # Solution vector to populate
  x <- matrix(0, nrow = n+1)
  x[1] <- 0
  # Solution vector to populate
  s <- matrix(0, nrow = n+1)
  s[1] <- 0
  # Fixed time-step size
  k <- t/n
  # Time grid
  tt <- seq(0, t, k)
  # Loop through time grid
  for(i in 2:(n+1))
  {

    # Simulate number of jumps in single time-step
    m <- stats::rpois(1, lambda = lambda*k)
    # Compound Poisson sum of the random jump-amplitudes
    y <- stats::rnorm(m, alpha, beta)
    cps_pois <- sum(log(a*(exp(y)-1)+1))
    s[i] <- s[i-1] + (mu-0.5*(volat)^2-eta*lambda)*k+volat*stats::rnorm(1, 0, sd = sqrt(k))+sum(y)
    x[i] <- x[i-1] + ((mu-rate)*a+rate-0.5*(volat*a)^2-a*eta*lambda)*k+a*volat*stats::rnorm(1, 0, sd = sqrt(k))+cps_pois

  }
  X <- data.frame(t = tt, S = s, X = x)
  return(X)
}


# These are verifying the computation of these jump integrals E(g(J))
# monteCarloJint <- function(a, g, alpha, beta, n)
# {
#
#   d <- matrix(0, nrow = length(a))
#   for(i in 1:length(a))
#   {
#     y <- rnorm(n, alpha, beta)
#     j <- g(exp(y)-1, a[i])
#     d[i] <- mean(j)
#   }
#   return(d)
# }
#
# trapJint <- function(a, g, alpha, beta, n)
# {
#   d <- matrix(0, nrow = length(a))
#   for(i in 1:length(a))
#   {
#     # Truncate interval
#     lb <- qnorm(0.99, alpha, beta, FALSE)
#     ub <- qnorm(0.99, alpha, beta)
#     y <- seq(lb, ub, length.out = n+1)
#     h <- diff(y)[1]
#     integrand <- g(exp(y)-1, a[i])*dnorm(y, alpha, beta)
#     d[i] <- (integrand%*%c(1, rep(2, length(y)-2), 1))*h/2
#   }
#   return(d)
# }






