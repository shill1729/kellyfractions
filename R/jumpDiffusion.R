#' Jump integral term for log utility
#'
#' @param a control variate
#' @param distr distribution of the jump-sizes
#' @param jump_param parameters for distribution a list
#' @param nstd number of standard deviations for integral region
#'
#' @return numeric
jump_value_integral <- function(a, distr = "norm", jump_param, nstd = 4)
{
  if(distr == "norm")
  {
    jm <- jump_param$mean
    jv <- jump_param$sd
    lb <- jm-nstd*jv
    ub <- jm+nstd*jv
  } else if(distr == "unif")
  {
    lb <- jump_param$min
    ub <- jump_param$max
  } else{
    stop("TO DO: implement kou bounds for integral region")
  }
  # Get the PDF of Y = log J
  dist <- paste("d", distr, sep = "")
  f_Y<- function(y) do.call(what = dist, args = list(y, unlist(jump_param)))
  integrand <- function(y, a) {
    (log(1+a*(exp(y)-1))-a*(exp(y)-1))*f_Y(y)
  }
  stats::integrate(f = integrand, lower = lb, upper = ub, subdivisions = 200, a = a)$value
}

#' Jump entropy term
#'
#' @param a fraction input
#' @param m drift rate
#' @param r bond rate
#' @param v volatility
#' @param lambda mean number of jumps per year
#' @param distr distribution of jump sizes
#' @param jump_param parameters for distribution
#'
#' @return numeric
jump_entropy_term <- function(a, m, r, v, lambda, distr, jump_param)
{

  (m-r)*a+r-0.5*(v^2)*a^2+lambda*jump_value_integral(a, distr, jump_param)
}

#' Entropy for jump diffusion log optimal control
#'
#' @param a control fraction
#' @param t time at start
#' @param x bankroll
#' @param mu drift rate
#' @param rate bond rate
#' @param volat volatility
#' @param lambda mean number of jumps per year
#' @param distr distribution of jump size
#' @param jump_param parameters of above distribution
#'
#' @return numeric
#' @export entropy_jdf
entropy_jdf <- function(a, t, x, mu, rate, volat, lambda, distr, jump_param)
{
  log(x)+t*jump_entropy_term(a, mu, rate, volat, lambda, distr, jump_param)
}

#' Jump integral term for log utility
#'
#' @param a control variate
#' @param distr distribution of the jump-sizes
#' @param jump_param parameters for distribution a list
#' @param nstd number of standard deviations for integral region
#'
#' @return numeric
jump_integral <- function(a, distr = "norm", jump_param, nstd = 4)
{
  if(distr == "norm")
  {
    jm <- jump_param$mean
    jv <- jump_param$sd
    lb <- jm-nstd*jv
    ub <- jm+nstd*jv
  } else if(distr == "unif")
  {
    lb <- jump_param$min
    ub <- jump_param$max
  } else{
    stop("TO DO: implement kou bounds for integral region")
  }
  # Get the PDF of Y = log J
  dist <- paste("d", distr, sep = "")
  f_Y<- function(y) do.call(what = dist, args = list(y, unlist(jump_param)))
  integrand <- function(y, a) {
    ((a*(exp(y)-1)^2)/(1+a*(exp(y)-1)))*f_Y(y)
  }
  stats::integrate(f = integrand, lower = lb, upper = ub, subdivisions = 200, a = a)$value
}

#' Function to compute fixed point of for log optimal criterion under jump-diffusions
#'
#' @param a control variate
#' @param mu drift
#' @param rate discounting rate
#' @param volat volatility
#' @param lambda mean rate of jumps
#' @param distr distribution of jump sizes
#' @param jump_param parameters of distribution
#'
#' @return numeric
#' @export kjdf_fixed_point
kjdf_fixed_point <- function(a, mu, rate, volat, lambda, distr, jump_param)
{
  (mu-rate)/(volat^2)-(lambda/(volat^2))*jump_integral(a, distr, jump_param)
}

#' Function to find root of for log optimal criterion under jump-diffusions
#'
#' @param a control variate
#' @param mu drift
#' @param rate discounting rate
#' @param volat volatility
#' @param lambda mean rate of jumps
#' @param distr distribution of jump sizes
#' @param jump_param parameters of distribution
#'
#' @return numeric
#' @export kjdf_root
kjdf_root <- function(a, mu, rate, volat, lambda, distr, jump_param)
{
  mu-rate-a*volat^2-lambda*jump_integral(a, distr, jump_param)
}

#' Moment-generating function of a uniform RV
#'
#' @param t mgf argument
#' @param min minimum of interval
#' @param max maximum of interval
#'
#' @return numeric
mgfunif <- function(t, min, max)
{
  if(t == 0)
  {
    return(1)
  } else{
    return((exp(t*max)-exp(t*min))/(t*(max-min)))
  }

}

#' Moment-generating function of a Gaussian RV
#'
#' @param t mgf argument
#' @param mean mean of the distribution
#' @param sd standard deviation
#'
#' @return numeric
mgfnorm <- function(t, mean, sd)
{
  exp(t*mean+0.5*(t^2)*(sd^2))
}

#' Moment generating function for Displaced Kou distribution
#'
#' @param t mgf argument
#' @param prob probability of mixtures
#' @param alpha mean jump size up
#' @param beta mean jump size down
#'
#' @return numeric
mgfkou <- function(t, prob, alpha, beta)
{
  mgf_dkou(t, prob, alpha, beta, ku = 0, kd = 0)
}

#' Moment generating function for Displaced Kou distribution
#'
#' @param t mgf argument
#' @param prob probability of mixtures
#' @param alpha mean jump size up
#' @param beta mean jump size down
#' @param ku displacement up
#' @param kd displacement down
#'
#' @return numeric
mgfdkou <- function(t, prob, alpha, beta, ku, kd)
{
  mgf_dkou(t, prob, alpha, beta, ku, kd)
}

#' log optimal criterion under jump-diffusions via fixed-point iterations
#'
#' @param mu drift
#' @param rate discounting rate
#' @param volat volatility
#' @param lambda mean rate of jumps
#' @param distr distribution of jump sizes
#' @param jump_param parameters of distribution
#' @param iterations number of iterations in fixed-point scheme
#'
#' @return numeric
#' @export kellyJumpDiffusionFP
kellyJumpDiffusionFP <- function(mu, rate, volat, lambda, distr, jump_param, iterations = 100)
{
  # Get the PDF of Y = log J
  mgf <- paste("mgf", distr, sep = "")
  M_Y <- function(t) do.call(what = mgf, args = as.list(c(t, unlist(jump_param))))
  bound <- M_Y(2)-2*M_Y(1)+1

  if(bound > (volat^2)/lambda)
  {
    print(c(bound,volat^2/lambda))
    warning("Criterion fails")

  }
  fp <- matrix(0, nrow = iterations)
  fp[1] <- kjdf_fixed_point(0.5, mu, rate, volat, lambda, distr, jump_param)
  for(i in 2:iterations)
  {
    fp[i] <- kjdf_fixed_point(fp[i-1], mu, rate, volat, lambda, distr, jump_param)
  }
  graphics::plot(fp, type = "l")
  return(fp[iterations])
}

#' log optimal criterion under jump-diffusions
#'
#' @param mu drift
#' @param rate discounting rate
#' @param volat volatility
#' @param lambda mean rate of jumps
#' @param distr distribution of jump sizes
#' @param jump_param parameters of distribution
#'
#' @return numeric
#' @export kellyJumpDiffusion
kellyJumpDiffusion <- function(mu, rate, volat, lambda, distr, jump_param)
{
  if(mu-rate<=0)
  {
    stop("mu must be greater than rate")
  }
  if(distr == "norm")
  {

    eta1 <- exp(jump_param$mean+0.5*jump_param$sd^2)-1
    eta2 <- exp(-jump_param$mean+0.5*jump_param$sd^2)-1
    mgfs <- eta1+eta2
  } else if(distr == "unif")
  {
    max <- jump_param$max
    min <- jump_param$min

    eta1 <- (exp(max)-exp(min))/(max-min)-1
    eta2 <- -(exp(-max)-exp(-min))/(max-min)-1
    mgfs <- eta1+eta2
  } else if(distr == "kou")
  {
    prob <- jump_param$prob
    alpha <- jump_param$alpha
    beta <- jump_param$beta

    eta1 <- mgfkou(1, prob, alpha, beta)-1
    eta2 <- mgfkou(-1, prob, alpha, beta)-1
    mgfs <- eta1+eta2
  } else if(distr == "dkou")
  {
    prob <- jump_param$prob
    alpha <- jump_param$alpha
    beta <- jump_param$beta
    ku <- jump_param$ku
    kd <- jump_param$kd
    eta1 <- mgfdkou(1, prob, alpha, beta, ku, kd)-1
    eta2 <- mgfdkou(-1, prob, alpha, beta, ku, kd)-1
    mgfs <- eta1+eta2
  }
  if(mu-rate > volat^2+lambda*mgfs)
  {
    print("Excess return")
    print(mu-rate)
    print("Variance+jump risk")
    print(volat^2+lambda*mgfs)
    warning("mu-rate must be less than variance plus mean jumps times sum of mgf at one and negative one")
  }
  stats::uniroot(kjdf_root, interval = c(0,1), mu = mu, rate = rate, volat = volat, lambda = lambda, distr = distr, jump_param)$root
}
