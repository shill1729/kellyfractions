
#' Kelly fraction for tandemn moneyline bets
#'
#' @param p true probability of favorite winning
#' @param a winnings if favorite wins
#' @param b wager to throw down for favorite \eqn{-b}
#' @param u winnings if underdog wins \eqn{+u}
#' @param v wager to throw down for underdog
#'
#' @description {A Kelly-criterion for a bet on two wagers for favorite
#' vs underdog team to win.}
#' @return numeric
#' @export kellyTandemMoneyline
kellyTandemMoneyline <- function(p, a, b, u, v)
{
  region_check <- (v-1)*(u+b)/((1+u)*(a+v))
  if(region_check >= 1)
  {
    stop("Feasible region DNE")
  }
  x <- ((1+u)*(a+v)*p-(1-v)*(b+u)*(1-p))/((b+u)*(a+v))
  return(x)
}

#' Entropy-growth for log-optimal tandem moneyline bets
#'
#' @param p true probability of favorite winning
#' @param a winnings if favorite wins
#' @param b wager to throw down for favorite \eqn{-b}
#' @param u winnings if underdog wins \eqn{+u}
#' @param v wager to throw down for underdog
#'
#' @description {The entropy-growth for the log-optimal bet on two wagers for favorite
#' vs underdog team to win.}
#' @return numeric
#' @export entropyTandemMoneyline
entropyTandemMoneyline <- function(p, a, b, u, v)
{
  x <- kellyTandemMoneyline(p, a, b, u, v)
  p*log(1+(a+v)*x-v)+(1-p)*log(1+u-(b+u)*x)
}


#' Simulate a series of IID trials of moneyline under log-optimal growth
#'
#' @param bankroll initial bankroll to bet
#' @param p true chance of outcome
#' @param wagers the vector of wagers, see details
#' @param trials the number of trials to simulate
#'
#' @description {Simulate a series of moneyline trials under log-optimal growth.}
#' @details {The argument \code{wagers} must contain four wagers, the first two the odds
#' for the favorite, the latter two the odds for the underdog.}
#' @return vector
#' @export simulateMoneyline
simulateMoneyline <- function(bankroll, p, wagers, trials = 100)
{
  a <- wagers[1]
  b <- wagers[2]
  u <- wagers[3]
  v <- wagers[4]

  # Initial wealth process
  z <- matrix(data = 0, nrow = trials+1)
  z[1] <- bankroll
  # Optimal strategy
  x <- kellyTandemMoneyline(p, a, b, u, v)
  y <- 1-x

  # Simulate the Bernoulli outcomes 0,1 (loss, win)
  outcomes <- stats::rbinom(n = trials, size = 1, prob = p)
  # Convert to bet outcomes:
  epsilon <- a*outcomes-b*(1-outcomes) # Favorite wins
  delta <- u*(1-outcomes)-v*outcomes # Underdog wins
  for(i in 1:(trials))
  {
    z[i+1] <- z[i]+epsilon[i]*z[i]*x+delta[i]*z[i]*y
  }
  return(z)
}
