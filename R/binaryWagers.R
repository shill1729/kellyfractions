#' The optimal bet for a binary wager
#'
#' @param p true chance of the event being bet on
#' @param a dollars to win
#' @param b dollars to risk
#'
#' @description {The classical Kelly fraction that optimizes long-term log-growth.}
#' @return numeric
#' @export kellyBinary
kellyBinary <- function(p, a, b)
{
  return(p/b-(1-p)/a)
}



#' Growth rate for Kelly fraction for arbitrary binary (non-random) payoffs
#'
#' @param p true chance of the event being bet on
#' @param a dollars to win
#' @param b dollars to risk
#'
#' @description {The growth rate or entropy of the Kelly strategy}
#' @return numeric
#' @export entropyBinary
entropyBinary <- function(p, a, b)
{
  p*log(p)+(1-p)*log(1-p)+log(a+b)+p*log(a/b)-log(a)
}



#' Simulate Kelly strategy for generalized coin toss
#'
#' @param bankroll initial bankroll to bet with
#' @param p probability of winning each trial
#' @param a payout on winning toss on top of wager
#' @param b wager lost on each trial on losing toss
#' @param trials number of trials to simulate
#'
#' @description {Simulating the generalized coin tossing gambling game.}
#' @details {Start off with a given value of wealth and bet the Kelly fraction each round.}
#' @return vector
#' @export optimalBinary
optimalBinary <- function(bankroll, p, a, b, trials = 100)
{
  # Initial wealth process
  z <- matrix(data = 0, nrow = trials+1)
  z[1] <- bankroll
  y <- kellyBinary(p, a, b)

  # Simulate the Bernoulli outcomes 0,1 (loss, win)
  outcomes <- stats::rbinom(n = trials, size = 1, prob = p)
  # Convert to bet outcomes:
  epsilon <- a*outcomes-b*(1-outcomes)
  for(i in 1:(trials))
  {
    z[i+1] <- z[i]+epsilon[i]*z[i]*y
  }
  return(z)
}
