#' The optimal bet for a binary wager
#'
#' @param chance true chance of event
#' @param win dollars to win
#' @param risk dollars to risk
#'
#' @description {The classical Kelly fraction that optimizes long-term log-growth.}
#' @return numeric
#' @export optimalBinary
optimalBinary <- function(chance, win, risk)
{
  return(chance/risk-(1-chance)/win)
}


