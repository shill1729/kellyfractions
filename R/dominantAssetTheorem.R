#' Estimate the expectation E((1+R_j)/(1+Ri)) for every pair of assets
#'
#' @param r the periodic arithmetic returns of the assets
#' @description {This is simply a sample average of the ratio of returns
#' between each asset in the data-set.}
#' @details {If the row entries are all less than one (except for the diagonal), then
#' this is the dominant asset.}
#'
#' @return matrix
#' @export dominant_asset_matrix
dominant_asset_matrix <- function(r)
{
  n <- ncol(r)
  # Estimating the dominant asset matrix
  M <- matrix(0, n, n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      M[i,j] = mean((1+r[,j])/(1+r[,i]))
    }
  }
  colnames(M) <- colnames(r)
  rownames(M) <- colnames(r)
  return(M)
}

#' Find the dominant asset via the Kelly-criterion
#'
#' @param r the periodic arithmetic returns of the assets
#' @description {Returns in the column index of the dominant asset in a
#' given dataset of periodic arithemetic returns.}
#' @details {This will reutrn an error if no dominant asset is found.}
#'
#' @return matrix
#' @export dat_strategy
dat_strategy <- function(r)
{
  M <- dominant_asset_matrix(r)
  dom_index <- which(apply(M<=1, 1, all))
  if(length(dom_index) == 0)
  {
    warning("No single dominant asset")
    approxDom <- M-1
    approxDom <- which.min(apply(approxDom, 1, sum))
    return(approxDom)
  }
  return(dom_index)
}
