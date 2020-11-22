
# kellyfractions

<!-- badges: start -->
<!-- badges: end -->

The goal of kellyfractions is to ...

## Installation

You can install the released version of kellyfractions from [CRAN](https://CRAN.R-project.org) with:

``` r
devtools::install_github("shill1729/kellyfractions")
```

## Examples of optimal sample paths
How to do simulations:
### Binary wagers

``` r
library(kellyfractions)
# +/- $1 bets at 0.6 chance
samplepath <- optimalBinary(bankroll = 1000, p = 0.6, a = 1, b = 1)
plot(samplepath, type = "l", main = "Binary wagers")
```

## Testing solution accuarcy
Here we compare three methods that should coincide for the special case of
cev=1.
```r
# library(kellyfractions)
# mu <- 0.10
# volat <- 0.25
# rate <- 0
# cev <- 1
# dynamics <- list(function(t, s) mu,
#                  function(t, s) volat*s^(cev-1)
# )
# kito <- kellyItoProcess(0, s = 100, dynamics, rate)
# kgbm <- kellyGBM(mu, volat, rate)
# kcev <- kellyCEV(spot = 100, rate, parameters = c(mu, volat, cev))
# gito <- entropyItoProcess(t = 1, s = 100, dynamics = dynamics, rate)
# ggbm <- entropyGBM(mu, volat, rate)
# gcev <- entropyCEV(t = 1, spot = 100, rate = rate, parameters = c(mu, volat, cev))
# # These should be approximately equal
# output <- rbind(c(kgbm, kito, kcev),
#                 c(ggbm, gito, gcev)
#                 )
# colnames(output) <- c("gbm", "ito", "cev")
# rownames(output) <- c("fraction", "growth")
# print(output)
```
