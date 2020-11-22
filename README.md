
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

