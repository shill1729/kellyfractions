# library(kellyfractions)
# bankroll <- 100
# distr <- "stable"
# param <- c(1.53, 0.20, 0.024, -0.0001)
# rate <- 0
# spot <- 100
# n <- 30
# # set.seed(17)
# print(param)
# input <- data.frame(bankroll, distr, rate)
# w <- kellyDTFM(distr, param, rate)
# g <- entropyDTFM(distr, param, rate)
# s <- optimalDTFM(n, spot, bankroll, distr, param, rate)
# output <- data.frame(fraction = w, growth_rate = g)
# print(list(intput = input, output = output))
