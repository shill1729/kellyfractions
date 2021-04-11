# library(kellyfractions)
# bankroll <- 100
# distr <- "unif"
# param <- c(-0.02, 0.0201)
# rate <- 0
# spot <- 100
# n <- 1000
# input <- data.frame(bankroll, distr, min = param[1], max = param[2], rate)
# w <- kellyDTFM(distr, param, rate)
# g <- entropyDTFM(distr, param, rate)
# s <- optimalDTFM(n, spot, bankroll, distr, param, rate)
# output <- data.frame(fraction = w, growth_rate = g)
# print(list(intput = input, output = output))
#
