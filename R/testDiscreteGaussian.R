# library(kellyfractions)
# bankroll <- 100
# distr <- "norm"
# param <- c(0.005, 0.08)
# rate <- 0
# spot <- 100
# n <- 1000
# input <- data.frame(bankroll, distr, mean = param[1], sd = param[2], rate)
# w <- kellyDTFM(distr, param, rate)
# g <- entropyDTFM(distr, param, rate)
# s <- optimalDTFM(n, spot, bankroll, distr, param, rate)
# output <- data.frame(fraction = w, growth_rate = g)
# print(list(intput = input, output = output))
