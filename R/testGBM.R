# library(kellyfractions)
# bankroll <- 100
# drift <- 0.1
# volat <- 0.2
# rate <- 0.09
# spot <- 10
# maturity <- 20
# input <- data.frame(bankroll, maturity, spot, drift, rate, volat)
# w <- kellyGBM(drift, volat, rate)
# g <- entropyGBM(drift, volat, rate)
# s <- optimalGBM(bankroll, maturity, spot, rate, c(drift, volat))
# output <- data.frame(fraction = w, growth_rate = g)
# print(list(intput = input, output = output))
