# library(kellyfractions)
# bankroll <- 100
# distr <- "gmm"
# param <- rbind(c(0.1, 0.9),
#                c(0.01, 0),
#                c(0.07, 0.03)
#                )
# rownames(param) <- c("probs", "means", "sds")
# rate <- 0
# spot <- 100
# n <- 1000
# print(param)
# input <- data.frame(bankroll, distr, rate)
# w <- kellyDTFM(distr, param, rate)
# g <- entropyDTFM(distr, param, rate)
# s <- optimalDTFM(n, spot, bankroll, distr, param, rate)
# output <- data.frame(fraction = w, growth_rate = g)
# print(list(intput = input, output = output))
