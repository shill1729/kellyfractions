# library(kellyfractions)
# bankroll <- 100
# p <- 0.51
# a <- 1
# b <- 1
# trials <- 500
# input <- data.frame(bankroll, p, a, b, trials)
# w <- kellyBinary(p, a, b)
# g <- entropyBinary(p, a, b)
# s <- optimalBinary(bankroll, p, a, b, trials)
# output <- data.frame(fraction = w, growth_rate = g)
# # Plot sample path and print optimal bet and growth-rate
# plot(s, type = "s", main = "Sample path", xlab = "trials", ylab = "wealth")
# print(list(intput = input, output = output))
