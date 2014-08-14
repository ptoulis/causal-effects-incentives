## Copyright 2013 Panos Toulis, David C. Parkes.
# Author: Panos Toulis(ptoulis@fas.harvard.edu)

## Simple binomial test.
binomial.test <- function(pA=0.25, pB=0.2, nreps=1, ntotal=10^4,
                          levels=c(0, 1), probs=c(1,1)) {
  
  data = sample(levels, size=ntotal, prob=probs, replace=T)
  print(table(data))
  ntest.units <- as.integer(seq(1, ntotal, length.out=ntotal/10))
  Wins = matrix(NA, nrow=nreps, ncol=length(ntest.units))
  Loss = matrix(NA, nrow=nreps, ncol=length(ntest.units))
  
  for(j in 1:nreps) {
    playerA = rbinom(length(data), size=data, prob=pA)
    playerB = rbinom(length(data), size=data, prob=pB)
    
    order.A = order(playerA, runif(length(playerA)))
    order.B = order(playerB, runif(length(playerB)))
    
 
    Awins = sapply(ntest.units, function(i) {
      set.A = tail(order.A, i)
      set.B = tail(order.B, i)
      sum(data[set.A]) > sum(data[set.B])
    })
    
    Aloss = sapply(ntest.units, function(i) {
      set.A = tail(order.A, i)
      set.B = tail(order.B, i)
      sum(data[set.A]) < sum(data[set.B])
    })
    

    Wins[j, ] <- Awins
    Loss[j, ] <- Aloss
  }
  
  plot(ntest.units, colMeans(Loss), xlab=sprintf("#test units (total %d)", ntotal), ylab="probability", col="red", type="l", ylim=c(0, 1))
  lines(ntest.units, colMeans(Wins), col="green")
  
}