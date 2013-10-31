## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
## Causal inference models.

source("../kidney-exchange/terminology.R")

kGiDensityFile <- "cache/gi.density.Rdata"
kUFrameFileRCM <- "cache/Uframe-rCM.Rdata"
kUFrameFileXCM <- "cache/Uframe-xCM.Rdata"

load(file=kGiDensityFile)
load(file=kUframeFileRCM)
load(file=kUframeFileXCM)


h.contrast <- function(Y1, Y0) {
  # constrast functions.
  CHECK_EQ(length(Y0), length(Y1))
  CHECK_MEMBER(Y0, c(0,1))
  CHECK_MEMBER(Y1, c(0,1))
  return(mean(Y1-Y0))
}

log.fi <- function(Rij) {
  # Log density of true types
  CHECK_rke(Rij)
  tab <- table(Rij$pairs$desc)
  probs.t <- sapply(names(tab), function(ps) subset(kPairs, desc==ps)$prob)
  x = as.numeric(tab)
  fi.log = dmultinom(x, size=sum(x), prob=probs.t, log=T)
  return(fi.log)
}

log.gi <- function(Rij) {
  # log density under deviation
  CHECK_rke(Rij)
  tab <- table(Rij$pairs$desc)
  probs.c <- sapply(names(tab), function(ps) gi.density[[ps]])
  x = as.numeric(tab)
  gi.log = dmultinom(x, size=sum(x), prob=probs.c, log=T)
  return(gi.log)
}

gi.over.fi <- function(Rij, verbose=F) {
  # This gives the odds:
  #   P(Rij | Yij=0) /  P(Rij | Yij=1)
  gi.log = log.gi(Rij)
  fi.log = log.fi(Rij)
  return (exp(gi.log - fi.log))
}

## Causal models/estimators beyond.
bootstrap.estimator <- function(Z, R0.obs, R1.obs, nboots=1000) {
  # Rj.obs = LIST of RKE objects (reports)
  # Z = assignment, i.e. data frame
  #     hid   mid
  #     1     0 
  #     25    1
  #       ...
  # Assume:  subset(Z, mid==1)  gives the hospital ids H that are 
  # in the SAME order as in the R0.obs list
  CHECK_EQ(length(R0.obs), length(R1.obs))
  m = length(R1.obs)
  H1.ids <- subset(Z, mid==1)$hid  # hospital ids in M1
  H0.ids <- subset(Z, mid==0)$hid  # hospital ids in M0
  
  p1.obs <- sapply(1:m, function(i) 1 / (1 + gi.over.fi(R1.obs[[i]])))
  p0.obs <- sapply(1:m, function(i) 1 / (1 + gi.over.fi(R0.obs[[i]])))
  boot.samples <- c()
  
  for (t in 1:nboots) {
    y1.obs <- rbinom(n=length(p1.obs), size=1, prob=p1.obs)
    y0.obs <- rbinom(n=length(p0.obs), size=1, prob=p0.obs)
    
    y1.mis <- sample(y1.obs, size=m, replace=T)
    y0.mis <- sample(y0.obs, size=m, replace=T)
    
    Y1 = c(y1.mis, y1.obs)
    Y0 = c(y0.obs, y0.mis)
    
    boot.samples[t] <- h.contrast(Y1, Y0)
  }
  return(boot.samples)
}



