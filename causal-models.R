## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
## Causal inference models.

print("Loading Cache...")
load(file=kFiDensityFile)  # loads fi.density
load(file=kGiDensityFile)  # loads gi.density
load(file=kUFrameFileRCM)  # loads  Urcm
load(file=kUFrameFileXCM)  # loads Uxcm


h.contrast <- function(Y1, Y0) {
  # constrast functions.
  CHECK_EQ(length(Y0), length(Y1))
  CHECK_MEMBER(Y0, c(0,1))
  CHECK_MEMBER(Y1, c(0,1))
  return(mean(Y1-Y0))
}

log.fi <- function(Rij) {
  # log density under deviation
  CHECK_rke(Rij)
  tab <- summary.Rij(Rij)
  dens <- fi.density$density
  probs.t <- sapply(names(tab), function(ps) dens[[ps]])
  x = as.numeric(tab)
  fi.log = dmultinom(x, size=sum(x), prob=probs.t, log=T)
  return(fi.log)
}

log.gi <- function(Rij) {
  # log density under deviation
  CHECK_rke(Rij)
  tab <- summary.Rij(Rij)
  dens = gi.density$density
  probs.c <- sapply(names(tab), function(ps) dens[[ps]])
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
bootstrap.estimator <- function(ucb.out, nboots=1000) {
  # Rj.obs = LIST of RKE objects (reports)
  # Z = assignment, i.e. data frame
  #     hid   mid
  #     1     0 
  #     25    1
  #       ...
  # Assume:  subset(Z, mid==1)  gives the hospital ids H that are 
  # in the SAME order as in the R0.obs list
  R0.obs <- ucb.out$R0.obs
  R1.obs <- ucb.out$R1.obs
  CHECK_EQ(length(R0.obs), length(R1.obs))
  m = length(R1.obs)
  
  p1.obs <- sapply(1:m, function(i) 1 / (1 + gi.over.fi(R1.obs[[i]])))
  p0.obs <- sapply(1:m, function(i) 1 / (1 + gi.over.fi(R0.obs[[i]])))
  p0 <- c(round(p0.obs, 2), rep("-", m))
  p1 <- c(rep("-", m), round(p1.obs, 2))
  
  for (i in 1:length(R0.obs)) {
    logdebug(sprintf("Report of hospital %d in M0. Real strategy=%s p0i=%.3f",
                     i, ucb.out$s0[i], p0.obs[i]))
    logdebug(summary.Rij(R0.obs[[i]]))
    logdebug(sprintf("Report of hospital %d in M1. Real strategy=%s p1i=%.3f",
                     i, ucb.out$s1[i], p1.obs[i]))
    logdebug(summary.Rij(R1.obs[[i]]))
  }
  boot.samples <- c()
  
  logdebug("s0, p0")
  logdebug(c(ucb.out$s0, rep("-", m)))
  logdebug(p0)
  logdebug("s1, p1")
  logdebug(c(rep("-", m), ucb.out$s1))
  logdebug(p1)
  logdebug(sprintf("expected diff=%.3f", mean(p1.obs) - mean(p0.obs)))
  
  for (t in 1:nboots) {
    y1.obs <- rbinom(n=length(p1.obs), size=1, prob=p1.obs)
    y0.obs <- rbinom(n=length(p0.obs), size=1, prob=p0.obs)
    
    y1.mis <- sample(y1.obs, size=m, replace=T)
    y0.mis <- sample(y0.obs, size=m, replace=T)
    
    Y1 = c(y1.mis, y1.obs)
    Y0 = c(y0.obs, y0.mis)
    logfine("Imputation of Y0")
    logfine(Y0)
    logfine("Imputation of Y1")
    logfine(Y1)
    logfine(h.contrast(Y1, Y0))
    boot.samples[t] <- h.contrast(Y1, Y0)
  }
  return(boot.samples)
}

## Estimation based on the game-theoretic prior.
gtprior.estimator <- function(ucb.out, mcmc.niters=1000) {
  # in the SAME order as in the R0.obs list
  R0.obs <- ucb.out$R0.obs
  R1.obs <- ucb.out$R1.obs
  CHECK_EQ(length(R0.obs), length(R1.obs))
  m = length(R1.obs)
  
  gifi.obs1 <- sapply(1:m, function(i) gi.over.fi(R1.obs[[i]]))
  gifi.obs0 <- sapply(1:m, function(i) gi.over.fi(R0.obs[[i]]))
    
  logdebug("s0, g0f0")
  logdebug(c(ucb.out$s0, rep("-", m)))
  logdebug(gifi.obs0)
  logdebug("s1, g1f1")
  logdebug(c(rep("-", m), ucb.out$s1))
  logdebug(gifi.obs1)
  
  # Run the MCMC chain
  D = 2 * m  # dimension
  beta = 1.5  # picking this correctly will improve performance.
  CHECK_EQ(D, kNoHospitals)
  
  Y0.mcmc <- matrix(0, nrow=D, ncol=mcmc.niters)
  Y1.mcmc <- matrix(0, nrow=D, ncol=mcmc.niters)
  Y0.mcmc[, 1] <- rep(0, D)
  Y1.mcmc[, 1] <- rep(0, D)
  
  rcm.prob <- function(nt, factor=1) {
    alpha.c <- Urcm$c[nt + 1]
    alpha.t <- Urcm$t[nt + 2]
    return(1 / (1 + factor * exp(beta * (alpha.c - alpha.t))))
  }
  xcm.prob <- function(nt, factor=1) {
    alpha.c <- Uxcm$c[nt + 1]
    alpha.t <- Uxcm$t[nt + 2]
    return(1 / (1 + factor * exp(beta * (alpha.c - alpha.t))))
  }
  
  for (t in 2:mcmc.niters) {
    Y0.t <- Y0.mcmc[, t-1]
    Y1.t <- Y1.mcmc[, t-1]
    for (d in 1:D) {
      # Gibbs.  Update each component
      # 1. How many truthful in all but d-th component?
      N0t <- sum(Y0.t[-d])
      N1t <- sum(Y1.t[-d])
      # 2. Compute conditional Bernoulli distribution
      r0.prob <- 0
      r1.prob <- 0
      if (d <= m){
        # This hospital is in M0
        r0.prob <- rcm.prob(N0t, factor=gifi.obs0[d])
        r1.prob <- xcm.prob(N1t)                   
      } else {
        # this hospital was in M1
        r0.prob <- rcm.prob(N0t)
        r1.prob <- xcm.prob(N1t, factor=gifi.obs1[d-m])
      }
#       logdebug("Y0, Y1=")
#       logdebug(Y0.t)
#       logdebug(Y1.t)
#       logdebug(sprintf("#truthful N0t=%d   N1t=%d", N0t, N1t))
#       logdebug(sprintf("Truthful prob pRCM.t=%2f, pxCM.t=%2f", r0.prob, r1.prob))
      # 3. Sample from conditional
      Y0.t[d] <- rbinom(1, size=1, prob=r0.prob)
      Y1.t[d] <- rbinom(1, size=1, prob=r1.prob)
    }
    #  4. Update the Gibbs chains
    Y0.mcmc[, t] <- Y0.t
    Y1.mcmc[, t] <- Y1.t
  }
  burnin <- as.integer(0.1 * mcmc.niters)
  m0 = colMeans(Y0.mcmc[, -c(1:burnin)])
  m1 = colMeans(Y1.mcmc[, -c(1:burnin)])
  return(m1 - m0)
}

