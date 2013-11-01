## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
## Experiments:
##
## This loads all KE code.
rm(list=ls())

source("../r-toolkit/checks.R")
source("../kidney-exchange/terminology.R")
source("../kidney-exchange/rke.R")
source("../kidney-exchange/matching.R")
source("../kidney-exchange/mechanisms.R")

source("terminology.R")
source("offline.R")  # will load offline the gi(.) density + "Utility frame"
source("causal-models.R")

logReset()
kLogFile <- "out/logs.txt"
kCurrentLevel  <- 5 # debug
#  FINE DEBUG  INFO  WARNING   ERROR
#   0     1      2     3        4
library(stringr)
if (file.exists(file=kLogFile)) {
  file.remove(file=kLogFile)
}

logthis <- function(x, level=kCurrentLevel) {
  wline <- function(vec) {
    write(str_c(vec, collapse="\t"), file=kLogFile, append=T)
  }
  if (level >= kCurrentLevel) {
    preamble <- sprintf("%s::", date())
    if(is.character(x)) {
      wline(c(preamble, x))
    } else if(is.data.frame(x)) {
      wline(c(preamble, "Data frame "))
      wline(c("\t", names(x)))
      logthis(as.matrix(x), level=level)
    } else if(is.matrix(x)) {
      for(i in 1:nrow(x))
        wline(c("\t", as.vector(x[i, ])))
    } else if (is.table(x)) {
      wline(c(preamble, "Table"))
      wline(names(x))
      wline(as.numeric(x))
    } else if(is.vector(x)) {
      wline(c(preamble, x))
    }
  }
}

logfine <- function(x) logthis(x, level=0)
logdebug <- function(x) logthis(x, level=1)
loginfo <- function(x) logthis(x, level=2)
logwarning <- function(x) logthis(x, level=3)
logerror <- function(x) logthis(x, level=4)



# Note:
# m=6 hospitals, n=20 pairs -> all truthful, rCM clears in 0.5sec
# 
Run.Mechanisms.UCB <- function(game, stopping.T, verbose=F) {
  # This will run the mechanism for niters steps, with agents following UCB policy
  #
  # game = GAME object. Look at terminology
  # T = time to stop and return agent reports.
  # nhospital.size = # pairs for each hospital
  pb <- NA
  if(verbose)
    pb <- txtProgressBar(style=3)
  CHECK_game(game)
  out = list()
  all.mechanisms <- c("rCM", "xCM")
  for (mech.id in c(0,1)) {
    mech = all.mechanisms[mech.id + 1]
    hset <- subset(game, mid==mech.id)
    hids <- hset$hid
    m = nrow(hset)  # no. of hospitals
    CHECK_SETEQ(diff(hset$hsize), c(0))  # each hospital be the same size
    hsize = head(hset$hsize, 1) 
    if(verbose) {
      print("")
      print(sprintf("Mechanism %s iters=%d, m=%d, size=%d", mech, stopping.T, m, hsize))
      print(game)
    }
    kpd <- NA
    strategies <- NA
    for (t in 1:stopping.T) {
      # 1. Sample a RKE pool
      rke.pool = rrke.pool(m=m, n=hsize, uniform.pra=T)
      # 2. Run UCB to get strategies from hospitals
      #   IMPORTANT note:  vector (hids) is mapped to  1,2,3,...m
      strategies <- sapply(hids, function(id) ucb.strategy(game=game, hospital.id=id))
      strategy.str <- paste(strategies, collapse="")
      # loginfo(sprintf("Strategies %s", strategy.str))
      # 3. Create the KPD and then run the mechanism.
      kpd <- kpd.create(rke.pool, strategy.str=strategy.str)
      matching <- Run.Mechanism(kpd, mech=mech, include.3way=F)
      # 4. Compute the utilities and update the game for each hospital.
      U <- get.matching.hospital.utilities(matching, m)
      for(i in 1:m) {
        # loginfo(sprintf("Updating for agent %d util=%.2f", hids[i], U[i]))
        game <- update.game(game, hids[i], ucb.choice=strategies[i], utility=U[i])
      }
      if(verbose) setTxtProgressBar(pb, value=t/stopping.T)
    }
    
    if (mech.id==1) {
      out$R1.obs <- kpd$reported.pool$rke.list
      out$h1.ids <- hids
      out$M1 <- mech
      out$s1 <- strategies
    } else {
      out$R0.obs <- kpd$reported.pool$rke.list
      out$h0.ids <- hids
      out$M0 <- mech
      out$s0 <- strategies
    }
  } # for all mechanisms
  out$Z <- subset(game, select=c("hid", "mid"))
  out$game <- game
  return(out)  
}


Gi.Fi.Experiment <- function(nsamples=100) {
  ## Compute the statistics for gi(R) / fi(R)
  ## where R is the report, when agents are truthful and when not.
  GiFi.t <- c() # under truthful
  GiFi.c <- c()
  pb = txtProgressBar(style=3)
  for (t in 1:nsamples) {
    Rt <- sample.Rij(kHospitalSize, Yij=1)
    Rc <- sample.Rij(kHospitalSize, Yij=0)
    
    GiFi.t = c(GiFi.t, gi.over.fi(Rt))
    GiFi.c = c(GiFi.c, gi.over.fi(Rc))
    setTxtProgressBar(pb, value=t / nsamples)
  }
  loginfo("Done")
  print(summary(GiFi.t))
  print(summary(GiFi.c))
  GiFi <- list(t=GiFi.t, c=GiFi.c)
  save(GiFi, file="out/GiFi-experiment.Rdata")
}

gifi.exp <- function(thres.upper=1000, thres.lower=0, breaks=30) { 
  load(file="out/GiFi-experiment.Rdata")
  x = GiFi$t
  y = GiFi$c
  x = 1 / (1 + x[x >= thres.lower & x < thres.upper])
  y = 1 / (1 + y[y >= thres.lower & y < thres.upper])
  hist(x, breaks=breaks, freq=F, main="P(Y=1 | R)", xlab="gi/fi", col=rgb(1, 0, 0, 0.5)); 
  hist(y, breaks=breaks, freq=F, col=rgb(0, 1, 0, 0.3), add=T)
}

# Experiments to measure incentives under the UCB dynamics.
Tau.UCBEstimands.Experiment <- function(stopping.T, nsamples=10) {
  ucb.estimand = list(M0=rep(0,0), M1=rep(0, 0))
  pb = txtProgressBar(style=3)
  
  for (i in 1:nsamples) {
    game <- random.game(N=2 * kNoHospitals, hospital.size=kHospitalSize)
    out <- Run.Mechanisms.UCB(game=game, stopping.T=stopping.T)
    g0 = subset(out$game, mid==0)
    g1 <- subset(out$game, mid==1)
    ucb.estimand$M0 <- c(ucb.estimand$M0, with(g0, nt / (nt + nc)))
    ucb.estimand$M1 <- c(ucb.estimand$M1, with(g1, nt / (nt + nc)))
    setTxtProgressBar(pb, value=i / nsamples)
    save(ucb.estimand, file=sprintf("out/UCBEstimand-t%d-experiment.Rdata", stopping.T))
    cat(sprintf("\nM1=%.3f  M0=%.3f  D=%.3f, se=%.3f",
                  mean(ucb.estimand$M1), mean(ucb.estimand$M0),
                  mean(ucb.estimand$M1) - mean(ucb.estimand$M0),
                  bootstrap.mean(ucb.estimand$M1)))
  }
  return(ucb.estimand)
}

# Experiments for the bootstrap estimator.
BootstrapVsGT.Experiment <- function(stopping.T=5,
                                 samples.per.T=10,
                                 nboots=1000,
                                 filename.suffix="default") {
  boot.estimates <- c()
  gt.estimates <- c()
  # Checks whether we are using the correct size of hospitals
  CHECK_TRUE(gi.density$size == kHospitalSize)
  CHECK_TRUE(fi.density$size == kHospitalSize)
  
  pb <- txtProgressBar(style=3)

  for (j in 1:samples.per.T) {
    game <- random.game()
    logdebug(sprintf("T=%d  j=%d ", stopping.T, j))
    logdebug("Game before running")
    logdebug(game)
    # watch out. These need to match to our computed cache.
    # 1. Run the UCB dynamics
    ucb.out <- Run.Mechanisms.UCB(game, stopping.T=stopping.T)
    logdebug("Game after running")
    logdebug(ucb.out$game)
    # 2a. Run the bootstrap estimator  
    tau.boot <- bootstrap.estimator(ucb.out, nboots=nboots)
    # 2b. Run the GT estimator
    tau.gt <- gtprior.estimator(ucb.out, 1000)
    
    # 3. Save the estimates
    boot.estimates <- c(boot.estimates, mean(tau.boot))
    gt.estimates <- c(gt.estimates, mean(tau.gt))
    
    logdebug(sprintf("Current estimates Boot=%.3f  GT=%.3f", mean(boot.estimates), mean(gt.estimates)))
    setTxtProgressBar(pb, value=j/samples.per.T)
  }
  tau.boot <- mean(boot.estimates)
  se <- bootstrap.mean(boot.estimates)
  tau.gt <- mean(gt.estimates)
  se.gt <- bootstrap.mean(gt.estimates)
  cat(sprintf("\nFinal bootstrap CI (%.3f, %.3f)", tau.boot - 2 * se, tau.boot + 2 * se))
  cat(sprintf("\nFinal gt CI (%.3f, %.3f)", tau.gt - 2 * se.gt, tau.gt + 2 * se.gt))
  
  save(boot.estimates, file=sprintf("out/Estimation-Bootstrap-experiment-%s.Rdata", filename.suffix))
  save(gt.estimates, file=sprintf("out/Estimation-GT-experiment-%s.Rdata", filename.suffix))
  
}

GTPrior.Experiment <- function() {
  ## Game-theoretic prior.
  gt.estimates <- c()
  
  # Checks whether we are using the correct size of hospitals
  CHECK_TRUE(gi.density$size == kHospitalSize)
  CHECK_TRUE(fi.density$size == kHospitalSize)
  
  pb <- txtProgressBar(style=3)
  
  for (j in 1:samples.per.T) {
    game <- random.game()
    logdebug(sprintf("T=%d  j=%d ", stopping.T, j))
    logdebug("Game before running")
    logdebug(game)
    # watch out. These need to match to our computed cache.
    # 1. Run the UCB dynamics
    ucb.out <- Run.Mechanisms.UCB(game, stopping.T=stopping.T)
    logdebug("Game after running")
    logdebug(ucb.out$game)
    # 2. Run the bootstrap estimator  
    tau.boot <- bootstrap.estimator(ucb.out, nboots=nboots)
    est = mean(tau.boot)
    se = bootstrap.mean(tau.boot)
    logdebug(sprintf("j=%d Bootstrap (%.3f, %.3f)", j, est - 2 * se, est + 2 * se))
    boot.estimates <- c(boot.estimates, est)
    setTxtProgressBar(pb, value=j/samples.per.T)
  }
  tau.boot <- mean(boot.estimates)
  se <- bootstrap.mean(boot.estimates)
  cat(sprintf("\nFinal bootstrap estimte (%.3f, %.3f)", tau.boot - 2 * se, tau.boot + 2 * se))
  save(boot.estimates, file=sprintf("out/Estimation-Bootstrap-experiment-%s.Rdata", filename.suffix))
}




