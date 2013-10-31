## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
## Experiments:
##
## This loads all KE code.
rm(list=ls())
library(logging)
basicConfig()
removeHandler("basic.stdout")

source("../r-toolkit/checks.R")
source("../kidney-exchange/terminology.R")
source("../kidney-exchange/rke.R")
source("../kidney-exchange/matching.R")
source("../kidney-exchange/mechanisms.R")

source("terminology.R")
source("offline.R")  # will load offline the gi(.) density + "Utility frame"
source("causal-models.R")

# Note:
# m=6 hospitals, n=20 pairs -> all truthful, rCM clears in 0.5sec
# 
Run.Mechanisms.UCB <- function(game, niters) {
  # This will run the mechanism for niters steps, with agents following UCB policy
  #
  # game = GAME object. Look at terminology
  # niters = steps to run the mechanism
  # nhospital.size = # pairs for each hospital
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
    print("")
    print(sprintf("Mechanism %s iters=%d, m=%d, size=%d", mech, niters, m, hsize))
    print(game)
    kpd <- NA
    strategies <- NA
    for (t in 1:niters) {
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
      setTxtProgressBar(pb, value=t/niters)
    }
    
    if (mech.id==1) {
      out$R1.obs <- kpd$reported.pool$rke.list
      out$h1.ids <- hids
      out$M1 <- mech
    } else {
      out$R0.obs <- kpd$reported.pool$rke.list
      out$h0.ids <- hids
      out$M0 <- mech
    }
  } # for all mechanisms
  out$Z <- subset(game, select=c("hid", "mid"))
  out$game <- game
  return(out)  
}



Bootstrap.Experiment <- function(nboots=100) {
  g <- random.game(N=kNoHospitals, hospital.size=kHospitalSize)
  # watch out. These need to match to our computed cache.
  print("Running Mechanisms + UCB")
  ucb.out <- Run.Mechanisms.UCB(g, niters=40)
 
  H0.ids <- ucb.out$h0.ids
  H1.ids <-ucb.out$h1.ids
  CHECK_SETEQ(H0.ids, subset(g, mid==0)$hid)
  
  R0.obs <- ucb.out$R0.obs
  R1.obs <- ucb.out$R1.obs
  
  Z = ucb.out$Z
  tau.boot <- bootstrap.estimator(ucb.out, nboots=nboots)
  return(tau.boot)
}