## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
## Experiments:
##
## This loads all KE code.
library(logging)
basicConfig()
removeHandler("basic.stdout")

source("terminology.R")
source("../kidney-exchange/terminology.R")
source("../kidney-exchange/rke.R")
source("../kidney-exchange/matching.R")
source("../kidney-exchange/mechanisms.R")

Run.Mechanism.UCB <- function(game, mech, mech.id, niters,
                              uniform.pra=T, include.3way=F) {
  # This will run the mechanism for niters steps, with agents following UCB policy
  #
  # game = GAME object. Look at terminology
  # mech = name of mechanism, e.g. "xCM", "rCM"
  # mech.id = index of mechanism. We will take hospitals from GAME that match this id
  # niters = steps to run the mechanism
  # nhospital.size = # pairs for each hospital
  # uniform.pra/include.3way = self-explanatory. Don't change.
  pb <- txtProgressBar(style=3)
  CHECK_game(game)
  hset <- subset(game, mid==mech.id)
  hids <- hset$hid
  m = nrow(hset)  # no. of hospitals
  CHECK_SETEQ(diff(hset$hsize), c(0))  # each hospital be the same size
  hsize = head(hset$hsize, 1)  #
  
  loginfo(sprintf("Running for %d iterations.", niters))
  loginfo(sprintf("No. of hospitals=%d  size=%d", m, hsize))
  kpd <- NA
  strategies <- NA
  for (t in 1:niters) {
    
    # 1. Sample a RKE pool
    rke.pool = rrke.pool(m=m, n=hsize, uniform.pra=uniform.pra)
    # 2. Run UCB to get strategies from hospitals
    #   IMPORTANT note:  vector (hids) is mapped to  1,2,3,...m
    strategies <- sapply(hids, function(id) ucb.strategy(game=game, hospital.id=id))
    strategy.str <- paste(strategies, collapse="")
    # loginfo(sprintf("Strategies %s", strategy.str))
    # 3. Create the KPD and then run the mechanism.
    kpd <- kpd.create(rke.pool, strategy.str=strategy.str)
    matching <- Run.Mechanism(kpd, mech=mech, include.3way=include.3way)
    # 4. Compute the utilities and update the game for each hospital.
    U <- get.matching.hospital.utilities(matching, m)
    for(i in 1:m) {
      # loginfo(sprintf("Updating for agent %d util=%.2f", hids[i], U[i]))
      game <- update.game(game, hids[i], ucb.choice=strategies[i], utility=U[i])
    }
    setTxtProgressBar(pb, value=t/niters)
  }
  # Return the KPD object.
  # This has reports and 
  return(list(reports=kpd, Z=subset(game, select=c("hid", "mid")), 
              strategies=strategies, game=game))
}



Bootstrap.Experiment <- function() {
  
}