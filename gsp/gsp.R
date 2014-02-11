# Panos Toulis,
# Causal effects of mechanism interventions
# Copyright (C) 2014
source("terminology.R")
bb.strategy <- function(agent.id, bids.prev, auctionConfig) {
  warning("TEST for bb.strategy()")
  CHECK_EQ(length(bids.prev), auctionConfig$nagents, msg="Correct #agents")
  bids.minus.j = rev(sort(bids.prev[-c(agent.id)])) # all bids but j (sorted)
  # Compute price paid to get slot "s"
  get.min.bid.for.slot <- function(s) {
    if(s > length(bids.minus.j))
      return(0)
    return(bids.minus.j[s])
  }
  # 1. All slots 1:K
  all.slots = auctionConfig.allSlots(auc)
  # 2. Uj = utility for agent j
  Uj = auctionConfig$valuations[agent.id]
  # 3. Utilities for all slots.
  slot.utils = sapply(all.slots, function(s) {
    auctionConfig$CTR[s] * (Uj - get.min.bid.for.slot(s))
  })
  loginfo(sprintf("Agent %d.", agent.id))
  loginfo("bids")
  loginfo(bids.prev)
  loginfo("Slot utilities")
  loginfo(slot.utils)
  # 4. Compute s* = best slot and the payment.
  target.slot = which.max(slot.utils)
  loginfo(sprintf("Best slot = %d", target.slot))
  target.payment = get.min.bid.for.slot(target.slot)
  # 5. Compute b' = best bid for that target slot.
  target.bid = NA
  if(target.slot==1) {
    target.bid = 0.5 * (Uj + target.payment)
  } else {
    Wj = auctionConfig$CTR[target.slot] / auctionConfig$CTR[target.slot-1]
    target.bid = (1-Wj) * Uj + Wj * target.payment
  }
  
  bids = bids.prev
  bids[agent.id] <- target.bid
  return(bids)
}
 
run.GSP <- function(bids, auctionConfig) {
  # Returns the allocation and the payment vectors
  # Both are length equal to #bids = #agents in config
  # 
  # Returns LIST(allocation, payments)
  # Working example. K=2 slots
  #            agents
  # bids =   1      2      3      4     5 
  #          0.5   0.3    0.6    0.9   0.7
  #           
  winners = rev(order(bids))  # 4 5 3 1 2
  nagents = length(bids)
  CHECK_TRUE(nagents==auctionConfig$nagents, msg="Correct #agents")
  # Current Slot allocation = 1 2 3 3 3
  slot.allocation = sapply(1:nagents, function(i) {
    slot.i = which(winners==i)
    if(slot.i > auctionConfig$nslots) return(NA)
    return(slot.i)
  })
  # Current slot allocation = NA NA NA 1 2 (a_i = slot allocation of agent i)
  payments = sapply(slot.allocation, function(s) {
    ifelse(is.na(s), 0, ifelse(s+1 <= nagents, bids[winners[s+1]], 0))
  })
  return(list(allocation=slot.allocation, payments=payments))
}

run.GSP.rounds <- function(auctionConfig, nrounds) {
  auctionData <- auctionData.empty(auctionConfig, nrounds=nrounds)
  current.bids = rep(0, auctionConfig$nagents)
  for(round in 1:nrounds) {
    agent.id = sample(1:auctionConfig$nagents, size=1)
    current.bids <- bb.strategy(agent.id, bids.prev=current.bids, auctionConfig)
    mech.out = run.GSP(current.bids, auctionConfig)
    auctionData <- auctionData.add(round, auctionData,
                                   current.bids, mech.out$allocation, mech.out$payments)
  }
  return(auctionData)
}
