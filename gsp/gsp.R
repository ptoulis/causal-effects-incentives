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
    mbid = 0
    if(s > length(bids.minus.j)) {
      mbid = 0
    } else {
      mbid = bids.minus.j[s]
    }
    return(max(auctionConfig$reservePrice, mbid))
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
  target.payment = get.min.bid.for.slot(target.slot)
  loginfo(sprintf("Best slot = %d  Target payment=%.2f", 
                  target.slot, target.payment))
  # 5. Compute b' = best bid for that target slot.
  target.bid = NA
  if(target.slot==1) {
    target.bid = 0.5 * (Uj + target.payment)
  } else {
    Wj = auctionConfig$CTR[target.slot] / auctionConfig$CTR[target.slot-1]

    target.bid = (1-Wj) * Uj + Wj * target.payment
    loginfo(sprintf("Uj =%.3f Wj=%.3f target bid=%.3f", Uj, Wj, target.bid))
  }
  
  bids = bids.prev
  bids[agent.id] <- min(Uj, target.bid)  # don't bid more than valuation
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
  slot.allocation = sapply(1:nagents, function(i) {
    slot.i = which(winners==i)
    if(slot.i > auctionConfig$nslots) return(NA)
    return(slot.i)
  })
  # Current slot allocation = NA NA NA 1 2 (a_i = slot allocation of agent i)
  # Apply the reserve price. Make NA if bid < reserve
  slot.allocation = sapply(1:length(slot.allocation), function(aid) {
    # aid = agent id
    slot = slot.allocation[aid]
    ifelse(bids[aid] >= auctionConfig$reservePrice, slot, NA)
  })

  payments = sapply(slot.allocation, function(s) {
    ifelse(is.na(s), 0, ifelse(s+1 <= nagents, bids[winners[s+1]], 0))
  })
  return(list(allocation=slot.allocation, payments=payments))
}

auc = auctionConfig.example()
kCurrentLogLevel <- 4

run.GSP.rounds <- function(auctionConfig, nrounds) {
  CHECK_auctionConfig(auctionConfig)
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

reserve.price.revenue <- function(auctionConfig, nReserves, nrounds) {
  reserve = seq(0, max(auctionConfig$valuations)+1, length.out=nReserves)
  cols = rainbow(length(reserve))
  for(i in 1:length(reserve)) {
    r = reserve[i]
    auctionConfig$reservePrice <- r
    x = run.GSP.rounds(auctionConfig, nrounds)
    revenue = auctionData.revenue(x)
    print(sprintf("Plotting for reserve=%.3f", r))
    if(i==1) 
      plot(revenue, type="l", ylim=c(0, 2*max(revenue)))
    else
      lines(revenue, col=cols[i])
  }
}





