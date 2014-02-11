# Panos Toulis,
# Causal effects of mechanism interventions
# Copyright (C) 2014
#
# AuctionConfig
#   Configuration of an auction (GSP).
#   It is a LIST{nslots, nagents, valuations, CTR, reservePrice} where
#     nslots = #of available slots (INT) -- K for brevity.
#     nagents = # of agents (INT) -- N for brevity
#     valuations = vector of valuations (length N)
#                  SHOULD be in decreasing order, U1 > U2 > ...
#     CTR = vector of click-through rates (length K)
#           SHOULD be in decreasing order, CTR1 > CTR2 ...
#     reservePrice = $$$ of reserve prices. Only bids > reserve are considered.
# 
# AuctionData
#   Data structure for holding data from running an GSP auction
#   IT is a LIST{bids, outcomes}.
#   Denote with M=#rounds in the auction. Then:
#     bids = (MxN matrix) bids per round. bids[i, j] = bid at round i from agent j
#     outcomes = (MxN matrix) outcomes. outcomes[i,j] = slot acquired from agent j at round i
#                 IF not allocated SHOULD BE = NA
#     payments = (MxN) matrix of payments per agent.
#
rm(list=ls())
source("../../r-toolkit/checks.R")
source("../../r-toolkit/logs.R")

CHECK_auctionConfig <- function(auctionConfig) {
  CHECK_MEMBER(c("nslots", "nagents", "valuations", "CTR", "reservePrice"),
               names(auctionConfig))
  CHECK_TRUE(all(auctionConfig$valuations > 0), msg="valuations are > 0")
  CHECK_TRUE(all(auctionConfig$CTR >= 0))
  CHECK_EQ(rev(sort(auctionConfig$valuations)),
           auctionConfig$valuations, msg="Config vals should be sorted")
  CHECK_EQ(rev(sort(auctionConfig$CTR)),
           auctionConfig$CTR, msg="Config should be sorted")
  CHECK_EQ(length(auctionConfig$valuations), auctionConfig$nagents)
  CHECK_EQ(length(auctionConfig$CTR), auctionConfig$nslots)
}

CHECK_auctionData <- function(auctionData) {
  CHECK_MEMBER(c("bids", "outcomes", "payments"), names(auctionData))
  CHECK_EQ(nrow(auctionData$bids), nrow(auctionData$outcomes))
  CHECK_EQ(ncol(auctionData$bids), ncol(auctionData$outcomes))
  CHECK_TRUE(all(auctionData$payments >= 0))
}

auctionData.revenue <- function(auctionData) {
  # Computes the revenue for every round
  # Returns a Mx1 vector of revenues
  #
  apply(auctionData$payments, 1, sum)
}

outcomes.CTR <- function(outcomes, auctionConfig) {
  # Computes MxN matrix of CTR achieved by each agent in each outcome
  # outcomes = MxN vector of outcomes
  #
  # Returns MxN matrix of CTR quantities.
  warning("Test for CTR")
  C = outcomes
  for(s in auctionConfig.allSlots()) {
    C[C==s] <- auctionConfig$CTR[s] 
  }
  C[is.na(C)] <- 0
  return(C)
}

auctionData.utility <- function(auctionData, auctionConfig) {
  # Computes a (MxN) matrix of agent utilities.
  # Uij = utility of agent j at round i
  #
  warning("Test for utility()")
  M = nrow(auctionData$bids)
  CTR = outcomes.CTR(auctionData$outcomes, auctionConfig)
  V = matrix(auctionConfig$valuations,
             nrow=M, ncol=auctionConfig$nagents,
             byrow=T)
  CHECK_TRUE(all(diff(V[,1])==0))
  return(theta * (V - auctionData$payments))
}

auctionConfig.example <- function() {
  q = 0.8
  nAgents = 50
  nSlots = 10
  return(list(nagents=nAgents,
              nslots=nSlots,
              valuations=rev(seq(0.1, 5, length.out=nAgents)),
              CTR=q^seq(1,nSlots),
              reservePrice=0))
}

auctionConfig.allSlots <- function(auctionConfig) {
  # Returns the vector of slot ids.
  return(1:length(auctionConfig$CTR))
}

auctionData.empty <- function(auctionConfig, nrounds) {
  out = list(bids=matrix(NA, nrow=nrounds, ncol=auctionConfig$nagents),
             outcomes=matrix(NA, nrow=nrounds, ncol=auctionConfig$nagents),
             payments=matrix(NA, nrow=nrounds, ncol=auctionConfig$nagents))
  return(out)
}

auctionData.add <- function(round, auctionData, bid, outcome, payment) {
  CHECK_TRUE(round <= nrow(auctionData$bids))
  CHECK_TRUE(all(bid >= 0), msg="bids > 0")
  CHECK_TRUE(all(payment >=0), msg="payments > 0")
  auctionData$bids[round, ] <- bid
  auctionData$outcomes[round, ] <- outcome
  auctionData$payments[round, ] <- payment
  return(auctionData)
}

is.static.equilbrium <- function(auctionData, auctionConfig) {
  # Checks whether the auction has reached the VCG equilibrium
  CHECK_auctionConfig(auctionConfig)
  bids = auctionData$bids[nrow(auctionData$bids),]
  CHECK_EQ(rev(sort(bids)), bids, msg="bids should be ordered in the equilibrium")
  almost.zero = function(x) {
    abs(x) < 1e-5
  }
  is.eq = sapply(1:auctionConfig$nagents, function(i) {
    if(i==1) {
      return(bids[i] > bids[i+1])
    } else if(i <= auctionConfig$nslots) {
      almost.zero(auctionConfig$CTR[i-1] * (auctionConfig$valuations[i] - bids[i]) -
        auctionConfig$CTR[i] * (auctionConfig$valuations[i] - bids[i+1]))
    } else {
      almost.zero(bids[i] - auctionConfig$valuations[i])
    }
  })
  return(all(is.eq))
}