## Causal Effects on Incentives.
## Copyright, Panos Toulis ptoulis@fas.harvard.edu
##
##  Terminology files describe the terms and lingo we use in the code.
##
##  Agent:
##    Agent is an entity that is assigned a unique ID > 0.
#     In a randomization an agent is assigned to a mechanism.
#     There are N total agents (assumed even)
##  Game:
##    Game is a data frame with the following:
##      nt,  nc, Ut, Uc, hid, hsize, mid
##    nt = #times strategy=t was tried
##    nc = #times strategy=c was tried
##    Ut = total utility from trying "t"
##    Uc = total utility from trying "c"
##    hid= Hospital id. Unique integer
##    hsize = hospital size
##    mid = mechanism id > 0
#
#   # Z = assignment, i.e. data frame
#     hid   mid
#     1     0 
#     25    1
#       ...

# checks-package. Get it from github
source("../r-toolkit/checks.R")

kFields <- c("nt", "nc", "Ut", "Uc", "hid", "hsize", "mid")
kGiDensityFile <- "cache/gi.density.Rdata"
kUFrameFileRCM <- "cache/Uframe-rCM.Rdata"
kUFrameFileXCM <- "cache/Uframe-xCM.Rdata"

kNoHospitals <- 6
kHospitalSize <- 20

random.game <- function(N, hospital.size) {
  ## Creates an empty game
  out <- matrix(0, nrow=N, ncol=length(kFields))
  colnames(out) <- kFields
  out = data.frame(out)
  out$nt <- rep(1, N)
  out$nc <- rep(1, N)
  out$hid <- 1:N
  out$hsize = rep(hospital.size, N)
  out$Ut <- rep(0.5 * hospital.size, N)
  out$Uc <- rep(0.5 * hospital.size, N)
  out$mid <- c(rep(0, N/2), rep(1, N/2))  # no need to randomize.
  return(out)
}

ucb.strategy <- function(game, hospital.id) {
  # Given a specific GAME and hospital id
  # picks the best strategy based on the UCB 1 rule.
  #
  CHECK_game(game)
  CHECK_MEMBER(hospital.id, game$hid)  # valid hospital id?
  hentry <- subset(game, hid==hospital.id)
  mean.ut <- hentry$Ut / hentry$nt
  mean.uc <- hentry$Uc / hentry$nc
  n <- hentry$nt + hentry$nc
  choice.t <- mean.ut + sqrt(2 * log(n) / hentry$nt)
  choice.c <- mean.uc + sqrt(2 * log(n) / hentry$nc)
  if (choice.t > choice.c) {
    return ("t")
  } else if (choice.t < choice.c) {
    return ("c")
  } else {
    return (sample(c("c", "t"), size=1))
  }
}

update.game <- function(game, hospital.id, ucb.choice, utility) {
  # Given the choice and the realized utility, it will update the 
  # GAME for this particular hospital.
  # i.e. sum up the utilities and increment counters.
  CHECK_game(game)
  CHECK_MEMBER(hospital.id, game$hid)  # valid hospital id?
  ucbs <- replicate(10, ucb.strategy(game, hospital.id))
  CHECK_MEMBER(ucb.choice, ucbs)
  #
  # Now update the utility
  row.id <- which(game$hid==hospital.id)
  if (ucb.choice=="t") {
    game[row.id, ]$nt <- game[row.id, ]$nt + 1
    game[row.id, ]$Ut <- game[row.id, ]$Ut + utility
  } else {
    game[row.id, ]$nc <- game[row.id, ]$nc + 1
    game[row.id, ]$Uc <- game[row.id, ]$Uc + utility
  }
  return(game)
}

CHECK_game <- function(game) {
  CHECK_GE(nrow(game), 0)
  CHECK_SETEQ(names(game), kFields)
  # nt, nc>0  (should be initialized to 1.)
  CHECK_TRUE(all(game$nt > 0))
  CHECK_TRUE(all(game$nc >  0))
  CHECK_UNIQUE(game$hid)
}