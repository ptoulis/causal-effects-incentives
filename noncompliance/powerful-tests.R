# Based on (DBR 1998, "More powerful Randomization-based p-values in double-blind trials with noncompliance")
source("../../r-toolkit/checks.R")
TYPE_COMPLIER = 1
TYPE_NONCOMPLIER = 0

N = 100
Science = data.frame(Unit=c(1:N), C=rep(NA, N), Z=rep(NA, N),
                     Y0=rep(NA, N), Y1=rep(NA, N))
Science$C = sample(c(rep(TYPE_COMPLIER, 0.03*N), rep(TYPE_NONCOMPLIER, 0.97 * N))) # compliers
Science$Z = sample(c(rep(1,N/2), rep(0, N/2)))

CHECK_EQ(sum(Science$Z), N/2)
CHECK_EQ(length(Science$C), N)

science.noncompliers <- function(sc) {
  subset(sc, C==TYPE_NONCOMPLIER)$Unit
}

science.compliers <- function(sc) {
  subset(sc, C==TYPE_COMPLIER)$Unit
}

generate.outcomes <- function(nunits, levels, distr) {
  ## Generates a vector of "length" with outcomes that have "levels"
  ## at the specified distribution.
  ## e.g. generate.outcomes(10, 3, c(0.5, 0.4, 0.1)) 
  ##    gives outcomes 0,1,2 on 10 units, with a distribution 50%, 40%, 10% resp.
  distr = distr / sum(distr)
  CHECK_EQ(length(distr), length(levels), msg="#units=#distribution points")
  freq = as.integer(nunits * distr)
  # freq = (15, 34, 44)
  x = c()
  for(i in 1:length(freq)) {
    y = levels[i]
    f = freq[i]
    x = c(x, rep(y, f))
  }
  return(sample(x))
}



