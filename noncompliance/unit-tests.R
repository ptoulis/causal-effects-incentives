## unit-tests
source("../../r-toolkit/checks.R")

test.generate.outcomes <- function() {
  K = rpois(1, lambda=5)
  S = 1:K
  m = 10
  N = sum(S) * m
  x = generate.outcomes(N, levels=1:K, distr=S)
  i = sample(S, size=1)
  print(sprintf("Checking K=%d, N=%d units, i=%d", K, N, i))
  print(x)
  CHECK_EQ(length(which(x==i)), i * m)
}