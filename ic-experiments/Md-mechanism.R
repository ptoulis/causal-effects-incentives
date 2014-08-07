## Copyright 2013 Panos Toulis, Donald B. Rubin
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Incentive compatible experiments.
# M(d) mechanism is a very simple one. Assume we have V total units,
# and the mechanism picks S of them at random as seeds.
# The superset of test units is V-S.
# M(d) simply determines that the agents will pick d out of them
# as the set of units to be tested.
#
rm(list=ls())
source("../../r-toolkit/checks.R")
logistic <- function(x) exp(x) / (1 + exp(x))

sample.network <- function(V=1000, S=100) {
  # Samples the binary outcomes Y on the test population
  #
  # Returns: data frame   Y     I    L
  # where Y={0, 1} where the unit adopted/bought or not
  #       I = {0, 1, 2...}  #incoming links to a unit i
  #       L = probabilities of adoption given I
  #
  # We assume iid P(Yi=1 | Ii) = exp(b0 + b1 Ii) / (1+exp(b0+b1 Ii)
  # where Ii = #incoming links to unit i
  #
  # TODO(ptoulis): Definition of the network and the model 
  # of Y can be more complicated.
  test.units = V-S
  # sample the incoming links
  Incoming = rpois(test.units, lambda=0.4)
  # Sample the Y
  # determine the parameters of influence
  b0 = -3.1  # corresponds to 0.04% baseline rate
  b1 = 0.9
  L = logistic(b0 + b1 * Incoming)
  Y = rbinom(test.units, size=1, prob=L)
  df = data.frame(Y=Y, I=Incoming, L=L)
  return(df)
}

agent.decision <- function(true.network, d, theta) {
  ## Given some parameters, agents observe a noisy version of the network
  # In this implementation we assume theta=P(edge is right)
  #
  # Args:
  #   true.network = the actual network information (Y, I, L)
  #   d = max # of test units for the agent to pick
  #   theta = parameter of the model of the agent
  # 
  # Returns:
  #   Vector of length d with unit ids chosen as the test set.
  #
  # TODO(ptoulis): More options are available for the noise model.
  # Also agents might have different models to choose the test set.
  CHECK_INTERVAL(theta, 0, 1)
  N = nrow(true.network)
  new.I = rbinom(N, size=true.network$I, prob=theta)
  # print(mean((new.I - true.network$I)^2))
  CHECK_TRUE(d <= length(new.I))
  units = tail(order(new.I), d)
}

run.Md <- function(network, d, pA, pB) {
  # Returns the scores of the agents.
  # print(sprintf("Checking agent A pA=%.3f", pA))
  units.A = agent.decision(network, d, pA)
  # print(sprintf("Checking agent B pB=%.3f", pB))
  units.B = agent.decision(network, d, pB)
  
  score.A = sum(network$Y[units.A])
  score.B = sum(network$Y[units.B])
  return(list(score.A=score.A, score.B=score.B))
}

optimal.Md <- function(V=10^5, S=10^3, nsamples=100) {
  pA = 0.25
  pB = 0.2
  V = 10^5
  S = 10^3
  G = sample.network(V, S)
  
  d.values = 10^4 * c(0.001, 0.005, 0.01, 0.03, 0.05, 0.08, 
                      0.1, 0.12, 0.15, 0.18, 0.21, 
                      0.3, 0.35, 0.4, 0.5, 0.7, 0.9, 1.0,
                      7, 7.5, 8.1, 8.7, 9, 9.15, 9.2, 9.5, 9.8)
 
  prob.win.A = c()#rep(NA, length(d.values))
  prob.lose.A = c()#rep(NA, length(d.values))
  pb = txtProgressBar(style=3)
  
  for(i in 1:length(d.values)) {
    d = d.values[i]
    A.wins = c()
    A.loss = c()
    for(j in 1:nsamples) {
      scores = run.Md(G, d, pA, pB)
      A.wins = c(A.wins, scores$score.A > scores$score.B)
      A.loss = c(A.loss, scores$score.A < scores$score.B)
    }
    prob.win.A = c(prob.win.A, mean(A.wins))
    prob.lose.A = c(prob.lose.A, mean(A.loss))
    setTxtProgressBar(pb, value=i/length(d.values))
    plot(head(d.values/10^4, i), prob.win.A, type="l", col="green", lwd=1.5, ylim=c(0,1))
    lines(head(d.values/10^4, i), prob.lose.A, col="red", lwd=1.2)
  }
  
}




