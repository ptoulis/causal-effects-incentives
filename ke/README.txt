## Notes for experiments.

1) To run the UCB dynamics and get the incentives estimand for a specific T run:
   Tau.UCBEstimands.Experiment(nsamples=100, stopping.T=....)
   
   NOTE: Saves in "out/UCBEstimand..."
   
2) To produce the "separability" gi/fi plot run
  Gi.Fi.Experiment(nsamples=1000)
  gifi.exp()

  NOTE: this is using the gi.density/fi.density globs for computing gi/fi(R) ratios
  Saves in "out/GiFi-experiment.."
  
3) To produce the U-frames for rCM or xCM run:
    compute.U.frame(mech="...", hospital.size=kHospitalSize, nhospitals=kNoHospitals, nsims=200, verbose=T)
   
   NOTE: Saves "Urcm/Uxcm" objects in "cache/"

4) To produce constrast estimates of Game-Theoretic (GT) and Boostrap approach
  BootstrapVsGT.Experiment(samples.per.T=100)
  
  NOTE: saves "boot.estimates/gt.estimates" in "out/Estimation-..."
  
5) To run the  contrast with BAD gi/fi run
  messup.gi.fi()
  BootstrapVsGT.Experiment(samples.per.T=100)
