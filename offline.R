compute.gi.density <- function(hospital.size, nsims=100) {
  # Computes the pdf of deviation
  gi.density <- NA
  pb = txtProgressBar(style=3)
  for (i in 1:nsims) {
    rke = rrke(hospital.size)
    m = max.matching(rke)
    rke.left <- rke.remove.pairs(rke, m$match$pair.id)
    if (i == 1) {
      gi.density <- table(rke.left$pairs$desc)
    } else {
      gi.density <- table(rke.left$pairs$desc) + gi.density
    }
    setTxtProgressBar(pb, value=i / nsims)
  }
  gi.density <- gi.density / sum(gi.density)
  save(gi.density, file="cache/gi.density.Rdata")
}

compute.U.frame <- function(mech, hospital.size, nhospitals, nsims=100, verbose=F) {
  # Creates the data frame U
  # where each row shows U$t = utility of truthful when have U$nt=truthful agents.
  # Nt should be from 0 to nhospitals
  #
  #    U =      t     c     Nt
  #                 ...
  #            8.5   7.5    3
  #           9.5    5.4    4 
  #           10.1    5.1   5
  #  This table shows that when Nt=3 truthful hospitals, the truthful hospital gets 8.5 on average
  # and the deviating gets 7.5. Inc Nt, i.e. having more truthful gets more util
  # and so the mechanism is incentive compatible.
  out = matrix(0, nrow=nhospitals+1, ncol=3)
  colnames(out) <- c("t", "c", "Nt")
  out <- data.frame(out) 
  out$Nt <- seq(0, nhospitals)
  pb = txtProgressBar(style=3)
  Nol = nrow(out) * nsims  # total iters
  count = 0  # simple counter
  for (nt.index in 1:nrow(out)) {
    for (i in 1:nsims) {
      nt = out$Nt[nt.index]  # no. of truthful hospitals
      # 1. Set nt agents to being truthful
      strategies = c(rep("t", nt), rep("c", nhospitals-nt))  # nt truthful, rest is dev.
      strategy = paste(strategies, collapse="")
      # 2. Sample a RKE pool, create KPD based on the strategy profile.
      rke.pool = rrke.pool(m=nhospitals, n=hospital.size, uniform.pra=T)
      kpd = kpd.create(rke.pool, strategy=strategy)
      # 3. Run mechanism and get utilities.
      matching = Run.Mechanism(mech=mech, kpd=kpd, include.3way=F)
      U = get.matching.hospital.utilities(matching, nhospitals)
      truthful.ids <- which(strategies=="t")
      dev.ids <- which(strategies=="c")
      # 4. Calculate avg(utility_t) and avg(utility_c)
      add.ut <- ifelse(length(truthful.ids) > 0, mean(U[truthful.ids]), 0)
      add.uc <- ifelse(length(dev.ids) > 0, mean(U[dev.ids]), 0)
      if(verbose)
        print(sprintf("strategy=%s mech=%s nt=%d ut=%.3f  uc=%.3f i=%d", 
                      strategy, mech, nt, add.ut, add.uc, i))
      out$t[nt.index] = out$t[nt.index] + add.ut
      out$c[nt.index] = out$c[nt.index] + add.uc
      count <- count + 1
      setTxtProgressBar(pb, value = count / Nol)
      print(out/nsims)
    }
  }
  out$t = out$t / nsims
  out$c = out$c / nsims
  if(mech == "rCM") {
    Urcm <- out
    save(Urcm, file=kUFrameFileRCM)
  }  else {
    Uxcm <- out
    save(Uxcm, file=kUFrameFileXCM)
  }
}



if(! file.exists(kGiDensityFile)) {
  print("File for Gi density does not exist. Need to compute it.")
  compute.gi.density(hospital.size=20, nsims=1000)
}

if(! file.exists(kUFrameFileRCM)) {
  print("File for U frame does not exist. Need to compute it.")
  compute.U.frame(mech="rCM", hospital.size=20, nhospitals=6, nsims=1000)
}

if(! file.exists(kUFrameFileXCM)) {
  print("File for U frame does not exist. Need to compute it.")
  compute.U.frame(mech="xCM", hospital.size=20, nhospitals=6, nsims=1000)
}
