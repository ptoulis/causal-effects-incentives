summary.Rij <- function(Rij) {
  # computes summary of the report (Needs to be a TABLE)
  table(factor(Rij$pairs$pair.type, levels=c("O", "U", "S", "R")))
}

sample.Rij <- function(nHospitalSize, Yij) {
  if (Yij==1) {
    return(rrke(nHospitalSize))
  } else {
    rke = rrke(nHospitalSize) 
    m = max.matching(rke)
    matched.ids <- m$match$pair.id
    CHECK_TRUE(! is.null(matched.ids))
    return(rke.remove.pairs(rke, matched.ids))
  }
}

compute.gi.density <- function(hospital.size, nsims=100) {
  # Computes the pdf of deviation
  gi.density <- NA
  pb = txtProgressBar(style=3)
  for (i in 1:nsims) {
    rke = sample.Rij(hospital.size, Yij=0)
    if (i == 1) {
      gi.density <- summary.Rij(rke)
    } else {
      gi.density <- summary.Rij(rke) + gi.density
    }
    setTxtProgressBar(pb, value=i / nsims)
  }
  gi.density <- gi.density / sum(gi.density)
  gi.density <- list(sims=nsims, size=hospital.size, density=gi.density)
  save(gi.density, file=kGiDensityFile)
}

compute.fi.density <- function(hospital.size, nsims=100) {
  # Computes the pdf of deviation
  fi.density <- NA
  pb = txtProgressBar(style=3)
  for (i in 1:nsims) {
    rke = sample.Rij(hospital.size, Yij=1)
    if (i == 1) {
      fi.density <- summary.Rij(rke)
    } else {
      fi.density <- summary.Rij(rke) + fi.density
    }
    setTxtProgressBar(pb, value=i / nsims)
  }
  fi.density <- fi.density / sum(fi.density)
  fi.density <- list(sims=nsims, size=hospital.size, density=fi.density)
  save(fi.density, file=kFiDensityFile)
}


compute.U.frame <- function(mech, hospital.size, nhospitals,
                            nsims=100, verbose=F) {
  # Creates the data frame U
  # where each row shows U$t = utility of truthful when have U$nt=truthful agents.
  # Nt should be from 0 to nhospitals
  #
  #    U =      t     c     Nt
  #                 ...
  #            8.5   7.5    3
  #            9.5    5.4   4 
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
  ETA = 1  # completion time in minutes
  t0 = proc.time()["elapsed"]
  
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
      
      out$t[nt.index] = out$t[nt.index] + add.ut
      out$c[nt.index] = out$c[nt.index] + add.uc
      count <- count + 1
      setTxtProgressBar(pb, value = count / Nol)
      t1 <- proc.time()["elapsed"]
      mins <- (t1-t0) / 60
      speed <- count / mins  # examples / min
      ETA <- (Nol-count) / speed
      if(verbose) {
        print(sprintf("strategy=%s mech=%s nt=%d ut=%.3f  uc=%.3f i=%d ETA=%.2f mins", 
                      strategy, mech, nt, add.ut, add.uc, i, ETA))
        print(out/nsims)
      }
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


kNDensitySamples <- 1000
kVerboseUF <- T
kNSimsUF <- 200

if(! file.exists(kGiDensityFile)) {
  print("File for Gi density does not exist. Need to compute it.")
  compute.gi.density(hospital.size=kHospitalSize, nsims=kNDensitySamples)
}

if(! file.exists(kFiDensityFile)) {
  print("File for Fi density does not exist. Need to compute it.")
  compute.fi.density(hospital.size=kHospitalSize, nsims=kNDensitySamples)
}

if(! file.exists(kUFrameFileRCM)) {
  print("")
  print("File for rCM U-Frame does not exist. Need to compute it.")
  compute.U.frame(mech="rCM", hospital.size=kHospitalSize, nhospitals=kNoHospitals,
                  nsims=kNSimsUF, verbose=kVerboseUF)
}

if(! file.exists(kUFrameFileXCM)) {
  print("")
  print("File for xCM U-frame does not exist. Need to compute it.")
  compute.U.frame(mech="xCM", hospital.size=kHospitalSize, nhospitals=kNoHospitals,
                  nsims=kNSimsUF, verbose=kVerboseUF)
}

create.template.UFrames <- function() {
  # for rcm
  nr = kNoHospitals + 1
  s = kHospitalSize
  epsilon = 0.3
  Urcm <- data.frame(t=rep(0, nr), c=rep(0, nr))
  Urcm$t <- c(0, rep(0.8 * s/2, nr -1)) + seq(0, kNoHospitals) * epsilon
  Urcm$c <- c(rep(0.87 * s/2, nr - 1), 0) + seq(0, kNoHospitals) * epsilon
  Urcm$c[nr] <- 0
  Urcm$nt <- seq(0, kNoHospitals)
  
  Uxcm <- data.frame(t=rep(0, nr), c=rep(0, nr))
  Uxcm$t <- c(0, rep(0.65 * s/2, nr -1)) + seq(0, kNoHospitals) * epsilon
  Uxcm$c <-  Uxcm$t 
  Uxcm$c[nr] <- 0
  Uxcm$c[1] <- 0.65 * s/2
  Uxcm$nt <- seq(0, kNoHospitals)
  
  return(list(Urcm=Urcm, Uxcm=Uxcm))
}


