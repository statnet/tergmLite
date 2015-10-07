
### EpiModelHIV Modules to check

new.initialize.hiv <- function(x, param, init, control, s) {

  # Restarted simulation
  reinit <- ifelse(class(x) == "netsim", TRUE, FALSE)

  if (reinit == FALSE) {
    dat <- list()
    dat$temp <- list()
    if (class(x$fit) == "network") {
      dat$nw <- simulate(x$formation,
                         basis = x$fit,
                         coef = x$coef.form.crude,
                         constraints = x$constraints,
                         control = control.simulate.formula(MCMC.burnin = 1e6))
    } else {
      dat$nw <- simulate(x$fit,
                         control = control.simulate.ergm(MCMC.burnin = 1e6))
    }

    if (class(dat$nw)[1] == "networkDynamic") {
      dat$nw <- network.collapse(dat$nw, at = 1)
    }

    dat$nw <- activate.vertices(dat$nw, onset = 1, terminus = Inf)
    if (control$delete.nodes == TRUE) {
      dat$nw <- network.extract(dat$nw, at = 1)
    }

    ## Network Model Parameters
    dat$nwparam <- list(x[-which(names(x) == "fit")])

    ## Simulation Parameters
    dat$param <- param
    dat$param$modes <- 1

    dat$init <- init
    dat$control <- control

    ## Nodal Attributes
    dat$attr <- list()

    dat$attr$male <- get.vertex.attribute(dat$nw, "male")

    n <- network.size(dat$nw)
    dat$attr$active <- rep(1, n)
    dat$attr$entTime <- rep(1, n)
    dat$attr$deathTime <- rep(NA, n)
    dat$attr$deathCause <- rep(NA, n)

    ## Initialize HIV related attributes
    dat <- initStatus(dat)
    dat <- initAge(dat)
    dat <- initInfTime(dat)
    dat <- initDx(dat)
    dat <- initTx(dat)
    dat <- circ(dat, at = 1)

    ## Stats List
    dat$stats <- list()

    ## Final steps
    dat$epi <- list()
    dat <- prevalence.hiv(dat, at = 1)
    dat <- simnet.hiv(dat, at = 1)

  } else {
    dat <- list()
    dat$nw <- x$network[[s]]
    dat$param <- param
    dat$param$modes <- 1
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- list()
    dat$stats$nwstats <- x$stats$nwstats[[s]]
    dat$temp <- list()
  }

  dat$param$modes <- 1
  class(dat) <- "dat"

  return(dat)
}


new.deaths.hiv <- function(dat, at) {

  # Susceptible Deaths ------------------------------------------------------

  ## Variables
  active <- dat$attr$active
  male <- dat$attr$male
  age <- dat$attr$age
  cd4Count <- dat$attr$cd4Count

  di.cd4.aids <- dat$param$di.cd4.aids
  ds.exit.age <- dat$param$ds.exit.age

  ## Eligible are: active uninf, pre-death infected, unhealthy old
  idsEligSus <- which(active == 1 & (is.na(cd4Count) |
                                       cd4Count > di.cd4.aids |
                                       (cd4Count <= di.cd4.aids & age > ds.exit.age)))
  nEligSus <- length(idsEligSus)

  # Set age-sex specific rates
  ds.rates <- dat$param$ds.rates
  if (nEligSus > 0) {
    rates <- ds.rates$mrate[100*male[idsEligSus] + age[idsEligSus]]
  }


  ## Process
  nDeathsSus <- 0; idsDeathsSus <- NULL
  if (nEligSus > 0) {
    vecDeathsSus <- which(rbinom(nEligSus, 1, rates) == 1)
    nDeathsSus <- length(vecDeathsSus)
  }


  ## Update Attributes
  if (nDeathsSus > 0) {
    idsDeathsSus <- idsEligSus[vecDeathsSus]
    dat$attr$active[idsDeathsSus] <- 0
    dat$attr$deathTime[idsDeathsSus] <- at
    dat$attr$deathCause[idsDeathsSus] <- "s"
  }


  # Infected Deaths ---------------------------------------------------------

  ## Variables
  active <- dat$attr$active
  di.cd4.rate <- dat$param$di.cd4.rate

  ## Process
  nDeathsInf <- 0; idsDeathsInf <- NULL

  cd4Count <- dat$attr$cd4Count
  stopifnot(length(active) == length(cd4Count))

  idsEligInf <- which(active == 1 & cd4Count <= di.cd4.aids)
  nEligInf <- length(idsEligInf)

  if (nEligInf > 0) {
    vecDeathsInf <- which(rbinom(nEligInf, 1, di.cd4.rate) == 1)
    if (length(vecDeathsInf) > 0) {
      idsDeathsInf <- idsEligInf[vecDeathsInf]
      nDeathsInf <- length(idsDeathsInf)
    }
  }

  idsDeathsDet <- which(active == 1 & cd4Count <= 0)
  if (length(idsDeathsDet) > 0) {
    idsDeathsInf <- c(idsDeathsInf, idsDeathsDet)
    nDeathsInf <- nDeathsInf + length(idsDeathsDet)
  }


  ## Update Attributes
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
    dat$attr$deathTime[idsDeathsInf] <- at
    dat$attr$deathCause[idsDeathsInf] <- "i"
  }


  # Update Network ----------------------------------------------------------
  idsDeaths <- c(idsDeathsSus, idsDeathsInf)
  if (length(idsDeaths) > 0) {
    dat$nw <- networkDynamic::deactivate.vertices(dat$nw, onset = at, terminus = Inf,
                                                  v = idsDeaths, deactivate.edges = TRUE)
  }


  # Output ------------------------------------------------------------------
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}


new.births.hiv <- function(dat, at) {

  # Variables ---------------------------------------------------------------
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active
  currNwSize <- network::network.size(dat$nw)


  # Process -----------------------------------------------------------------
  nBirths <- 0
  if (b.rate.method == "stgrowth") {
    exptPopSize <- dat$epi$num[1] * (1 + b.rate*at)
    numNeeded <- exptPopSize - sum(active == 1)
    if (numNeeded > 0) {
      nBirths <- rpois(1, numNeeded)
    }
  }
  if (b.rate.method == "totpop") {
    nElig <- dat$epi$num[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }
  if (b.rate.method == "fpop") {
    nElig <- dat$epi$num.feml[at - 1]
    if (nElig > 0) {
      nBirths <- rpois(1, nElig * b.rate)
    }
  }


  # Update Attr -------------------------------------------------------------
  if (nBirths > 0) {
    dat <- setBirthAttr(dat, at, nBirths)
  }


  # Update Network ----------------------------------------------------------
  if (nBirths > 0) {
    newIds <- (currNwSize + 1):(currNwSize + nBirths)

    stopifnot(unique(sapply(dat$attr, length)) == (currNwSize + nBirths))

    dat$nw <- networkDynamic::add.vertices.active(x = dat$nw, nv = nBirths,
                                                  onset = at, terminus = Inf)

    dat$nw <- network::set.vertex.attribute(x = dat$nw,
                                            attrname = c("male", "age", "agecat"),
                                            value = list(male = dat$attr$male[newIds],
                                                         age = dat$attr$age[newIds],
                                                         agecat = dat$attr$agecat[newIds]),
                                            v = newIds)

  }


  # Output ------------------------------------------------------------------
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}


new.simnet.hiv <- function(dat, at) {

  resim.int <- dat$control$resim.int
  if (at > 1 & at %% resim.int > 0) {
    return(dat)
  }

  # Delete Nodes ------------------------------------------------------------
  if (at > 1 & dat$control$delete.nodes == TRUE) {
    dat$nw <- network.extract(dat$nw, at = at)
    inactive <- which(dat$attr$active == 0)
    dat$attr <- deleteAttr(dat$attr, inactive)
  }

  # Resimulation ------------------------------------------------------------
  nwparam <- get_nwparam(dat)
  if (at == 1) {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.crude)
  } else {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.adj)
  }

  suppressWarnings(
    dat$nw <- simulate(
      dat$nw,
      formation = nwparam$formation,
      dissolution = nwparam$dissolution,
      coef.form = as.numeric(nwparam$coef.form),
      coef.diss = coef.diss,
      constraints = nwparam$constraints,
      time.start = at,
      time.slices = 1 * resim.int,
      time.offset = 0,
      output = "networkDynamic",
      monitor = dat$control$nwstats.formula))

  if (at == 1) {
    dat$stats$nwstats <- as.data.frame(attributes(dat$nw)$stats)
  } else {
    dat$stats$nwstats <- rbind(dat$stats$nwstats,
                               tail(attributes(dat$nw)$stats, 1 * resim.int))
  }


  return(dat)
}

## within the prevalence module
discord_edgelist.hiv <- function(dat, at) {

  status <- dat$attr$status
  active <- dat$attr$active

  idsInft <- which(active == 1 & status == "i")
  nInft <- length(idsInft)

  del <- NULL

  if (nInft > 0) {

    el <- get.dyads.active(dat$nw, at = at)
    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]
      stat <- matrix(status[el], ncol = 2)
      isInf <- matrix(stat == "i", ncol = 2)
      isSus <- matrix(stat == "s", ncol = 2)
      SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
      ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
      pairs <- rbind(SIpairs, ISpairs[,2:1])

      del <- list()
      del$at <- rep(at, nrow(pairs))
      del$sus <- numeric()
      del$inf <- numeric()
      if (nrow(pairs) > 0) {
        del$sus <- pairs[, 1]
        del$inf <- pairs[, 2]
      }
    }

  }

  return(del)
}
