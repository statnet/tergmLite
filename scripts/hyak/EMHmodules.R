
### EpiModelHIV Modules to check

new.initialize.hiv <- function(x, param, init, control, s) {

  # browser()
  dat <- list()
  dat$temp <- list()
  if (class(x$fit) == "network") {
    nw <- simulate(x$formation,
                   basis = x$fit,
                   coef = x$coef.form.crude,
                   constraints = x$constraints)
  } else {
    nw <- simulate(x$fit, control = control.simulate.ergm(MCMC.burnin = 1e6))
  }

  if (class(nw)[1] == "networkDynamic") {
    nw <- network.collapse(nw, at = 1)
  }

  dat$el <- as.edgelist(nw)
  attributes(dat$el)$vnames <- NULL
  p <- ergm_prep(nw, x$formation, x$coef.diss$dissolution, x$coef.form,
                 x$coef.diss$coef.adj, x$constraints, control = control.simulate.network())
  p$model.form$formula <- NULL
  p$model.diss$formula <- NULL
  dat$p <- p

  ## Network Model Parameters
  dat$nwparam <- list(x[-which(names(x) == "fit")])

  ## Simulation Parameters
  dat$param <- param
  dat$param$modes <- 1

  dat$init <- init
  dat$control <- control

  ## Nodal Attributes
  dat$attr <- list()

  dat$attr$male <- get.vertex.attribute(nw, "male")

  n <- network.size(nw)
  dat$attr$active <- rep(1, n)
  dat$attr$entTime <- rep(1, n)
  dat$attr$deathTime <- rep(NA, n)
  dat$attr$deathCause <- rep(NA, n)

  ## Initialize HIV related attributes
  dat <- EpiModelHIV:::initStatus(dat)

  dat$attr$age <- get.vertex.attribute(nw, "age")
  dat$attr$agecat <- get.vertex.attribute(nw, "agecat")

  dat <- EpiModelHIV:::initInfTime(dat)
  dat <- EpiModelHIV:::initDx(dat)
  dat <- EpiModelHIV:::initTx(dat)
  dat <- circ(dat, at = 1)

  ## Stats List
  dat$stats <- list()

  ## Final steps
  dat$epi <- list()
  dat <- prevalence.hiv(dat, at = 1)
  dat <- new.simnet.hiv(dat, at = 1)

}


reinit.hiv <- function(x, param, init, control) {
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

  dat$param$modes <- 1
  class(dat) <- "dat"

  return(dat)
}

new.aging.hiv <- function(dat, at) {

  ## Parameters
  time.unit <- dat$param$time.unit
  agecat.cutoff <- dat$param$agecat.cutoff

  ## Attributes
  age <- dat$attr$age
  male <- dat$attr$male
  active <- dat$attr$active
  agecat <- dat$attr$agecat

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  agecat[active == 1 & male == 0 & age >= agecat.cutoff] <- 1
  agecat[active == 1 & male == 1 & age >= agecat.cutoff] <- 3

  ## Save out
  dat$attr$age <- age
  dat$attr$agecat <- agecat

  return(dat)
}


new.deaths.hiv <- function(dat, at) {

  ### 1. Susceptible Deaths ###

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


  ### 2. Infected Deaths ###

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


  ### 3. Update Attributes ###
  if (nDeathsInf > 0) {
    dat$attr$active[idsDeathsInf] <- 0
    dat$attr$deathTime[idsDeathsInf] <- at
    dat$attr$deathCause[idsDeathsInf] <- "i"
  }

  ## 4. Update Population Structure ##
  inactive <- which(dat$attr$active == 0)
  dat$el <- delete_vertices(dat$el, inactive)
  dat$attr <- deleteAttr(dat$attr, inactive)

  if (unique(sapply(dat$attr, length)) != attributes(dat$el)$n) {
    stop("mismatch between el and attr length in death mod")
  }

  ### 5. Summary Statistics ###
  dat$epi$ds.flow[at] <- nDeathsSus
  dat$epi$di.flow[at] <- nDeathsInf

  return(dat)
}


new.births.hiv <- function(dat, at) {

  # Variables
  b.rate.method <- dat$param$b.rate.method
  b.rate <- dat$param$b.rate
  active <- dat$attr$active


  # Process
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


  # Update Population Structure
  if (nBirths > 0) {
    dat <- setBirthAttr(dat, at, nBirths)
    attributes(dat$el)$n <- attributes(dat$el)$n + nBirths
  }

  if (unique(sapply(dat$attr, length)) != attributes(dat$el)$n) {
    stop("mismatch between el and attr length in births mod")
  }

  # Output
  dat$epi$b.flow[at] <- nBirths

  return(dat)
}


update_nwp <- function(dat, at) {

  mf <- dat$p$model.form
  md <- dat$p$model.diss
  mhf <- dat$p$MHproposal.form
  mhd <- dat$p$MHproposal.diss

  n <- attributes(dat$el)$n
  maxdyads <- choose(n, 2)

  ## 1. Update model.form ##

  # edges
  # inputs <- c(0, 1, 0) # not changed
  mf$terms[[1]]$maxval <- maxdyads

  # concurrent by male
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u)
  ui <- seq(along = u)
  inputs <- c(ui, nodecov)

  outputs <- c(0, 2, length(inputs), inputs)
  mf$terms[[2]]$inputs <- outputs
  mf$terms[[2]]$maxval <- n

  # absdiffby
  nodecov <- dat$attr$age
  nodeby <- dat$attr$male
  offsetVal <- mf$terms[[3]]$inputs[4]
  inputs <- c(offsetVal, nodecov, nodeby)
  mf$terms[[3]]$inputs <- c(0, 1, length(inputs), inputs)

  # nodematch
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  inputs <- nodecov
  mf$terms[[4]]$inputs <- c(0, 1, length(inputs), inputs)

  ## Update combined maxval here
  mf$maxval <- c(maxdyads, rep(n, 2), Inf, Inf)


  ## 2. Update model.diss ##
  md$terms[[1]]$maxval <- maxdyads
  md$maxval <- maxdyads


  ## 3. Update MHproposal.form ##
  mhf$arguments$constraints$bd$attribs <-
    matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxout <-
    matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxin <- matrix(rep(n, n), ncol = 1)
  mhf$arguments$constraints$bd$minout <-
    mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)


  ## 4. Update MHproposal.diss ##
  mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd


  ## 5. Output ##
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)
  dat$p <- p

  return(dat)

}

update_nwp_lim <- function(dat, at) {

  mf <- dat$p$model.form
  md <- dat$p$model.diss
  mhf <- dat$p$MHproposal.form
  mhd <- dat$p$MHproposal.diss

  n <- attributes(dat$el)$n
  maxdyads <- choose(n, 2)

  ## 1. Update model.form ##

  # edges
  # inputs <- c(0, 1, 0) # not changed
  mf$terms[[1]]$maxval <- maxdyads

  # nodematch
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  inputs <- nodecov
  mf$terms[[2]]$inputs <- c(0, 1, length(inputs), inputs)

  ## Update combined maxval here
  mf$maxval <- c(maxdyads, Inf)


  ## 2. Update model.diss ##
  md$terms[[1]]$maxval <- maxdyads
  md$maxval <- maxdyads


  ## 3. Update MHproposal.form ##
  mhf$arguments$constraints$bd$attribs <-
    matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxout <-
    matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxin <- matrix(rep(n, n), ncol = 1)
  mhf$arguments$constraints$bd$minout <-
    mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)


  ## 4. Update MHproposal.diss ##
  mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd


  ## 5. Output ##
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)
  dat$p <- p

  return(dat)

}


new.simnet.hiv <- function(dat, at) {

  nwparam <- get_nwparam(dat)
  if (at == 1) {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.crude)
  } else {
    coef.diss <- as.numeric(nwparam$coef.diss$coef.adj)
  }

  dat$el <- simulate_network(p = dat$p,
                             el = dat$el,
                             coef.form = nwparam$coef.form,
                             coef.diss = coef.diss,
                             time.start = at)

  if (at == 1) {
    dat$stats$nwstats <- matrix(NA, ncol = 5, nrow = dat$control$nsteps)
    colnames(dat$stats$nwstats) <- c("edges", "meandeg", "deg0", "deg1", "concurrent")
  }
  n <- attributes(dat$el)$n
  dat$stats$nwstats[at, 1] <- nrow(dat$el)
  dat$stats$nwstats[at, 2] <- nrow(dat$el)/n
  dat$stats$nwstats[at, 4] <- sum(table(dat$el) == 1)/n
  dat$stats$nwstats[at, 5] <- sum(table(dat$el) > 1)/n
  dat$stats$nwstats[at, 3] <- (n - sum(table(dat$el) == 1) - sum(table(dat$el) > 1))/n

  return(dat)
}


new.infect.hiv <- function(dat, at) {

  ## Discordant Edgelist
  del <- new.discord_edgelist.hiv(dat, at)

  nInf <- 0
  idsInf <- idsTrans <- NULL

  if (!is.null(del)) {

    ## Acts
    del <- new.acts(dat, del, at)

    ## Transmission
    dat <- new.trans(dat, del, at)
    del <- dat$del
    dat$del <- NULL

    ## Update Nodal Attr
    idsInf <- unique(del$sus)
    idsTrans <- unique(del$inf)
    nInf <- length(idsInf)

    if (nInf > 0) {
      dat$attr$status[idsInf] <- 1
      dat$attr$infTime[idsInf] <- at
      dat$attr$ageInf[idsInf] <- dat$attr$age[idsInf]
      dat$attr$dxStat[idsInf] <- 0
      dat$attr$vlLevel[idsInf] <- 0
      dat$attr$txCD4min[idsInf] <-
        pmin(rnbinom(nInf,
                     size = nbsdtosize(dat$param$tx.init.cd4.mean,
                                       dat$param$tx.init.cd4.sd),
                     mu = dat$param$tx.init.cd4.mean),
             dat$param$tx.elig.cd4)
    }

    ## Transmission data frame
    if (dat$control$save.transmat == TRUE) {
      if (nInf > 0) {
        if (at == 2) {
          dat$stats$transmat <- as.data.frame(del)
        } else {
          dat$stats$transmat <- rbind(dat$stats$transmat, as.data.frame(del))
        }
      }
    }

  }

  ## Incidence vector
  dat$epi$si.flow[at] <- nInf
  dat$epi$si.flow.male[at] <- sum(dat$attr$male[idsInf] == 1, na.rm = TRUE)
  dat$epi$si.flow.feml[at] <- sum(dat$attr$male[idsInf] == 0, na.rm = TRUE)

  ## Supplemental incidence stats
  if (!is.null(dat$control$getincid.infect)) {
    dat <- do.call(dat$control[["getincid.infect"]], list(dat, at, idsInf, idsTrans))
  }

  return(dat)
}

## within the prevalence module
new.discord_edgelist.hiv <- function(dat, at) {

  status <- dat$attr$status
  active <- dat$attr$active

  idsInft <- which(active == 1 & status == 1)
  nInft <- length(idsInft)

  del <- NULL

  if (nInft > 0) {

    if (is.null(dat$el)) {
      el <- get.dyads.active(dat$nw, at = at)
    } else {
      el <- dat$el
    }

    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]

      disc <- which(abs(status[el[, 1]] - status[el[, 2]]) == 1)
      if (length(disc) > 0) {
        tmp.del <- el[disc, ]
        tmp.del[status[tmp.del[, 2]] == 1, ] <- tmp.del[status[tmp.del[, 2]] == 1, 2:1]

        del <- list()
        del$sus <- tmp.del[, 2]
        del$inf <- tmp.del[, 1]
      }
    }

  }

  return(del)
}


new.acts <- function(dat, del, at) {

  # Variables
  nedges <- length(del[[1]])

  act.rate.early <- dat$param$act.rate.early
  act.rate.late <- dat$param$act.rate.late
  act.rate.cd4 <- dat$param$act.rate.cd4

  cd4Count <- dat$attr$cd4Count[del$inf]

  isLate <- which(cd4Count < act.rate.cd4)

  rates <- rep(act.rate.early, nedges)
  rates[isLate] <- act.rate.late


  # Process
  act.rand <- dat$param$acts.rand
  if (act.rand == TRUE) {
    numActs <- rpois(nedges, rates)
  } else {
    numActs <- rates
  }

  cond.prob <- dat$param$cond.prob
  cond.prob <- rep(cond.prob, nedges)

  del$numActs <- numActs

  if (act.rand == TRUE) {
    del$protActs <- rbinom(nedges, rpois(nedges, numActs), cond.prob)
  } else {
    del$protActs <- numActs * cond.prob
  }

  del$protActs <- pmin(numActs, del$protActs)
  del$unprotActs <- numActs - del$protActs

  stopifnot(all(del$unprotActs >= 0))

  return(del)
}


new.trans <- function(dat, del, at) {

  nedges <- length(del[[1]])
  if (nedges == 0) {
    return(del)
  }

  # Base transmission probability
  vlLevel <- dat$attr$vlLevel[del$inf]
  males <- dat$attr$male[del$sus]
  ages <- dat$attr$age[del$sus]
  circs <- dat$attr$circStat[del$sus]
  prop.male <- dat$epi$propMale[at - 1]
  base.tprob <- hughes_tp(vlLevel, males, ages, circs, prop.male)

  # Acute and aids stage multipliers
  acute.stage.mult <- dat$param$acute.stage.mult
  aids.stage.mult <- dat$param$aids.stage.mult

  isAcute <- which(at - dat$attr$infTime[del$inf] <
                     (dat$param$vl.acute.topeak + dat$param$vl.acute.toset))
  isAIDS <- which(dat$attr$cd4Count[del$inf] < 200)

  base.tprob[isAcute] <- base.tprob[isAcute] * acute.stage.mult
  base.tprob[isAIDS] <- base.tprob[isAIDS] * aids.stage.mult


  # Condoms
  # Probability as a mixture function of protected and unprotected acts
  cond.eff <- dat$param$cond.eff
  prob.stasis.protacts <- (1 - base.tprob*(1 - cond.eff)) ^ del$protActs
  prob.stasis.unptacts <- (1 - base.tprob) ^ del$unprotActs
  prob.stasis <- prob.stasis.protacts * prob.stasis.unptacts
  finl.tprob <- 1 - prob.stasis

  # Transmission
  del$base.tprob <- base.tprob
  del$finl.tprob <- finl.tprob

  stopifnot(length(unique(sapply(del, length))) == 1)

  # Random transmission given final trans prob
  idsTrans <- which(rbinom(nedges, 1, del$finl.tprob) == 1)

  if (at == 2) {
    dat$epi$del.length <- length(del[[1]])
    dat$epi$del.numActs <- mean(del$numActs)
    dat$epi$del.unprotActs <- mean(del$unprotActs)
    dat$epi$del.base.tprob <- mean(del$base.tprob)
    dat$epi$del.finl.tprob <- mean(del$finl.tprob)
  } else {
    dat$epi$del.length[at] <- length(del[[1]])
    dat$epi$del.numActs[at] <- mean(del$numActs)
    dat$epi$del.unprotActs[at] <- mean(del$unprotActs)
    dat$epi$del.base.tprob[at] <- mean(del$base.tprob)
    dat$epi$del.finl.tprob[at] <- mean(del$finl.tprob)
  }

  # Subset discord edgelist to transmissions
  del <- keep.attr(del, idsTrans)

  dat$del <- del

  return(dat)
}

hughes_tp <- function(vls, susmales, susages, suscircs, prop.male, fmat = FALSE) {

  suscircs[is.na(suscircs)] <- 0

  sus.hsv2 <- 0.59*prop.male + 0.86*(1 - prop.male)
  sus.gud <- 0.039*prop.male + 0.053*(1 - prop.male)
  sus.tvagin <- 0.068*prop.male + 0.12*(1 - prop.male)
  sus.cerv <- 0.066*(1 - prop.male)

  interc <- -8.3067
  coef.vl <- 1.062566
  coef.male <- 0.6430989
  coef.age <- -0.0403451
  coef.hsv2 <- 0.7625081
  coef.circ <- -0.6377294
  coef.gud <- 0.9749536
  coef.vagin <- 0.9435334
  coef.cerv <- 1.288279

  tp.full <- exp(interc + coef.vl*(vls - 4) +
                   coef.male*susmales + coef.age*(susages - 35) +
                   coef.hsv2*sus.hsv2 + coef.circ*susmales*suscircs +
                   coef.gud*sus.gud + coef.vagin*sus.tvagin +
                   coef.cerv*sus.cerv)

  if (fmat == TRUE) {
    tp.full <- data.frame(tp.full, vls, susmales, susages, suscircs)
  }

  return(tp.full)
}


nbsdtosize <- function(mu, sd) {
  mu ^ 2 / (sd ^ 2 - mu)
}

get_attr <- function(x, sim = 1) {
  if (is.null(x$attr)) {
    stop("No attr on x")
  } else {
    x$attr[[1]]
  }
}

cut_age <- function(age, breaks = c(0, 29, 39, Inf)) {
  cut(age, breaks = breaks, labels = FALSE)
}

keep.attr <- function(attrList, keep) {
  lapply(attrList, function(x) x[keep])
}
