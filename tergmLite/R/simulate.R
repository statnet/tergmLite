
#' @export
simulate_network <- function(p,
                             el,
                             coef.form,
                             coef.diss,
                             time.slices = 1,
                             time.start = NULL,
                             time.burnin = 0,
                             time.interval = 1,
                             time.offset = 0,
                             control = control.simulate.network()) {

  control$changes <- TRUE

  n <- attributes(el)$n

  eta.form <- ergm.eta(coef.form, p$model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, p$model.diss$etamap)
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- FALSE
  control$collect.diss <- FALSE

  z <- stergm_getMCMCsample(el, p$model.form, p$model.diss,
                            p$MHproposal.form, p$MHproposal.diss,
                            eta.form, eta.diss, control)

  out <- z[[1]]
  attributes(out)$n <- n
  attributes(out)$changes <- z[[2]]

  return(out)
}


#' @export
stergm_getMCMCsample <- function(el, model.form, model.diss,
                                 MHproposal.form, MHproposal.diss,
                                 eta.form, eta.diss, control) {

  verbose <- FALSE
  model.mon <- NULL

  Clist.form <- ergm_Cprepare(el = el, m = model.form)
  Clist.diss <- ergm_Cprepare(el = el, m = model.diss)

  Clist.mon <- NULL
  collect.form <- control$collect.form
  collect.diss <- control$collect.diss

  maxedges <- control$MCMC.init.maxedges
  maxchanges <- control$MCMC.init.maxchanges

  repeat {
    z <- .C("MCMCDyn_wrapper",
            as.integer(Clist.form$tails),
            as.integer(Clist.form$heads),
            time = if (is.null(Clist.form$time)) as.integer(0)
                       else as.integer(Clist.form$time),
            lasttoggle = as.integer(NVL(Clist.form$lasttoggle,
                                        Clist.diss$lasttoggle, Clist.mon$lasttoggle, 0)),
            as.integer(Clist.form$nedges),
            as.integer(Clist.form$n),
            as.integer(Clist.form$dir),
            as.integer(Clist.form$bipartite),
            as.integer(Clist.form$nterms),
            as.character(Clist.form$fnamestring),
            as.character(Clist.form$snamestring),
            as.character(MHproposal.form$name),
            as.character(MHproposal.form$pkgname),
            as.double(Clist.form$inputs),
            as.double(ergm:::.deinf(eta.form)),
            as.integer(Clist.diss$nterms),
            as.character(Clist.diss$fnamestring),
            as.character(Clist.diss$snamestring),
            as.character(MHproposal.diss$name),
            as.character(MHproposal.diss$pkgname),
            as.double(Clist.diss$inputs),
            as.double(ergm:::.deinf(eta.diss)),
            if (!is.null(model.mon)) as.integer(Clist.mon$nterms) else as.integer(0),
            if (!is.null(model.mon)) as.character(Clist.mon$fnamestring) else character(0),
            if (!is.null(model.mon)) as.character(Clist.mon$snamestring) else character(0),
            if (!is.null(model.mon)) as.double(Clist.mon$inputs) else double(0),
            as.integer(MHproposal.form$arguments$constraints$bd$attribs),
            as.integer(MHproposal.form$arguments$constraints$bd$maxout),
            as.integer(MHproposal.form$arguments$constraints$bd$maxin),
            as.integer(MHproposal.form$arguments$constraints$bd$minout),
            as.integer(MHproposal.form$arguments$constraints$bd$minin),
            as.integer(MHproposal.form$arguments$constraints$bd$condAllDegExact),
            as.integer(length(MHproposal.form$arguments$constraints$bd$attribs)),
            as.integer(control$time.samplesize),
            as.integer(control$MCMC.burnin.min),
            as.integer(control$MCMC.burnin.max),
            as.double(control$MCMC.burnin.pval),
            as.double(control$MCMC.burnin.add),
            as.integer(control$time.burnin),
            as.integer(control$time.interval),
            collect.form = as.integer(collect.form),
            s.form = if (collect.form) double(Clist.form$nstats *
                                                (control$time.samplesize + 1))
                       else double(0),
            collect.diss = as.integer(collect.diss),
            s.diss = if (collect.diss) double(Clist.diss$nstats *
                                                (control$time.samplesize + 1))
                       else double(0),
            s.mon = if (!is.null(model.mon)) double(Clist.mon$nstats *
                                                      (control$time.samplesize + 1))
                       else double(0),
            as.integer(maxedges),
            newnwtails = integer(maxchanges),
            newnwheads = integer(maxchanges),
            as.integer(maxchanges),
            as.integer(control$changes),
            diffnwtime = if (control$changes) integer(maxchanges) else integer(0),
            diffnwtails = if (control$changes) integer(maxchanges) else integer(0),
            diffnwheads = if (control$changes) integer(maxchanges) else integer(0),
            diffnwdirs = if (control$changes) integer(maxchanges) else integer(0),
            as.integer(verbose),
            status = integer(1),
            PACKAGE = "tergm")
    if (z$status == 0) {
      break
    }
    if (z$status == 1) {
      maxedges <- 5 * maxedges
      message("Too many edges encountered. Increasing capacity to ", maxedges)
    }
    if (z$status == 3) {
      maxchanges <- 5 * maxchanges
      message("Too many changes elapsed. Increasing capacity to ", maxchanges)
    }
  }

  nedges <- z$newnwtails[1]
  newedgelist <- if (nedges > 0) {
    cbind(z$newnwtails[2:(nedges + 1)], z$newnwheads[2:(nedges + 1)])
  } else {
    matrix(0, ncol = 2, nrow = 0)
  }

  return(newedgelist)
}


ergm_Cprepare <- function(el, m, response = NULL) {

  n <- attributes(el)$n

  dir <- FALSE

  Clist <- list(n = n, dir = dir)
  bip <- FALSE
  if (is.null(bip)) {
    bip <- 0
  }
  Clist$bipartite <- bip
  Clist$ndyads <- choose(n, 2)

  e <- el

  if (length(e) == 0) {
    Clist$nedges <- 0
    Clist$tails <- NULL
    Clist$heads <- NULL
    if (!is.null(response))
      Clist$weights <- numeric(0)
  } else {
    if (!is.matrix(e)) {
      e <- matrix(e, ncol = 2 + (!is.null(response)))
    }
    if (!is.null(response)) {
      e <- e[e[, 3] != 0, , drop = FALSE]
    }

    Clist$nedges <- dim(e)[1]
    Clist$tails <- e[, 1]
    Clist$heads <- e[, 2]
    if (!is.null(response)) {
      Clist$weights <- e[, 3]
    }

  }

  Clist$lasttoggle <- NULL
  Clist$time <- NULL
  mo <- m$terms
  Clist$nterms <- length(mo)
  Clist$nstats <- 0
  Clist$fnamestring <- ""
  Clist$snamestring <- ""
  Clist$inputs <- numeric(0)
  if (Clist$nterms > 0) {
    for (i in 1:Clist$nterms) {
      term_i <- mo[[i]]
      Clist$fnamestring <- paste(Clist$fnamestring, term_i$name)
      Clist$snamestring <- paste(Clist$snamestring, if (!is.null(term_i$soname)) {
        term_i$soname
      }
      else if (!is.null(term_i$pkgname)) {
        term_i$pkgname
      }
      else stop("ERGM term specifying C function `", term_i$name,
                "' is missing C library or package name."))
      Clist$inputs <- c(Clist$inputs, term_i$inputs)
      Clist$nstats <- Clist$nstats + term_i$inputs[2]
    }
  }
  while (substring(Clist$fnamestring, 1, 1) == " ") {
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  }
  while (substring(Clist$snamestring, 1, 1) == " ") {
    Clist$snamestring <- substring(Clist$snamestring, 2)
  }

  Clist$diagnosable <- !m$etamap$offsetmap
  names(Clist$diagnosable) <- m$coef.names

  return(Clist)
}
