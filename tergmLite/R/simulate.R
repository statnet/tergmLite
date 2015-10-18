
#' @export
simulate_network <- function(p, el, coef.form, coef.diss, save.changes = FALSE) {

  control <- tergm::control.simulate.network()
  control$changes <- TRUE

  eta.form <- ergm.eta(coef.form, p$model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, p$model.diss$etamap)

  control$time.burnin <- 0
  control$time.interval <- 1
  control$time.samplesize <- 1
  control$collect.form <- FALSE
  control$collect.diss <- FALSE

  z <- stergm_getMCMCsample(el, p$model.form, p$model.diss,
                            p$MHproposal.form, p$MHproposal.diss,
                            eta.form, eta.diss, control, save.changes)

  attributes(z)$n <- attributes(el)$n

  return(z)
}


#' @export
stergm_getMCMCsample <- function(el, model.form, model.diss,
                                 MHproposal.form, MHproposal.diss,
                                 eta.form, eta.diss, control, save.changes) {

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
  out <- if (nedges > 0) {
    cbind(z$newnwtails[2:(nedges + 1)], z$newnwheads[2:(nedges + 1)])
  } else {
    matrix(0, ncol = 2, nrow = 0)
  }

  if (save.changes == TRUE) {
    if (z$diffnwtime[1] > 0) {
      changes <- cbind(z$diffnwtails[2:(z$diffnwtails[1] + 1)],
                       z$diffnwheads[2:(z$diffnwheads[1] + 1)],
                       z$diffnwdirs[2:(z$diffnwdirs[1] + 1)])

    } else {
      changes <- matrix(0, ncol = 3, nrow = 0)
    }
    colnames(changes) <- c("tail", "head", "to")
    attributes(out)$changes <- changes
  }

  return(out)
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


#' @export
simulate_ergm_formula <- function(object,
                                  nsim = 1,
                                  seed = NULL,
                                  coef,
                                  response = NULL,
                                  reference = ~Bernoulli,
                                  constraints = ~.,
                                  monitor = NULL,
                                  basis = NULL,
                                  statsonly = FALSE,
                                  esteq = FALSE,
                                  sequential = TRUE,
                                  control = control.simulate.formula(),
                                  verbose = FALSE, ...) {

  check.control.class(myname = "ERGM simulate.formula")

  control <- control.simulate.ergm.toplevel(control, ...)
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  if (is.null(nw <- basis)) {
    nw <- ergm.getnetwork(object)
  }

  nw <- as.network(nw)
  if (!is.network(nw)) {
    stop("A network object on the LHS of the formula or via",
         " the 'basis' argument must be given")
  }

  if (is.null(basis)) {
    basis <- nw
  }
  form <- ergm.update.formula(object, basis ~ ., from.new = "basis")
  if (!is.null(monitor)) {
    monitor <- ergm.update.formula(monitor, nw ~ ., from.new = "nw")
    monitor.m <- ergm.getmodel(monitor, basis, response = response)
    monitored.length <- coef.length.model(monitor.m)
    monitor <- term.list.formula(monitor[[3]])
    form <- append.rhs.formula(form, monitor)
  }
  else {
    monitored.length <- 0
  }

  m <- ergm.getmodel(form, basis, response = response, role = "static")
  if (missing(coef)) {
    coef <- c(rep(0, coef.length.model(m)))
    warning("No parameter values given, using Bernouli network\n\t")
  }

  coef <- c(coef, rep(0, monitored.length))
  if (coef.length.model(m) != length(coef)) {
    stop("coef has ", length(coef) - monitored.length, " elements, while the model requires ",
         coef.length.model(m) - monitored.length, " parameters.")
  }

  MHproposal <- MHproposal(constraints, arguments = control$MCMC.prop.args,
                           nw = nw, weights = control$MCMC.prop.weights, class = "c",
                           reference = reference, response = response)
  if (any(is.nan(coef) | is.na(coef))) {
    stop("Illegal value of coef passed to simulate.formula")
  }

  eta0 <- ergm.eta(coef, m$etamap)
  curstats <- summary(form, response = response)
  names(curstats) <- m$coef.names
  control$MCMC.init.maxedges <- 1 + max(control$MCMC.init.maxedges,
                                        network.edgecount(nw))

  if (verbose) {
    cat(paste("Starting MCMC iterations to generate ", nsim,
              " network", ifelse(nsim > 1, "s\n", "\n"), sep = ""))
  }

  nthreads <- max(if (inherits(control$parallel, "cluster"))
    nrow(summary(control$parallel)) else control$parallel, 1)

  if (sequential && statsonly) {
    control$MCMC.samplesize <- nsim
    z <- ergm.getMCMCsample(nw, m, MHproposal, eta0, control,
                            verbose = verbose, response = response)
    colnames(z$statsmatrix) <- m$coef.names
    out.mat <- sweep(z$statsmatrix[seq_len(nsim), , drop = FALSE], 2, curstats, "+")
  } else {
    if (!statsonly) {
      nw.list <- list()
    }
    out.mat <- matrix(nrow = 0, ncol = length(curstats),
                      dimnames = list(NULL, m$coef.names))
    if (nthreads > 1)
      curstats <- matrix(curstats, nrow = nthreads, ncol = length(curstats),
                         byrow = TRUE)
    for (i in 1:ceiling(nsim/nthreads)) {
      control$MCMC.samplesize <- nthreads
      control$MCMC.burnin <- if (i == 1 || sequential == FALSE) {
        control$MCMC.burnin
      } else {
        control$MCMC.interval
      }
      z <- ergm.getMCMCsample(nw, m, MHproposal, eta0,
                              control, verbose = verbose, response = response)
      out.mat <- rbind(out.mat, curstats + z$statsmatrix)
      if (!statsonly)
        if (nthreads > 1)
          nw.list[[length(nw.list) + 1]] <- z$newnetwork
      else nw.list <- c(nw.list, z$newnetworks)
      if (sequential) {
        nw <- if (nthreads > 1)
          z$newnetwork
        else z$newnetworks
        curstats <- curstats + z$statsmatrix
      }
      if (verbose) {
        cat(sprintf("Finished simulation %d of %d.\n", i, nsim))
      }
    }
  }
  out.mat <- out.mat[seq_len(nsim), , drop = FALSE]
  if (esteq)
    out.mat <- .ergm.esteq(coef, m, out.mat)
  if (statsonly)
    return(out.mat)
  if (nsim == 1) {
    return(nw.list[[1]])
  } else {
    nw.list <- nw.list[seq_len(nsim)]
    attributes(nw.list) <- list(formula = object, stats = out.mat,
                                coef = coef, control = control, constraints = constraints,
                                reference = reference, monitor = monitor, response = response)
    class(nw.list) <- "network.list"
    return(nw.list)
  }
}
