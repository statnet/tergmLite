
#' @export
simulate_network <- function(nw = NULL,
                             p = NULL,
                             el,
                             formation = NULL,
                             dissolution = NULL,
                             coef.form = NULL,
                             coef.diss = NULL,
                             constraints = ~.,
                             time.slices = 1,
                             time.start = NULL,
                             time.burnin = 0,
                             time.interval = 1,
                             time.offset = 0,
                             control = control.simulate.network(),
                             output = "network") {

  control$changes <- TRUE

  if (!is.null(nw) & is.null(p)) {
    p <- ergm_prep(nw, formation, dissolution, coef.form, coef.diss, constraints, control)
  }

  if (!is.null(el)) {
    n <- attributes(el)$n
  }

  eta.form <- ergm.eta(coef.form, p$model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, p$model.diss$etamap)
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- FALSE
  control$collect.diss <- FALSE

  z <- stergm_getMCMCsample(nw, el,
                            p$model.form, p$model.diss,
                            p$MHproposal.form, p$MHproposal.diss,
                            eta.form, eta.diss,
                            control, output)

  if (output == "network") {
    out <- z$newnetwork
  }
  if (output == "edgelist") {
    out <- z
    attributes(out)$n <- n
  }

  return(out)
}


#' @export
stergm_getMCMCsample <- function(nw = NULL, el = NULL,
                                 model.form, model.diss,
                                 MHproposal.form, MHproposal.diss,
                                 eta.form, eta.diss,
                                 control, output) {

  verbose <- FALSE
  model.mon <- NULL

  if (output == "edgelist") {
    Clist.form <- ergm_Cprepare(el = el, m = model.form)
    Clist.diss <- ergm_Cprepare(el = el, m = model.diss)
  }
  if (output == "network") {
    Clist.form <- ergm_Cprepare(nw = nw, m = model.form)
    Clist.diss <- ergm_Cprepare(nw = nw, m = model.diss)
  }


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

  out <- newnw_extract(nw, z, output)

  return(out)
}


#' @export
newnw_extract <- function(oldnw = NULL, z, output = "network") {

  nedges <- z$newnwtails[1]
  newedgelist <- if (nedges > 0) {
    cbind(z$newnwtails[2:(nedges + 1)], z$newnwheads[2:(nedges + 1)])
  } else {
    matrix(0, ncol = 2, nrow = 0)
  }

  if (output == "edgelist") {
    newnw <- newedgelist
  } else {
    newnw <- network.initialize(oldnw$gal$n, directed = FALSE)
    newnw <- add.edges(newnw, tail = newedgelist[, 1], head = newedgelist[, 2])
  }

  return(newnw)
}
