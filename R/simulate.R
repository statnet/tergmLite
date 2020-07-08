
# STERGM --------------------------------------------------------------

#' @title Fast Version of tergm::simulate.network for networkLite Object
#'
#' @description Resimulates a networkLite object given thin network data structure,
#'              edgelist, and STERGM model coefficients.
#'
#' @param p A list of network-related nodal covariates and related terms that
#'          is produced with \code{\link{stergm_prep}}.
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size.
#' @param coef.form Vector of coefficients associated with the formation formula.
#' @param coef.diss Vector of coefficients associated with the dissolution formula.
#' @param save.changes Logical, if \code{TRUE}, saves a matrix of changed edges
#'                     as an attribute of the output edgelist matrix.
#'
#' @details
#' This function is used within the network resimulation module in \code{EpiModel}
#' to update temporal ERGMs based on the model coefficients and current network
#' structure. If network structure (e.g., number of nodes) or nodal attributes
#' has changed since the last simulation, this network resimulation should be run
#' only after \code{\link{updateModelTermInputs}}.
#'
#' @export
#'
#' @return
#' Returns an updated network edgelist object, typically stored on the master
#' \code{dat} list object, based on the model simulation. If \code{save.changes}
#' is \code{TRUE}, also returns a list of new edges and disolved edges with the
#' resimulation.
#'
#' @examples
#' library("EpiModel")
#'
#' # Set seed for reproducibility
#' set.seed(1234)
#'
#' nw <- network_initialize(100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges + nodefactor("group")
#' target.stats <- c(15, 10)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.25)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # Full network structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' # Current network structure
#' dat$el[[1]]
#'
#' # New network structure
#' dat$el[[1]] <- simulate_network(p = dat$p[[1]],
#'                                 el = dat$el[[1]],
#'                                 coef.form = dat$nwparam[[1]]$coef.form,
#'                                 coef.diss = dat$nwparam[[1]]$coef.diss$coef.adj,
#'                                 save.changes = TRUE)
#' dat$el[[1]]
#'
#' # Specific changes listed under changes list
#' #    (new edges: to = 1; dissolved edges: to = 0):
#' attributes(dat$el[[1]])$changes
#'
simulate_network <- function(p,
                             el,
                             coef.form,
                             coef.diss,
                             save.changes = FALSE) {

  control <- control.simulate.network()
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
  # copy other attributes of the input el (size, directedness, etc) to the output
  attributes(z)$n <- attributes(el)$n
  attributes(z)$directed <- attributes(el)$directed
  attributes(z)$bipartite <- attributes(el)$bipartite
  attributes(z)$loops <- attributes(el)$loops
  attributes(z)$inverted <- attributes(el)$inverted
  class(z) <- c("edgelist", class(z))

  return(z)
}

stergm_getMCMCsample <- function(el,
                                 model.form,
                                 model.diss,
                                 MHproposal.form,
                                 MHproposal.diss,
                                 eta.form,
                                 eta.diss,
                                 control,
                                 save.changes) {

  Clist.form <- ergm_Cprepare(el = el, m = model.form)
  Clist.diss <- ergm_Cprepare(el = el, m = model.diss)
  Clist.mon <- NULL

  control.parallel <- control
  control.parallel$MCMC.samplesize <- 1

  z <- stergm_MCMC_slave(Clist.form = Clist.form, Clist.diss = Clist.diss,
                         Clist.mon = Clist.mon,
                         proposal.form = MHproposal.form,
                         proposal.diss = MHproposal.diss,
                         eta.form = eta.form, eta.diss = eta.diss,
                         control = control.parallel, verbose = FALSE)

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

ergm_Cprepare <- function(el,
                          m,
                          response = NULL) {

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


# ERGM ----------------------------------------------------------------

#' @title Fast Version of ergm::simulate.ergm for Edgelist-formatted Network
#'
#' @description Resimulates a networkLite object given thin network data structure,
#'              edgelist, and ERGM model coefficients.
#'
#' @param p A list of network-related nodal covariates and related terms that
#'          is produced with \code{\link{ergm_prep}}.
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size.
#' @param coef Vector of coefficients associated with the formation formula.
#'
#' @details
#' This function is used within the network resimulation module in \code{EpiModel}
#' to update cross-sectional ERGMs based on the model coefficients and current
#' network structure. If network structure (e.g., number of nodes) or nodal attributes
#' has changed since the last simulation, this network resimulation should be run
#' only after \code{\link{updateModelTermInputs}}.
#'
#' @export
#'
#' @return
#' Returns an updated network edgelist object, typically stored on the master
#' \code{dat} list object, based on the model simulation.
#'
#' @examples
#' library("EpiModel")
#'
#' # Set seed for reproducibility
#' set.seed(1234)
#'
#' nw <- network_initialize(100)
#' nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
#' formation <- ~edges + nodefactor("group")
#' target.stats <- c(15, 10)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.1)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # Full network structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' # Current network structure
#' dat$el[[1]]
#'
#' # New network structure (all edges are new)
#' dat$el[[1]] <- simulate_ergm(p = dat$p[[1]],
#'                              el = dat$el[[1]],
#'                              coef = dat$nwparam[[1]]$coef.form)
#' dat$el[[1]]
#'
simulate_ergm <- function(p, el, coef) {

  control <- control.simulate.formula()

  eta0 <- ergm.eta(coef, p$model.form$etamap)

  control$MCMC.samplesize <- 1

  z <- ergm_getMCMCsample(el, p$model.form, p$MHproposal, eta0, control)

  attributes(z)$n <- attributes(el)$n

  return(z)
}


ergm_getMCMCsample <- function(el, model, MHproposal, eta0, control) {

  Clist <- ergm_Cprepare(el, model, response = NULL)
  control.parallel <- control
  control.parallel$MCMC.samplesize <- 1

  z <- ergm_MCMC_slave(Clist = Clist, prev.run = NULL,
                      burnin = NULL, samplesize = NULL, interval = NULL,
                      maxedges = NULL, proposal = MHproposal, eta = eta0,
                      control = control.parallel, verbose = FALSE)

  nedges <- z$newnwtails[1]
  el <- if (nedges > 0) {
    cbind(z$newnwtails[2:(nedges + 1)], z$newnwheads[2:(nedges + 1)])
  } else {
    matrix(0, ncol = 2, nrow = 0)
  }

  return(el)
}
