
#' @export
simulate_network <- function(nw,
                             el,
                             formation,
                             dissolution,
                             coef.form,
                             coef.diss,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1,
                             time.start = NULL,
                             time.burnin = 0,
                             time.interval = 1,
                             time.offset = 0,
                             control = control.simulate.network(),
                             output = "network") {

  control$changes <- TRUE

  if (!is.network(nw)) {
    stop("A network object must be given")
  }

  formation <- ergm.update.formula(formation, nw ~ ., from.new = "nw")          ## NW
  dissolution <- ergm.update.formula(dissolution, nw ~ ., from.new = "nw")

  if (is.character(monitor)) {
    monitor <- formation
  }

  if (!is.null(monitor)) {
    monitor <- ergm.update.formula(monitor, nw ~ ., from.new = "nw")            ## NW
  }

  nw$gal$lasttoggle <- NULL

  model.form <- ergm_getmodel(formation, nw, role = "formation")                ## NW
  if (!missing(coef.form) && coef.length.model(model.form) != length(coef.form)) {
    stop("coef.form has ", length(coef.form), " elements, while the model requires ",
         coef.length.model(model.form), " parameters.")
  }
  model.diss <- ergm_getmodel(dissolution, nw, role = "dissolution")            ## NW
  if (!missing(coef.diss) && coef.length.model(model.diss) != length(coef.diss)) {
    stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",
         coef.length.model(model.diss), " parameters.")
  }

  model.mon <- if (!is.null(monitor)) {
    ergm.getmodel(monitor, nw, role = "target")                                 ## NW
  } else {
    NULL
  }


  ####### Model simulation #######

  MHproposal.form <- MHproposal(constraints,                                    ## NW
                                control$MCMC.prop.args.form,
                                nw,
                                weights = control$MCMC.prop.weights.form,
                                class = "f")
  MHproposal.diss <- MHproposal(constraints,
                                control$MCMC.prop.args.diss,
                                nw,
                                weights = control$MCMC.prop.weights.diss,
                                class = "d")

  eta.form <- ergm.eta(coef.form, model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, model.diss$etamap)
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- FALSE
  control$collect.diss <- FALSE

  z <- stergm_getMCMCsample(nw,
                            el,
                            model.form,
                            model.diss,
                            model.mon,
                            MHproposal.form,
                            MHproposal.diss,
                            eta.form,
                            eta.diss,
                            control,
                            output)

  if (output == "network") {
    out <- z$newnetwork
  }
  if (output == "edgelist") {
    out <- z
  }

  return(out)
}


#' @export
stergm_getMCMCsample <- function(nw,
                                 el,
                                 model.form,
                                 model.diss,
                                 model.mon,
                                 MHproposal.form,
                                 MHproposal.diss,
                                 eta.form,
                                 eta.diss,
                                 control,
                                 output) {
# browser()
  verbose <- FALSE

  # summary(microbenchmark(ergm.Cprepare(nw, model.form)), unit = "s") # med = 0.044 s
  # summary(microbenchmark(ergm.Cprepare(nw, model.diss)), unit = "s") # med = 0.043 s

  Clist.form <- ergm_Cprepare(nw, el, model.form)
  Clist.diss <- ergm_Cprepare(nw, el, model.diss)

  # str(Clist.form)

  Clist.mon <- NULL
  collect.form <- control$collect.form
  collect.diss <- control$collect.diss

  maxedges <- control$MCMC.init.maxedges
  maxchanges <- control$MCMC.init.maxchanges

  getz <- function() {
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
        message("Too many edges encountered in the simulation. Increasing capacity to ",
                maxedges)
      }
      if (z$status == 3) {
        maxchanges <- 5 * maxchanges
        message("Too many changes elapsed in the simulation. Increasing capacity to ",
                maxchanges)
      }
    }
    return(z)
  }

  # summary(microbenchmark(getz()), unit = "s") # med = 0.016 s
  z <- getz()

  # summary(microbenchmark(newnw_extract(nw, z)), unit = "s") # med = 0.141
  newnetwork <- newnw_extract(nw, z, output)

  return(newnetwork)
}


#' @export
newnw_extract <- function(oldnw, z, output = "network", response = NULL) {

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
