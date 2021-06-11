
#' @title Initializes EpiModel netsim Object for tergmLite Simulation
#'
#' @param dat A list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{netsim}.
#'
#' @details
#' This function is typically used within the initialization modules of
#' \code{EpiModel} to establish the necessary \code{networkLite} infrastructure
#' needed for \code{tergmLite} network resimulation. Specifically, this function
#' converts (and then removes) the \code{network} class objects into an edgelist
#' only format and prepares the ERGM structural information for simulation. The
#' example below demonstrates the specific information returned.
#'
#' Implemented terms are:
#' \itemize{
#'   \item edges
#'   \item nodematch
#'   \item nodefactor
#'   \item concurrent (including heterogenous by attribute)
#'   \item degree (including heterogenous by attribute)
#'   \item degrange
#'   \item absdiff
#'   \item absdiffby (in the EpiModel package)
#'   \item nodecov
#'   \item nodemix
#'   \item absdiffnodemix (in the EpiModel package)
#'   \item triangle
#'   \item gwesp(fixed=TRUE)
#'   \item mean.age
#' }
#' All other terms will return errors.
#'
#' @export
#'
#' @return
#' Returns the list object \code{dat} and adds two elements to the objects: \code{el}
#' is an edgelist representation of the network; and \code{p} is a list object
#' that contains all the relevant structural information for ERGM/TERGM simulation.
#' The function also removes the network class object on the \code{dat} object,
#' stored under \code{nw} because it is no longer needed.
#'
#' @examples
#' \dontrun{
#' library("EpiModel")
#' nw <- network_initialize(100)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
#'
#' # networkLite representation after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' # networkLite representation used by tergmLite
#' str(dat$p, max.level = 3)
#'
#' # Elements removed are nw (network class object)
#' # Elements added are el (edgelist representation of network)...
#' dat$el
#'
#' # ... and p (contains all relevant ERGM structural information for simulation)
#' str(dat$p, max.level = 3)
#' }

#'
init_tergmLite <- function(dat) {

  num_nw <- ifelse(inherits(dat$nw, "network"), 1, length(dat$nw))

  dat$el <- list()
  dat$p <- list()
  dat$control$mcmc.control <- list()
  dat$control$nwstats.formulas <- list()

  supported.terms <- c("edges", "nodematch", "nodefactor",
                       "concurrent", "concurrent_by_attr",
                       "degree", "degree_by_attr", "meandeg",
                       "absdiff", "absdiffby", "nodecov", "nodemix",
                       "absdiffnodemix", "degrange", "triangle", "gwesp",
                       "mean_age", "_lasttoggle", "_utp_wtnet",
                       "on_union_lt_net_Network", "on_intersect_lt_net_Network",
                       "_union_lt_net_Network", "_previous_lt_net_Network",
                       "_intersect_lt_net_Network")

  for (i in 1:num_nw) {
    dat$p[[i]] <- list()

    nwp <- dat$nwparam[[i]]
    is_tergm <- all(nwp$coef.diss$duration > 1)

    nw <- dat$nw[[i]]

    dat$el[[i]] <- as.edgelist(nw)
    attributes(dat$el[[i]])$vnames <- NULL

    nwstats_formula_name <- paste(c("nwstats.formula", if (num_nw > 1) i), collapse = ".")

    if (is_tergm) {
      mcmc_control_name <- paste(c("mcmc.control.tergm", if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula.tergm", "init_tergmLite", dat$control[[mcmc_control_name]])
      ## enforce some specific values appropriate for tergmLite/EpiModel netsim
      dat$control$mcmc.control[[i]]$MCMC.samplesize <- 1L
      dat$control$mcmc.control[[i]]$time.burnin <- 0L
      dat$control$mcmc.control[[i]]$time.interval <- 1L
      dat$control$mcmc.control[[i]]$time.samplesize <- 1L
      dat$control$mcmc.control[[i]]$changes <- TRUE
      dat$control$mcmc.control[[i]]$collect <- FALSE

      formation <- dat$nwparam[[i]]$formation
      dissolution <- dat$nwparam[[i]]$coef.diss$dissolution
      dat$nwparam[[i]]$tergm_formula <- trim_env(~Form(formation) + Diss(dissolution), keep = c("formation", "dissolution"))

      proposal <- ergm_proposal(nwp$constraints, hints = dat$control$mcmc.control[[i]]$MCMC.prop, arguments = dat$control$mcmc.control[[i]]$MCMC.prop.args, weights = dat$control$mcmc.control[[i]]$MCMC.prop.weights, nw = nw, class = "t")
      model <- ergm_model(dat$nwparam[[i]]$tergm_formula, nw = nw, term.options = dat$control$mcmc.control[[i]]$term.options, extra.aux=list(proposal=proposal$auxiliaries, system=trim_env(~.lasttoggle)), dynamic = TRUE)

      term_names <- unlist(c(lapply(model$terms[[1]]$submodel$terms, function(x) x$name), lapply(model$terms[[2]]$submodel$terms, function(x) x$name)))

      nwstats_formula <- NVL(dat$control[[nwstats_formula_name]], trim_env(~.))
      if (is.character(nwstats_formula)) {
        nwstats_formula <- switch(nwstats_formula,
                                  formation = formation,
                                  dissolution = dissolution,
                                  all = {formula_addition <- append_rhs.formula(~., dissolution, keep.onesided = TRUE);
                                         environment(formula_addition) <- environment(dissolution);
                                         nonsimp_update.formula(formation, formula_addition, from.new = TRUE)})
      }
      dat$control$nwstats.formulas[[i]] <- nwstats_formula

      if (dat$control$tergmLite.track.duration) {
        if (is.null(nw %n% "time")) nw %n% "time" <- 0
        if (is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- matrix(0L, nrow = 0, ncol = 3)
      }
    } else {
      mcmc_control_name <- paste(c("mcmc.control.ergm", if (num_nw > 1) i), collapse = ".")
      dat$control$mcmc.control[[i]] <- check.control.class("simulate.formula", "init_tergmLite", dat$control[[mcmc_control_name]])
      ## enforce some specific values appropriate for tergmLite/EpiModel netsim
      dat$control$mcmc.control[[i]]$MCMC.samplesize <- 1L

      proposal <- ergm_proposal(nwp$constraints, hints = dat$control$mcmc.control[[i]]$MCMC.prop, arguments = dat$control$mcmc.control[[i]]$MCMC.prop.args, weights = dat$control$mcmc.control[[i]]$MCMC.prop.weights, nw = nw, class = "c")
      model <- ergm_model(dat$nwparam[[i]]$formation, nw = nw, term.options = dat$control$mcmc.control[[i]]$term.options,  extra.aux=list(proposal=proposal$auxiliaries))

      term_names <- unlist(lapply(model$terms, function(x) x$name))

      nwstats_formula <- NVL(dat$control[[nwstats_formula_name]], trim_env(~.))
      if (is.character(nwstats_formula)) {
        nwstats_formula <- switch(nwstats_formula,
                                  formation = dat$nwparam[[i]]$formation,
                                  dissolution = trim_env(~.),
                                  all = dat$nwparam[[i]]$formation)
      }
      dat$control$nwstats.formulas[[i]] <- nwstats_formula
    }

    proposal$aux.slots <- model$slots.extra.aux$proposal
    dat$p[[i]]$state <- ergm_state(nw, model=model, proposal=proposal, stats=rep(0,nparam(model, canonical=TRUE)))

    model_mon <- ergm_model(dat$control$nwstats.formulas[[i]], nw = nw, term.options = dat$control$mcmc.control[[i]]$term.options, dynamic = TRUE)
    term_names <- c(term_names, unlist(lapply(model_mon$terms, function(x) x$name)))

    ## check for unsupported terms
    difference <- setdiff(term_names, supported.terms)
    if (length(difference) > 0) {
      ## special error message for gwesp without fixed=TRUE
      if (difference[1] == "esp") {
        stop("tergmLite does not support the esp term, and only supports the gwesp term with fixed=TRUE")
      }
      ## error message for all others
      stop("tergmLite does not know how to update the term ", difference[1])
    }
  }

  dat$nw <- NULL

  return(dat)
}
