
#' @title Prepare Network and STERGM Objects for tergmLite
#'
#' @description Converts network object, formation and dissolution formulas,
#'              formation and dissolution coefficients, and control settings to a
#'              thin list format for STERGM resimulation.
#'
#' @param nw An object of class \code{network}
#' @param formation Right-hand sided formation formula
#' @param dissolution Right-hand sided dissolution formula
#' @param coef.form Vector of coefficients associated with the formation formula
#' @param coef.diss Vector of coefficients associated with the dissolution formula
#' @param constraints Constraints for the formation model (only \code{bd})
#'                    constraints currently supported.
#' @param control Control settings passed to \code{tergm::control.simulate.network}
#'
#' @details
#' This is an internal function used within \code{\link{init_tergmLite}}.
#'
stergm_prep <- function(nw,
                        formation,
                        dissolution,
                        coef.form,
                        coef.diss,
                        constraints,
                        control = control.simulate.network()) {

  if (inherits(nw, "network") == FALSE) {
    stop("A network object must be given")
  }

  formation <- statnet.common::nonsimp_update.formula(formation, nw ~ ., from.new = "nw")
  dissolution <- statnet.common::nonsimp_update.formula(dissolution, nw ~ ., from.new = "nw")

  model.form <- ergm::ergm_model(formation, nw, role = "formation")
  if (!missing(coef.form) && ergm::nparam(model.form) != length(coef.form)) {
    stop("coef.form has ", length(coef.form), " elements, while the model requires ",
         ergm::nparam(model.form), " parameters.")
  }
  model.diss <- ergm::ergm_model(dissolution, nw, role = "dissolution")
  if (!missing(coef.diss) && ergm::nparam(model.diss) != length(coef.diss)) {
    stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",
         ergm::nparam(model.diss), " parameters.")
  }

  MHproposal.form <- ergm::ergm_proposal(constraints, control$MCMC.prop.args.form,
                                      nw, weights = control$MCMC.prop.weights.form,
                                      class = "f")
  MHproposal.diss <- ergm::ergm_proposal(constraints, control$MCMC.prop.args.diss,
                                      nw, weights = control$MCMC.prop.weights.diss,
                                      class = "d")

  MHproposal.form$arguments$constraints$.attributes <- NULL
  MHproposal.diss$arguments$constraints$.attributes <- NULL

  out <- list(model.form = model.form, model.diss = model.diss,
              MHproposal.form = MHproposal.form, MHproposal.diss = MHproposal.diss)

  return(out)
}

#' @title Prepare Network and ERGM Objects for tergmLite
#'
#' @description Converts network object, formation and dissolution formulas,
#'              formation and dissolution coefficients, and control settings to a
#'              thin list format for ERGM resimulation.
#'
#' @param nw An object of class \code{network}
#' @param formation Right-hand sided formation formula
#' @param coef Vector of coefficients associated with the formation formula
#' @param constraints Constraints for the formation model (only \code{bd})
#'                    constraints currently supported.
#' @param control Control settings passed to \code{ergm::control.simulate.ergm}
#'
#' @details
#' This is an internal function used within \code{\link{init_tergmLite}}.
#'
ergm_prep <- function(nw,
                      formation,
                      coef,
                      constraints,
                      control = ergm::control.simulate.ergm()) {

  form <- statnet.common::nonsimp_update.formula(formation, nw ~ ., from.new = "nw")
  m <- ergm::ergm_model(form, nw, response = NULL, role = "static")

  MHproposal <- ergm_proposal(constraints, arguments = control$MCMC.prop.args,
                           nw = nw, weights = control$MCMC.prop.weights, class = "c",
                           reference = ~Bernoulli, response = NULL)

  MHproposal$arguments$constraints$.attributes <- NULL

  out <- list(model.form = m, MHproposal = MHproposal)
  return(out)
}

#' @title Initializes EpiModel netsim Object for tergmLite Simulation
#'
#' @param dat A list object containing a \code{networkDynamic} object and other
#'        initialization information passed from \code{netsim}.
#'
#' @export
#'
#' @examples
#' library("EpiModel")
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, depend = TRUE)
#'
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' dat <- init_tergmLite(dat)
#' str(dat, max.level = 1)
#'
init_tergmLite <- function(dat) {

  num_nw <- ifelse(inherits(dat$nw, "network"), 1, length(dat$nw))

  dat$el <- list()
  dat$p <- list()

  for (i in 1:num_nw) {

    nwp <- dat$nwparam[[i]]
    is_tergm <- all(nwp$coef.diss$duration > 1)
    if (num_nw == 1) {
      nw <- dat$nw
    } else {
      nw <- dat$nw[[i]]
    }

    dat$el[[i]] <- as.edgelist(nw)
    attributes(dat$el[[i]])$vnames <- NULL

    if (is_tergm) {
      p <- stergm_prep(nw, nwp$formation, nwp$coef.diss$dissolution,
                       nwp$coef.form, nwp$coef.diss$coef.adj, nwp$constraints)
      p$model.form$formula <- NULL
      p$model.diss$formula <- NULL
    } else {
      p <- tergmLite::ergm_prep(nw, nwp$formation, nwp$coef.form, nwp$constraints)
      p$model.form$formula <- NULL
    }
    dat$p[[i]] <- p

  }

  dat$nw <- NULL

  return(dat)
}
