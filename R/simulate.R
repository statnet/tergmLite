
# STERGM --------------------------------------------------------------

#' @title A Version of tergm::simulate.network Tailored to tergmLite Simulation
#'
#' @description Resimulates an edgelist given ergm_state data representation
#'              and model coefficients.
#'
#' @param state An \code{ergm_state} object representing the starting state
#'              for the simulation.
#' @param coef Vector of coefficients for the generative model.
#' @param control A control list of class \code{control.simulate.network.tergm},
#'                with augmentations made in \code{init_tergmLite}.
#' @param save.changes Logical, if \code{TRUE}, saves a matrix of changed edges
#'                     as an attribute of the output edgelist matrix
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
#' Returns a named list with members \code{el} (the updated network edgelist
#' representation) and \code{state} (the \code{ergm_state} object returned by
#' \code{tergm_MCMC_slave}).  If \code{save.changes} is \code{TRUE}, also 
#' returns a list of new edges and dissolved edges with the resimulation,
#' attached to \code{el} as the \code{changes} attribute.
#'
#' @examples
#' \dontrun{
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
#' dat$el[[1]] <- simulate_network(state = dat$p[[1]]$state, 
#'                                 coef = c(dat$nwparam[[1]]$coef.form, 
#'                                          dat$nwparam[[1]]$coef.diss$coef.adj),
#'                                 control = dat$control$mcmc.control[[1]],
#'                                 save.changes = TRUE)$el
#' dat$el[[1]]
#'
#' # Specific changes listed under changes list
#' #    (new edges: to = 1; dissolved edges: to = 0):
#' attributes(dat$el[[1]])$changes
#' }
#'
simulate_network <- function(state, 
                             coef,
                             control,
                             save.changes = FALSE) {
  z <- tergm_MCMC_slave(state, ergm.eta(coef, state$model$etamap), control, verbose = FALSE)

  ## need matrix not tibble
  el <- as.matrix(z$state$el)
  storage.mode(el) <- "integer"
  
  if (save.changes == TRUE) {
    if (z$diffnwtime[1] > 0) {
      changes <- cbind(z$diffnwtails[2:(z$diffnwtails[1] + 1)],
                       z$diffnwheads[2:(z$diffnwheads[1] + 1)],
                       z$diffnwdirs[2:(z$diffnwdirs[1] + 1)])
    } else {
      changes <- matrix(0, ncol = 3, nrow = 0)
    }
    colnames(changes) <- c("tail", "head", "to")
    attributes(el)$changes <- changes
  }

  attributes(el)$n <- attributes(state$el)$n
  attributes(el)$directed <- attributes(state$el)$directed
  attributes(el)$bipartite <- attributes(state$el)$bipartite
  attributes(el)$loops <- attributes(state$el)$loops
  attributes(el)$inverted <- attributes(state$el)$inverted
  class(el) <- c("edgelist", class(el))
  
  return(list(el = el, state = z$state))
}


# ERGM ----------------------------------------------------------------

#' @title A Version of ergm::simulate.ergm Tailored to tergmLite Simulation
#'
#' @description Resimulates an edgelist given ergm_state data representation
#'              and model coefficients.
#'
#' @param state An \code{ergm_state} object representing the starting state
#'              for the simulation.
#' @param coef Vector of coefficients for the generative model.
#' @param control A control list of class \code{control.simulate.formula}.
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
#' Returns a named list with members \code{el} (the updated network edgelist
#' representation) and \code{state} (the \code{ergm_state} object returned by
#' \code{ergm_MCMC_slave}).
#'
#' @examples
#' \dontrun{
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
#' dat$el[[1]] <- simulate_ergm(state = dat$p[[1]]$state,
#'                              coef = dat$nwparam[[1]]$coef.form,
#'                              control = dat$control$mcmc.control[[1]])$el
#' dat$el[[1]]
#' }
#'
simulate_ergm <- function(state, coef, control) {
  z <- ergm_MCMC_slave(state, ergm.eta(coef, state$model$etamap), control, verbose = FALSE)
  el <- structure(as.matrix(z$state$el), n = attr(state$el, "n"), class = c("edgelist", "matrix"))
  storage.mode(el) <- "integer"
  list(el = el, state = z$state)
}
