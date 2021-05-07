
#' @title Methods for Computing and Updating ERGM/TERGM Simulation Inputs
#'
#' @description Function to appropriately update the ergm_state object used
#'              to start ERGM/TERGM simulation within the context of a tergmLite
#'              simulation.
#'
#' @param dat EpiModel dat object tracking simulation state
#' @param network Numeric number of network location for multi-network simulations.
#'
#' @details Calls \code{ergm_model} to update model inputs based on potential
#' exogenous changes to network structure (e.g., number of nodes) or nodal attributes
#' used within ERGM model (see example below). This function is typically used
#' within \code{EpiModel} module for network resimulation, immediately prior to
#' calling \code{\link{simulate_network}} or \code{\link{simulate_ergm}}. Both
#' \code{ergm_proposal} and \code{ergm_state} are also called to fully prepare
#' for the next round of simulation.
#'
#' @return
#' Returns an updated \code{dat} object with the network list structure inputs
#' used by \code{\link{simulate_network}} or \code{\link{simulate_ergm}} with changes
#' to network size or nodal covariates.
#'
#' @export
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
#' # Examine the network list structure for nodefactor term
#' dat$p[[1]]$state$model$terms[[2]]
#'
#' # inputs vector corresponds to group attribute stored here
#' dat$attr$group
#'
#' # As example of what could happen in EpiModel: randomly reshuffle group
#' #   attribute values of 100 nodes
#' dat$attr$group <- sample(dat$attr$group)
#' dat$attr$group
#'
#' # Update network list structure
#' dat <- updateModelTermInputs(dat)
#'
#' # Check that network list structure for nodefactor term has been updated
#' dat$p[[1]]$state$model$terms[[2]]
#' }
#'
updateModelTermInputs <- function(dat, network = 1) {
  nwL <- networkLite(dat$el[[network]], dat$attr)
  
  if (is(dat$control$MCMC_control[[network]], "control.simulate.network.tergm")) { # dynamic
    proposal <- ergm_proposal(dat$nwparam[[network]]$constraints, hints = dat$control$MCMC_control[[network]]$MCMC.prop, arguments = dat$control$MCMC_control[[network]]$MCMC.prop.args, weights = dat$control$MCMC_control[[network]]$MCMC.prop.weights, nw = nwL, class = "t")
    model <- ergm_model(dat$nwparam[[network]]$tergm_formula, nw = nwL, term.options = dat$control$MCMC_control[[network]]$term.options, extra.aux=list(proposal=proposal$auxiliaries, system=~.lasttoggle))
    
    if(dat$control$track_duration) {
      nwL %n% "time" <- dat$p[[network]]$state$nw0 %n% "time"
      nwL %n% "lasttoggle" <- dat$p[[network]]$state$nw0 %n% "lasttoggle"
    }
  } else { # static
    proposal <- ergm_proposal(dat$nwparam[[network]]$constraints, hints = dat$control$MCMC_control[[network]]$MCMC.prop, arguments = dat$control$MCMC_control[[network]]$MCMC.prop.args, weights = dat$control$MCMC_control[[network]]$MCMC.prop.weights, nw = nwL, class = "c")
    model <- ergm_model(dat$nwparam[[network]]$formation, nw = nwL, term.options = dat$control$MCMC_control[[network]]$term.options,  extra.aux=list(proposal=proposal$auxiliaries))
  }
  
  proposal$aux.slots <- model$slots.extra.aux$proposal
  dat$p[[network]]$state <- ergm_state(nwL, model=model, proposal=proposal, stats=rep(0,nparam(model, canonical=TRUE)))
  
  if(nparam(dat$p[[network]]$state_mon, canonical = TRUE) > 0) {
    model_mon <- ergm_model(dat$control$monitors[[network]], nw = nwL, term.options = dat$control$MCMC_control[[network]]$term.options)
    dat$p[[network]]$state_mon <- ergm_state(nwL, model=model_mon, proposal=NULL, stats=rep(0, nparam(model_mon, canonical=TRUE)))
  }
  
  return(dat)
}
