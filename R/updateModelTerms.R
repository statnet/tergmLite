
#' @title Methods for Computing and Updating ERGM/STERGM Term Inputs
#'
#' @description Function to appropriately update model inputs based on ERGM model
#'              terms when using networkLite representation.
#'
#' @param dat EpiModel dat object tracking simulation state
#' @param network Numberic number of network location for multi-network simulations.
#'
#' @details Calls \code{ergm_model} to update model inputs based on potential
#' exogenous changes to network structure (e.g., number of nodes) or nodal attributes
#' used within ERGM model (see example below). This function is typically used
#' within \code{EpiModel} module for network resimulation, immediately prior to
#' calling \code{\link{simulate_network}} or \code{\link{simulate_ergm}}.
#'
#' Implemented terms are:
#'  \itemize{
#'    \item \code{edges}
#'    \item \code{nodematch}
#'    \item \code{nodefactor}
#'    \item \code{concurrent} (including heterogenous by attribute)
#'    \item \code{degree} (including heterogenous by attribute)
#'    \item \code{degrange}
#'    \item \code{absdiff}
#'    \item \code{absdiffby} (in the EpiModel package)
#'    \item \code{nodecov}
#'    \item \code{nodemix}
#'    \item \code{absdiffnodemix} (in the EpiModel package)
#'    \item \code{triangle}
#'    \item \code{gwesp(fixed=TRUE)}
#'  }
#' All other ERGM terms will return errors.
#'
#' @return
#' Returns an updated \code{dat} object with the network list structure inputs
#' used by \code{\link{simulate_network}} or \code{\link{simulate_ergm}} with changes
#' to network size or nodal covariates.
#'
#' @export
#' @importFrom statnet.common NVL
#'
#' @examples
#' library("EpiModel")
#'
#' # Set seed for reproducibility
#' set.seed(12345)
#'
#' nw <- network.initialize(n = 100, directed = FALSE)
#' nw <- set.vertex.attribute(nw, "group", rep(0:1, each = 50))
#' formation <- ~edges + nodefactor("group")
#' target.stats <- c(15, 10)
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10, i.num.g2 = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, depend = TRUE)
#'
#' # Full network structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#' str(dat, max.level = 1)
#'
#' # networkLite representation used by tergmLite
#' dat <- init_tergmLite(dat)
#'
#' # Examine the network list structure for nodefactor term
#' dat$p[[1]]$model.form$terms[[2]]
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
#' dat$p[[1]]$model.form$terms[[2]]
#'
updateModelTermInputs <- function(dat, network = 1) {

  p <- dat$p[[network]]

  if ("model.diss" %in% names(p)) {
    dynamic <- TRUE
  } else {
    dynamic <- FALSE
  }

  if (dynamic == TRUE) {
    mf <- p$model.form
    md <- p$model.diss
    mhf <- p$MHproposal.form
    mhd <- p$MHproposal.diss
  } else {
    mf <- p$model.form
    mhf <- p$MHproposal
  }

  n <- attributes(dat$el[[network]])$n

  ## 1. Formation Model

  updates <- ergm_model(dat$nwparam[[network]]$formation, networkLite(dat$el[[network]], dat$attr))

  for (t in seq_along(mf$terms)) {

    # pull term name
    term <- mf$terms[[t]]

    supported.terms <- c("edges", "nodematch", "nodefactor",
                         "concurrent", "concurrent_by_attr",
                         "degree", "degree_by_attr",
                         "absdiff", "absdiffby", "nodecov", "nodemix",
                         "absdiffnodemix", "degrange", "triangle", "gwesp")
    if (!(term$name %in% supported.terms)) {
      ## special error message for gwesp without fixed=TRUE
      if (term$name == "esp") {
        stop("tergmLite does not support the esp term, and only supports the gwesp term with fixed=TRUE")
      }
      ## error message for all others
      stop("tergmLite does not know how to update the term ",
           term$name," in the formation model formula")
    }

    mf$terms[[t]] <- updates$terms[[t]]
  }

  mf$maxval <- updates$maxval


  ## 2. Dissolution Model

  # loop over dissolution model terms and update
  if (dynamic == TRUE) {
    for (t in seq_along(md$terms)) {

      term <- md$terms[[t]]

      supported.terms <- c("edges", "nodematch", "nodefactor", "nodemix")
      if (!(term$name %in% supported.terms)) {
        stop("tergmLite does not know how to update the term ",
             term$name," in the formation model formula")
      }

      md$terms[[t]] <- updates$terms[[t]]

      # in case offset is used in dissolution model but not formation model
      md$terms[[t]]$coef.names <- term$coef.names
    }
    for (j in seq_along(md$terms)) {
      if (j == 1) {
        new.maxval <- rep(if (!is.null(md$terms[[j]]$maxval)) md$terms[[j]]$maxval else +Inf,
                          length.out = length(md$terms[[j]]$coef.names))
      } else {
        new.maxval <- c(new.maxval,
                        rep(if (!is.null(md$terms[[j]]$maxval)) md$terms[[j]]$maxval else +Inf,
                            length.out = length(md$terms[[j]]$coef.names)))
      }
    }
    md$maxval <- new.maxval

  }


  ## 3. MHproposal Lists (constraints)

  ## Update MHproposal.form
  hasBD <- !is.null(mhf$arguments$constraints$bd$attribs[1])
  if (hasBD == TRUE) {
    bd <- mhf$arguments$constraints$bd
    bd$attribs <- matrix(rep(bd$attribs[1], n), ncol = 1)
    bd$maxout <- matrix(rep(bd$maxout[1], n), ncol = 1)
    bd$maxin <- matrix(rep(n - 1, n), ncol = 1)
    bd$minout <- matrix(rep(0, n), ncol = 1)
    bd$minin <- matrix(rep(0, n), ncol = 1)
    mhf$arguments$constraints$bd <- bd
  }

  # MHproposal.diss (currently matches mhf bd constraint)
  if (dynamic == TRUE) {
    mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd
  }


  ## 4. Update and return

  # update the elements of the parameter list and return
  if (dynamic == TRUE) {
    p <- list(model.form = mf, model.diss = md,
              MHproposal.form = mhf, MHproposal.diss = mhd)
  } else {
    p <- list(model.form = mf, MHproposal = mhf)
  }

  dat$p[[network]] <- p

  return(dat)
}
