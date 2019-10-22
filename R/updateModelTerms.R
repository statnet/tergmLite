
#' @title Methods for Computing and Updating ERGM/STERGM Term Inputs
#'
#' @description Function to appropriately update model params based on ergm model
#'              terms when using Edgelist-formatted representations
#'
#' @param dat EpiModel dat object tracking simulation state
#' @param network Numberic number of network location for multi-network simulations.
#'
#' @details Calls ergm_model.
#'
#' Implemented terms are:
#'  \itemize{
#'    \item edges
#'    \item nodematch
#'    \item nodefactor
#'    \item concurrent (including heterogenous by attribute)
#'    \item degree (including heterogenous by attribute)
#'    \item degrange
#'    \item absdiff
#'    \item absdiffby (in the EpiModelHIV package)
#'    \item nodecov
#'    \item nodemix
#'    \item absdiffnodemix (in the EpiModelHIV package)
#'    \item triangle
#'    \item gwesp(fixed=TRUE)
#'  }
#'  All other terms will return errors.
#'
#' @export
#' @importFrom statnet.common NVL
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
      if(term$name == "esp")
        stop("tergmLite only supports the gwesp term with fixed=TRUE")
      
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
