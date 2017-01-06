
#' @title Alternate Methods for Computing Ergm Term Inputs
#'
#' @description Function to appropriately update model params based on ergm model
#'              terms when using 'fast_edgelist' representations.
#'
#' @param dat EpiModel dat object tracking simulation state
#' @param network Numberic number of network location for multi-network simulations.
#'
#' @details Contains hard-coded implementations of some of the most commonly used
#' ergm terms called instead of InitErgmTerm.x, checkErgmTerm, etc.
#'
#' Implemented terms are:
#'  \itemize{
#'    \item edges
#'    \item nodematch
#'    \item nodefactor
#'    \item concurrent
#'    \item degree
#'    \item absdiff
#'    \item nodecov
#'    \item nodemix
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
  maxdyads <- choose(n, 2)

  # All input vectors are padded out as the would be by updatemodel.ErgmTerm
  # outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp),
  #                     length(outlist$coef.names),
  #                     length(outlist$inputs), outlist$inputs)


  ## 1. Formation Model ##

  for (t in seq_along(mf$terms)) {

    # pull term name
    term <- mf$terms[[t]]

    if (term$name == "edges") {

      ## Reference: ergm:::InitErgmTerm.edges

      # Only need to update maxval, no update to inputs
      mf$terms[[t]]$maxval <- maxdyads
    }

    else if (term$name == "nodematch") {

      ## Reference: ergm:::InitErgmTerm.nodematch

      # Get the formation formula to parse the params
      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form, t)
      # get the name of the attribute to be used for nodecov
      attrname <- args[[1]]
      diff <- args$diff
      # collect the values for the attribute
      nodecov <- dat$attr[[attrname]]
      u <- sort(unique(nodecov))
      # optionally remove values not indicated by "keep:
      if (!is.null(args$keep)) {
        u <- u[args$keep]
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      dontmatch <- nodecov == (length(u) + 1)
      nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
      ui <- seq(along = u)
      if (diff == TRUE) {
        inputs <- c(ui, nodecov)
      } else {
        inputs <- nodecov
      }
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)

    }

    else if (term$name == "nodefactor") {

      ## Reference: ergm:::InitErgmTerm.nodefactor

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form, t)

      attrname <- args[[1]]

      # Handles interaction terms: nodefactor("attr1", "attr2")
      if (length(attrname) == 1) {
        nodecov <- dat$attr[[attrname]]
      } else {
        nodecov <- do.call(paste, c(sapply(attrname,
                                           function(oneattr) dat$attr[[oneattr]],
                                           simplify = FALSE), sep = "."))
      }

      u <- sort(unique(nodecov))
      if (any(statnet.common::NVL(args$base, 0) != 0)) {
        u <- u[-args$base]
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)

      mf$terms[[t]]$inputs <- c(length(ui), length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    }

    else if (term$name == "concurrent") {

      # Reference: ergm:::InitErgmTerm.concurrent

      # Homogenous form only updates maxval
      mf$terms[[t]]$maxval <- n

    }

    else if (term$name == "concurrent_by_attr") {

      # Reference: ergm:::InitErgmTerm.concurrent

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form, t)

      # heterogeneous form updates input vector
      byarg <- args$by
      if (!is.null(byarg)) {
        nodecov <- dat$attr[[byarg]]
        u <- sort(unique(nodecov))
        if (any(is.na(nodecov))) {
          u <- c(u, NA)
        }
        nodecov <- match(nodecov, u)
        lu <- length(u)
        ui <- seq(along = u)
        inputs <- c(ui, nodecov)
        mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                  length(inputs), inputs)
      }

      mf$terms[[t]]$maxval <- n

    }

    else if (term$name == "degree") {

      ## TODO: check this
      ## Reference ergm:::InitErgmTerm.degree

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form,t)

      d <- args[[1]]
      byarg <- args$byarg
      homophily <- args$homophily
      emptynwstats <- NULL
      if (!is.null(byarg)) {
        nodecov <- dat$attr[[byarg]]
        u <- sort(unique(nodecov))
        if (any(is.na(nodecov))) {
          u <- c(u, NA)
        }
        nodecov <- match(nodecov, u)
        if (length(u) == 1) {
          stop("Attribute given to degree() has only one value",
               call. = FALSE)
        }
      }
      if (!is.null(byarg) && !homophily) {
        lu <- length(u)
        du <- rbind(rep(d, lu), rep(1:lu, rep(length(d), lu)))
        if (any(du[1, ] == 0)) {
          emptynwstats <- rep(0, ncol(du))
          tmp <- du[2, du[1, ] == 0]
          for (i in 1:length(tmp)) tmp[i] <- sum(nodecov ==
                                                   tmp[i])
          emptynwstats[du[1, ] == 0] <- tmp
        }
      }  else {
        if (any(d == 0)) {
          emptynwstats <- rep(0, length(d))
          emptynwstats[d == 0] <- attr(dat$el,'n') # network size
        }
      }
      if (is.null(byarg)) {
        if (length(d) == 0) {
          return(NULL)
        }
        inputs <- c(d)
      } else if (homophily) {
        if (length(d) == 0) {
          return(NULL)
        }
        inputs <- c(d, nodecov)
      } else {
        if (ncol(du) == 0) {
          return(NULL)
        }
        inputs <- c(as.vector(du), nodecov)
      }
      if (!is.null(emptynwstats)) {
        mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                  length(inputs), inputs)
        mf$terms[[t]]$emptynwstats <- emptynwstats

      } else {
        mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                  length(inputs), inputs)
        mf$terms[[t]]$maxval <- attr(dat$el, "n") # network size
      }


    }

    else if (term$name == "absdiff") {

      # Reference: ergm:::InitErgmTerm.absdiff

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form, t)
      attrname <- args[[1]]
      # Transformation function
      pow <- args$pow
      inputs <- c(pow, dat$attr[[attrname]])
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    }

    else if (term$name == "nodecov") {

      ## Reference: ergm:::InitErgmTerm.nodecov

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form,t)
      attrname <- args[[1]]
      # get the transformation function
      f <- args$transform
      nodecov <- dat$attr[[attrname]]
      # strange that the original version of the term doesn't require this
      #    logic to implement the default term..
      if (!is.null(f)) {
        inputs <- f(nodecov)
      } else {
        inputs <- nodecov
      }
      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)
    }

    else if (term$name == "nodemix") {

      # ergm:::InitErgmTerm.nodemix
      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form,t)
      attrname <- args[[1]]
      nodecov <- dat$attr[[attrname]]
      base <- args$base
      # ASSUMES NETWORK IS NOT BIPARTITE
      u <- sort(unique(nodecov))
      if (any(is.na(nodecov))) {
        u <- c(u, NA)
      }
      nodecov <- match(nodecov, u, nomatch = length(u) + 1)
      ui <- seq(along = u)
      uui <- matrix(1:length(ui)^2, length(ui), length(ui))
      urm <- t(sapply(ui, rep, length(ui)))
      ucm <- sapply(ui, rep, length(ui))
      uun <- outer(u, u, paste, sep = ".")
      uui <- uui[upper.tri(uui, diag = TRUE)]
      urm <- urm[upper.tri(urm, diag = TRUE)]
      ucm <- ucm[upper.tri(ucm, diag = TRUE)]
      uun <- uun[upper.tri(uun, diag = TRUE)]
      if (any(NVL(base, 0) != 0)) {
        urm <- as.vector(urm)[-base]
        ucm <- as.vector(ucm)[-base]
        uun <- as.vector(uun)[-base]
      }
      inputs <- c(urm, ucm, nodecov)
      attr(inputs, "ParamsBeforeCov") <- 2 * length(uun)
      mf$terms[[t]]$inputs <- c(2 * length(uun), length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)

    }

    else if (term$name == "absdiffnodemix") {

      # Reference: EpiModelHIV::InitErgmTerm.absdiffnodemix

      form <- dat$nwparam[[network]]$formation
      args <- get_formula_term_args_in_formula_env(form, t)

      nodecov <- dat$attr[[args[[1]]]]
      nodecovby <- dat$attr[[args[[2]]]]

      u <- sort(unique(nodecovby))
      if (any(is.na(nodecovby))) {
        u <- c(u, NA)
      }
      nodecovby <- match(nodecovby, u, nomatch = length(u) + 1)
      ui <- seq(along = u)

      uui <- matrix(1:length(ui) ^ 2, length(ui), length(ui))
      urm <- t(sapply(ui, rep, length(ui)))
      ucm <- sapply(ui, rep, length(ui))
      uun <- outer(u, u, paste, sep = ".")
      uui <- uui[upper.tri(uui, diag = TRUE)]
      urm <- urm[upper.tri(urm, diag = TRUE)]
      ucm <- ucm[upper.tri(ucm, diag = TRUE)]
      uun <- uun[upper.tri(uun, diag = TRUE)]

      inputs <- c(length(nodecov), length(urm), nodecov, nodecovby, urm, ucm)

      mf$terms[[t]]$inputs <- c(0, length(mf$terms[[t]]$coef.names),
                                length(inputs), inputs)

    }

    else {
      stop("tergmLite does not know how to update the term ",
           term$name," in the formation model formula")
    }

  }

  # update combinded maxval
  for (j in seq_along(mf$terms)) {
    if (j == 1) {
      new.maxval <- rep(if (!is.null(mf$terms[[j]]$maxval)) mf$terms[[j]]$maxval else +Inf,
                        length.out = length(mf$terms[[j]]$coef.names))
    } else {
      new.maxval <- c(new.maxval,
                      rep(if (!is.null(mf$terms[[j]]$maxval)) mf$terms[[j]]$maxval else +Inf,
                          length.out = length(mf$terms[[j]]$coef.names)))
    }
  }
  mf$maxval <- new.maxval


  ## 2. Dissolution Model ##

  # loop over dissolution model terms and update
  if (dynamic == TRUE) {
    for (t in seq_along(md$terms)) {
      term <- md$terms[[t]]

      if (term$name == "edges") {
        md$terms[[t]]$maxval <- maxdyads
      }

      else if (term$name == "nodematch") {

        ## Reference: ergm:::InitErgmTerm.nodematch

        diss <- dat$nwparam[[network]]$coef.diss$dissolution
        args <- get_formula_term_args_in_formula_env(diss, t)
        diff <- args$diff
        nodecov <- dat$attr[[attrname]]
        u <- sort(unique(nodecov))
        if (!is.null(args$keep)) {
          u <- u[args$keep]
        }
        nodecov <- match(nodecov, u, nomatch = length(u) + 1)
        dontmatch <- nodecov == (length(u) + 1)
        nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
        ui <- seq(along = u)
        if (diff == TRUE) {
          inputs <- c(ui, nodecov)
        } else {
          inputs <- nodecov
        }
        md$terms[[t]]$inputs <- c(0, length(md$terms[[t]]$coef.names),
                                  length(inputs), inputs)

      }

      else if (term$name == "nodefactor") {

        ## Reference: ergm:::InitErgmTerm.nodefactor

        diss <- dat$nwparam[[network]]$coef.diss$dissolution
        args <- get_formula_term_args_in_formula_env(diss, t)
        attrname <- args[[1]]
        # Handles interaction terms: nodefactor("attr1", "attr2")
        if (length(attrname) == 1) {
          nodecov <- dat$attr[[attrname]]
        } else {
          nodecov <- do.call(paste, c(sapply(attrname,
                                             function(oneattr) dat$attr[[oneattr]],
                                             simplify = FALSE), sep = "."))
        }

        u <- sort(unique(nodecov))
        if (any(statnet.common::NVL(args$base, 0) != 0)) {
          u <- u[-args$base]
        }
        nodecov <- match(nodecov, u, nomatch = length(u) + 1)
        ui <- seq(along = u)
        inputs <- c(ui, nodecov)
        md$terms[[t]]$inputs <- c(length(ui), length(md$terms[[t]]$coef.names),
                                  length(inputs), inputs)
      }

      else if (term$name == "nodemix") {

        # ergm:::InitErgmTerm.nodemix
        diss <- dat$nwparam[[network]]$coef.diss$dissolution
        args <- get_formula_term_args_in_formula_env(diss, t)
        attrname <- args[[1]]
        nodecov <- dat$attr[[attrname]]
        base <- args$base
        u <- sort(unique(nodecov))
        if (any(is.na(nodecov))) {
          u <- c(u, NA)
        }
        nodecov <- match(nodecov, u, nomatch = length(u) + 1)
        ui <- seq(along = u)
        uui <- matrix(1:length(ui)^2, length(ui), length(ui))
        urm <- t(sapply(ui, rep, length(ui)))
        ucm <- sapply(ui, rep, length(ui))
        uun <- outer(u, u, paste, sep = ".")
        uui <- uui[upper.tri(uui, diag = TRUE)]
        urm <- urm[upper.tri(urm, diag = TRUE)]
        ucm <- ucm[upper.tri(ucm, diag = TRUE)]
        uun <- uun[upper.tri(uun, diag = TRUE)]
        if (any(NVL(base, 0) != 0)) {
          urm <- as.vector(urm)[-base]
          ucm <- as.vector(ucm)[-base]
          uun <- as.vector(uun)[-base]
        }
        inputs <- c(urm, ucm, nodecov)
        md$terms[[t]]$inputs <- c(2 * length(uun), length(md$terms[[t]]$coef.names),
                                  length(inputs), inputs)

      }

      else {
        stop("tergmLite does not know how to update the term ",
             term$name, " in the dissolution model formula")
      }

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


  ## 3. MHproposal Lists (for constraints) ##

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
    mhd$arguments$constraints$bd <- bd
  }

  # MHproposal.diss (currently matches mhf bd constraint)
  if (dynamic == TRUE) {
    mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd
  }


  ## 4. Update and return ##

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



#' @title Evaluate Ergm Model Formula Terms
#' @description This is a work-around for evaluating model terms in the non-standard
#'              tergmLite sequence, as an alternative to ergm.getModel. Computes
#'              a list of the arguments to the terms in the formula with offsets
#'              removed, evaluated in the formula calling environment returns a
#'              list where the first element is the term name and subsequent (named)
#'              elements are the argument values named by the argument names.
#'
#' @param form an ergm model formula
#' @param termIndex an integer index for the formula term in form to be evaluated
#'
#' @export
#' @importFrom statnet.common term.list.formula
#'
get_formula_term_args_in_formula_env <- function(form, termIndex) {

  # get the calling environment of the formula in case
  # there are substitutions
  formula.env <- environment(form)
  args <- statnet.common::term.list.formula(form[[2]])[[termIndex]]

  # remove the offset term if it exists
  if (args[1] == "offset()") {
    args <- args[[-1]]
  }

  # term name
  tname <- args[[1]]

  # hack to convert from a call to a list when evaluated
  args[[1]] <- as.name("list")

  # evaluate in formula's calling environment
  outlist <- eval(args, formula.env)

  # Set default base to 1
  if (tname == "nodefactor" & is.null(outlist$base)) {
    outlist$base <- 1
  }

  # Set default pow to 1
  if (tname == "absdiff" & is.null(outlist$pow)) {
    outlist$pow <- 1
  }

  # set default fixed argument to FALSE
  # also needs default keep argument, set to NULL
  if (tname == "nodematch" & is.null(args$diff)) {
    outlist$diff <- FALSE
  }

  return(outlist)
}
