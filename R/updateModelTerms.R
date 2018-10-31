
#' @title Alternate Methods for Computing and Updating ERGM/STERGM Term Inputs
#'
#' @description Function to appropriately update model params based on ergm model
#'              terms when using Edgelist-formatted representations
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
#'    \item concurrent (including heterogenous by attribute)
#'    \item degree (including heterogenous by attribute)
#'    \item degrange
#'    \item absdiff
#'    \item absdiffby (in the EpiModelHIV package)
#'    \item nodecov
#'    \item nodemix
#'    \item absdiffnodemix (in the EpiModelHIV package)
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

  for (t in seq_along(mf$terms)) {

    # pull term name
    term <- mf$terms[[t]]

    supported.terms <- c("edges", "nodematch", "nodefactor",
                         "concurrent", "concurrent_by_attr",
                         "degree", "degree_by_attr",
                         "absdiff", "absdiffby", "nodecov", "nodemix",
                         "absdiffnodemix", "degrange")
    if (!(term$name %in% supported.terms)) {
      stop("tergmLite does not know how to update the term ",
           term$name," in the formation model formula")
    }

    func <- paste("updateInputs", term$name, sep = "_")
    mf$terms[[t]] <- do.call(func, list(dat = dat, network = network, term = term, t = t))
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

      func <- paste("updateInputs", term$name, sep = "_")
      md$terms[[t]] <- do.call(func, list(dat = dat, network = network, term = term, t = t))
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


updateInputs_edges <- function(dat, network, term, t) {

  # Only need to update maxval, no update to inputs
  n <- attributes(dat$el[[network]])$n
  maxdyads <- choose(n, 2)

  term$maxval <- maxdyads

  return(term)
}

updateInputs_nodematch <- function(dat, network, term, t) {

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
  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_nodefactor <- function(dat, network, term, t) {

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
  # nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  nodepos <- match(nodecov, u, nomatch = 0) - 1

  # inputs <- c(ui, nodecov)
  inputs <- nodepos

  # attr(inputs, "ParamsBeforeCov") <- length(ui)

  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_concurrent <- function(dat, network, term, t) {

  # Homogenous form only updates maxval
  n <- attributes(dat$el[[network]])$n
  term$maxval <- n

  return(term)
}

updateInputs_concurrent_by_attr <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)

  n <- attributes(dat$el[[network]])$n

  # heterogeneous form updates input vector
  byarg <- args$by
  if (!is.null(byarg)) {
    nodecov <- dat$attr[[byarg]]
    u <- sort(unique(nodecov))
    if (any(is.na(nodecov))) {
      u <- c(u, NA)
    }
    nodecov <- match(nodecov, u)
    ui <- seq(along = u)
    inputs <- c(ui, nodecov)
    term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
  }

  term$maxval <- n

  return(term)
}

updateInputs_degree <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)

  d <- args[[1]]
  byarg <- args$byarg
  if (!is.null(byarg)) {
    stop("wrong tergmLite term used")
  }
  homophily <- args$homophily
  if (!is.null(homophily)) {
    stop("degree homophily argument not supported in tergmLite")
  }
  emptynwstats <- NULL

  if (any(d == 0)) {
    emptynwstats <- rep(0, length(d))
    emptynwstats[d == 0] <- attr(dat$el,'n') # network size
  }

  if (length(d) == 0) {
    return(NULL)
  }
  inputs <- c(d)

  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
  term$maxval <- attr(dat$el[[network]], "n")

  return(term)
}

updateInputs_degree_by_attr <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)

  d <- args[[1]]
  byarg <- args$by
  homophily <- args$homophily
  if (is.null(homophily)) {
    homophily <- FALSE
  } else {
    stop("homophily = TRUE not supported in tergmLite")
  }
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
  }
  if (is.null(byarg)) {
    if (length(d) == 0) {
      return(NULL)
    }
    inputs <- c(d)
  } else {
    if (ncol(du) == 0) {
      return(NULL)
    }
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)) {
    term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
    term$emptynwstats <- emptynwstats

  } else {
    term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
    term$maxval <- attr(dat$el[[network]], "n")
  }

  return(term)
}

updateInputs_absdiff <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)

  attrname <- args[[1]]
  pow <- args$pow

  inputs <- c(pow, dat$attr[[attrname]])
  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_absdiffby <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)
  attrname <- args[[1]]
  byname <- args[[2]]
  offset <- args[[3]]
  values <- args[[4]]
  if (length(values) != 2) {
    stop(paste("\"by\" nodal attribute must be binary, and therefore should have 2 unique values.
               \nVector of values passed was of length ", length(values), sep = ""))
  }
  if (!all(values %in% unique(dat$attr[[byname]]))) {
    stop(paste("Values of binary nodal attribute do not match those in formula term argument:
               \nValues in formula term argument:",
               values[1], values[2], "\nValues of nodal attribute:", unique(dat$attr[[byname]])[1],
               unique(dat$attr[[byname]])[2], sep = " "))
  }
  if (!all(values == c(0, 1))) {
    nodeby <- 1 * (dat$attr[[byname]] == values[2])
  } else {
    nodeby <- dat$attr[[byname]]
  }
  inputs <- c(offset, dat$attr[[attrname]], nodeby)
  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_nodecov <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)
  attrname <- args[[1]]

  f <- args$transform
  nodecov <- dat$attr[[attrname]]

  if (!is.null(f)) {
    inputs <- f(nodecov)
  } else {
    inputs <- nodecov
  }
  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_nodemix <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form,t)
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
  attr(inputs, "ParamsBeforeCov") <- 2 * length(uun)
  term$inputs <- c(2 * length(uun), length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_absdiffnodemix <- function(dat, network, term, t) {

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

  term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)

  return(term)
}

updateInputs_degrange <- function(dat, network, term, t) {

  form <- dat$nwparam[[network]]$formation
  args <- get_formula_term_args_in_formula_env(form, t)

  n <- attributes(dat$el[[network]])$n

  from <- args$from
  to <- args$to
  byarg <- args$by
  homophily <- args$homophily
  to <- ifelse(to == Inf, n + 1, to)
  if (length(to) == 1 && length(from) > 1) {
    to <- rep(to, length(from))
  } else if (length(from) == 1 && length(to) > 1) {
    from <- rep(from, length(to))
  } else if (length(from) != length(to)) {
    stop("The arguments of term degrange must have arguments either of the same
         length, or one of them must have length 1.")
  } else if (any(from >= to)) {
    stop("Term degrange must have from<to.")
  }

  emptynwstats <- NULL
  if (!is.null(byarg)) {
    nodecov <- dat$attr[[byarg]]
    u <- sort(unique(nodecov))
    if (any(is.na(nodecov))) {
      u <- c(u, NA)
    }
    nodecov <- match(nodecov, u)
    if (length(u) == 1) {
      stop("Attribute given to degrange() has only one value",
           call. = FALSE)
    }
  }

  if (!is.null(byarg) && !homophily) {
    lu <- length(u)
    du <- rbind(rep(from, lu), rep(to, lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1, ] == 0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3, du[1, ] == 0]
      for (i in 1:length(tmp)) tmp[i] <- sum(nodecov == tmp[i])
      emptynwstats[du[1, ] == 0] <- tmp
    }
  } else {
    if (any(from == 0)) {
      emptynwstats <- rep(0, length(from))
      emptynwstats[from == 0] <- n
    }
  }
  if (is.null(byarg)) {
    if (length(from) == 0) {
      return(NULL)
    }
    coef.names <- ifelse(to >= n + 1,
                         paste("deg", from, "+", sep = ""),
                         paste("deg", from, "to", to, sep = ""))
    name <- "degrange"
    inputs <- c(rbind(from, to))
  } else if (homophily) {
    if (length(from) == 0) {
      return(NULL)
    }
    coef.names <- ifelse(to >= n + 1,
                         paste("deg", from, "+", ".homophily.", byarg, sep = ""),
                         paste("deg", from, "to", to, ".homophily.", byarg, sep = ""))
    name <- "degrange_w_homophily"
    inputs <- c(rbind(from, to), nodecov)
  } else {
    if (ncol(du) == 0) {
      return(NULL)
    }
    coef.names <- ifelse(du[2, ] >= n + 1,
                         paste("deg", du[1, ], "+.", byarg, u[du[3, ]], sep = ""),
                         paste("deg", du[1, ], "to", du[2, ], ".", byarg,
                               u[du[3, ]], sep = ""))
    name <- "degrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }

  if (!is.null(emptynwstats)) {
    term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
    term$emptynwstats <- emptynwstats

  } else {
    term$inputs <- c(0, length(term$coef.names), length(inputs), inputs)
    term$maxval <- attr(dat$el[[network]], "n")
  }

  return(term)
}
