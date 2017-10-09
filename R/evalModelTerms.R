
#' @title Evaluate ERGM Formula Terms
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
  if (tname == "nodefactor" & length(outlist) == 1) {
    outlist$base <- 1
  }
  if (tname == "nodefactor" & length(outlist) == 2 & is.null(names(outlist)[2])) {
    names(outlist)[2] <- "base"
  }

  # Set default pow to 1
  if (tname == "absdiff" & length(outlist) == 1) {
    outlist$pow <- 1
  }
  if (tname == "absdiff" & length(outlist) == 2 & is.null(names(outlist)[2])) {
    names(outlist)[2] <- "pow"
  }

  # Set default values in absdiffby to 0, 1
  if (tname == "absdiffby" & length(outlist) == 3) {
    outlist$values <- c(0, 1)
  }
  if (tname == "absdiffby" & length(outlist) == 4 & is.null(names(outlist)[4])) {
    names(outlist)[4] <- "values"
  }

  # set default fixed argument to FALSE
  # also needs default keep argument, set to NULL
  if (tname == "nodematch" & length(outlist) == 1) {
    outlist$diff <- FALSE
  }
  if (tname == "nodematch" & length(outlist) == 2 & is.null(names(outlist)[2])) {
    names(outlist)[2] <- "diff"
  }

  if (tname == "degrange" & length(outlist) == 1 & names(outlist) == "from") {
    outlist$to <- Inf
  }

  if (tname == "degrange") {
    if (!("from" %in% names(outlist))) {
      outlist$from <- 0
    }
    if (!("to" %in% names(outlist))) {
      outlist$to <- Inf
    }
    if (!("homophily" %in% names(outlist))) {
      outlist$homophily <- FALSE
    }
  }

  return(outlist)
}
