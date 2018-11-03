

#' Converts Network Object, Formation and Dissolution Formulas, Formation and
#' Dissolution Coefficients, and Control Settings to a Data-Thin List Format
#' for STERGM Resimulation
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
#' @export
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

  out <- list(model.form = model.form, model.diss = model.diss,
              MHproposal.form = MHproposal.form, MHproposal.diss = MHproposal.diss)

  return(out)
}


#' Converts Network Object, Formation Formulas, Formation Coefficients, and
#' Control Settings to a Data-Thin List Format for ERGM Resimulation
#'
#' @param nw An object of class \code{network}
#' @param formation Right-hand sided formation formula
#' @param coef Vector of coefficients associated with the formation formula
#' @param constraints Constraints for the formation model (only \code{bd})
#'                    constraints currently supported.
#' @param control Control settings passed to \code{ergm::control.simulate.ergm}
#'
#' @export
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

  out <- list(model.form = m, MHproposal = MHproposal)
  return(out)
}


# used for testing, not needed for simulation

# ergm_getmodel <- function(formula,
#                           nw,
#                           response = NULL,
#                           silent = FALSE,
#                           role = "static",
#                           ...) {
#
#   if ((dc <- data.class(formula)) != "formula") {
#     stop(paste("Invalid formula of class ", dc), call. = FALSE)
#   }
#
#   if (formula[[1]] != "~") {
#     stop("Formula must be of the form 'network ~ model'.", call. = FALSE)
#   }
#
#   if (length(formula) < 3) {
#     stop(paste("No model specified for network ", formula[[2]]), call. = FALSE)
#   }
#
#   v <- term.list.formula(formula[[3]])
#   formula.env <- environment(formula)
#   model <- structure(list(formula = formula,
#                           coef.names = NULL,
#                           offset = NULL,
#                           terms = NULL,
#                           networkstats.0 = NULL,
#                           etamap = NULL),
#                      class = "model.ergm")
#
#   termroot <- if (is.null(response)) {
#     "InitErgm"
#   } else {
#     "InitWtErgm"
#   }
#
#   for (i in 1:length(v)) {
#     if (is.call(v[[i]])) {
#       if (v[[i]][[1]] == "offset") {
#         if (length(v[[i]][[2]]) <= 1) {
#           v[[i]] <- as.call(v[[i]][2])
#         } else {
#           v[[i]] <- as.call(v[[i]][[2]])
#         }
#         model$offset <- c(model$offset, TRUE)
#       } else {
#         model$offset <- c(model$offset, FALSE)
#       }
#       args = v[[i]]
#       args[[1]] = as.name("list")
#       fname <- paste(termroot, "Term.", v[[i]][[1]], sep = "")
#       newInitErgm <- exists(fname, envir = formula.env,
#                             mode = "function")
#       v[[i]] <- call(ifelse(newInitErgm, fname, paste(termroot,
#                                                       ".", v[[i]][[1]], sep = "")))
#     } else {
#       fname <- paste(termroot, "Term.", v[[i]], sep = "")
#       newInitErgm <- exists(fname, envir = formula.env,
#                             mode = "function")
#       v[[i]] <- call(ifelse(newInitErgm, fname, paste(termroot,
#                                                       ".", v[[i]], sep = "")))
#       model$offset <- c(model$offset, FALSE)
#       args = list()
#     }
#
#     if (!newInitErgm) {
#       v[[i]][[2]] <- nw
#       names(v[[i]])[2] <- ""
#       v[[i]][[3]] <- model
#       names(v[[i]])[3] <- ""
#       v[[i]][[4]] <- args
#       dotdotdot <- c(if (!is.null(response)) list(response = response),
#                      list(role = role), list(...))
#       for (j in seq_along(dotdotdot)) {
#         if (is.null(dotdotdot[[j]]))
#           next
#         v[[i]][[4 + j]] <- dotdotdot[[j]]
#         names(v[[i]])[4 + j] <- names(dotdotdot)[j]
#       }
#       if (!exists(as.character(v[[i]][[1]]), envir = formula.env, mode = "function")) {
#         stop("The term ", substring(as.character(v[[i]][[1]]),
#                                     first = nchar(termroot) + 2),
#              " does not exist for this type of ERGM. Are you sure you have the right name?\n",
#              call. = FALSE)
#       }
#       if (silent) {
#         silentwarnings <- capture.output(model <- eval(v[[i]], formula.env))
#       } else {
#         model <- eval(v[[i]], formula.env)
#       }
#       if (is.null(model$terms[[length(model$terms)]]$pkgname)) {
#         model$terms[[length(model$terms)]]$pkgname <-
#           ergm:::which.package.InitFunction(v[[i]][[1]], formula.env)
#       }
#     } else {
#       v[[i]][[2]] <- nw
#       names(v[[i]])[2] <- ""
#       v[[i]][[3]] <- args
#       names(v[[i]])[3] <- ""
#       dotdotdot <- c(if (!is.null(response)) list(response = response), list(role = role), list(...))
#       for (j in seq_along(dotdotdot)) {
#         if (is.null(dotdotdot[[j]]))
#           next
#         v[[i]][[3 + j]] <- dotdotdot[[j]]
#         names(v[[i]])[3 + j] <- names(dotdotdot)[j]
#       }
#       outlist <- eval(v[[i]], formula.env)
#       if (is.null(outlist$pkgname)) {
#         outlist$pkgname <- ergm:::which.package.InitFunction(v[[i]][[1]], formula.env)
#       }
#       model <- updatemodel_ErgmTerm(model, outlist)
#     }
#   }
#
#   model$etamap <- ergm.etamap(model)
#   ergm:::ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
#
#   class(model) <- "ergm.model"
#   return(model)
# }

# updatemodel_ErgmTerm <- function(model, outlist) {
#   if (!is.null(outlist)) {
#     model$coef.names <- c(model$coef.names, outlist$coef.names)
#     termnumber <- 1 + length(model$terms)
#     tmp <- attr(outlist$inputs, "ParamsBeforeCov")
#     outlist$inputs <- c(ifelse(is.null(tmp), 0, tmp), length(outlist$coef.names),
#                         length(outlist$inputs), outlist$inputs)
#     model$minval <- c(model$minval, rep(if (!is.null(outlist$minval)) outlist$minval else -Inf,
#                                         length.out = length(outlist$coef.names)))
#     model$maxval <- c(model$maxval, rep(if (!is.null(outlist$maxval)) outlist$maxval else +Inf,
#                                         length.out = length(outlist$coef.names)))
#     model$duration <- c(model$duration, if (!is.null(outlist$duration)) outlist$duration else FALSE)
#     model$terms[[termnumber]] <- outlist
#   }
#
#   return(model)
# }
