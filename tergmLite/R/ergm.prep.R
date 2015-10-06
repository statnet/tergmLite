
#' @export
ergm_prep <- function(nw, formation, dissolution,
                      coef.form, coef.diss, constraints, control) {

  if (!is.network(nw)) {
    stop("A network object must be given")
  }

  formation <- ergm.update.formula(formation, nw ~ ., from.new = "nw")
  dissolution <- ergm.update.formula(dissolution, nw ~ ., from.new = "nw")

  model.form <- ergm_getmodel(formation, nw, role = "formation")
  if (!missing(coef.form) && coef.length.model(model.form) != length(coef.form)) {
    stop("coef.form has ", length(coef.form), " elements, while the model requires ",
         coef.length.model(model.form), " parameters.")
  }
  model.diss <- ergm_getmodel(dissolution, nw, role = "dissolution")
  if (!missing(coef.diss) && coef.length.model(model.diss) != length(coef.diss)) {
    stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",
         coef.length.model(model.diss), " parameters.")
  }

  MHproposal.form <- MHproposal(constraints,
                                control$MCMC.prop.args.form,
                                nw,
                                weights = control$MCMC.prop.weights.form,
                                class = "f")
  MHproposal.diss <- MHproposal(constraints,
                                control$MCMC.prop.args.diss,
                                nw,
                                weights = control$MCMC.prop.weights.diss,
                                class = "d")

  out <- list(model.form = model.form, model.diss = model.diss,
              MHproposal.form = MHproposal.form, MHproposal.diss = MHproposal.diss)

  return(out)
}
