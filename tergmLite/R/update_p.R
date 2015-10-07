
#' @export
update_p.hiv <- function(dat) {

  mf <- dat$p$model.form
  md <- dat$p$model.diss
  mhf <- dat$p$MHproposal.form
  mhd <- dat$p$MHproposal.diss

  n <- attributes(dat$el)$n
  maxdyads <- choose(n, 2)

  ## 1. Update model.form ##

  # edges
  # inputs <- c(0, 1, 0) # not changed
  mf$terms[[1]]$maxval <- maxdyads

  # concurrent by male
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u)
  ui <- seq(along = u)
  inputs <- c(ui, nodecov)

  outputs <- c(0, 2, length(inputs), inputs)
  mf$terms[[2]]$inputs <- outputs
  mf$terms[[2]]$maxval <- n

  # absdiffby
  nodecov <- dat$attr$age
  nodeby <- dat$attr$male
  offsetVal <- mf$terms[[3]]$inputs[4]
  inputs <- c(offsetVal, nodecov, nodeby)
  mf$terms[[3]]$inputs <- c(0, 1, length(inputs), inputs)

  # nodematch
  nodecov <- dat$attr$male
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov, u, nomatch = length(u) + 1)
  inputs <- nodecov
  mf$terms[[4]]$inputs <- c(0, 1, length(inputs), inputs)

  ## Update combined maxval here
  mf$maxval <- c(maxdyads, rep(n, 2), Inf, Inf)


  ## 2. Update model.diss ##
  md$terms[[1]]$maxval <- maxdyads
  md$maxval <- maxdyads


  ## 3. Update MHproposal.form ##
  mhf$arguments$constraints$bd$attribs <- matrix(rep(mhf$arguments$constraints$bd$attribs[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxout <- matrix(rep(mhf$arguments$constraints$bd$maxout[1], n), ncol = 1)
  mhf$arguments$constraints$bd$maxin <- matrix(rep(n, n), ncol = 1)
  mhf$arguments$constraints$bd$minout <- mhf$arguments$constraints$bd$minin <- matrix(rep(0, n), ncol = 1)


  ## 4. Update MHproposal.diss ##
  mhd$arguments$constraints$bd <- mhf$arguments$constraints$bd


  ## 5. Output ##
  p <- list(model.form = mf, model.diss = md,
            MHproposal.form = mhf, MHproposal.diss = mhd)
  dat$p <- p

  return(dat)
}
