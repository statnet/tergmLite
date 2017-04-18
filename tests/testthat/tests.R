
test_that("concurrent", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + concurrent,
                target.stats = c(50, 25),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  p <- stergm_prep(nw,
                   est$formation,
                   est$coef.diss$dissolution,
                   est$coef.form,
                   est$coef.diss$coef.adj,
                   est$constraints)

  dat <- list()
  dat$p <- dat$el <- dat$attr <- dat$nwparam <- list()
  dat$p[[1]] <- p
  dat$el[[1]] <- as.edgelist(simulate(est$fit))
  dat$nwparam[[1]] <- est

  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p[[1]], p)

})

test_that("concurrent_by", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 1, 0.5))

  est <- netest(nw = nw,
                formation = ~edges + concurrent(by = "riskg"),
                target.stats = c(50, 20, 10),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  p <- stergm_prep(nw,
                   est$formation,
                   est$coef.diss$dissolution,
                   est$coef.form,
                   est$coef.diss$coef.adj,
                   est$constraints)

  dat <- list()
  dat$p <- dat$el <- dat$attr <- dat$nwparam <- list()
  dat$p[[1]] <- p
  dat$el[[1]] <- as.edgelist(simulate(est$fit))
  dat$attr$riskg <- nw %v% "riskg"
  dat$nwparam[[1]] <- est

  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p[[1]], p)

})

