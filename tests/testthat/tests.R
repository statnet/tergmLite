
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

test_that("degrange", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4),
                target.stats = c(50, 0),
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


  ## from + to args

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4, to = 6),
                target.stats = c(50, 0),
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

test_that("nodefactor single", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  riskg <- rep(1:4, each = 25)
  nw <- set.vertex.attribute(nw, "riskg", riskg)

  est <- netest(nw = nw,
                formation = ~edges + nodefactor("riskg"),
                target.stats = c(50, 25, 25, 25),
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
  dat$attr <- list()
  dat$attr$riskg <- riskg

  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p[[1]], p)

  # p$model.form
  # p$model.form$terms[[2]]
  # dat$p[[1]]$model.form$terms[[2]]
  #
  # length(p$model.form$terms[[2]]$inputs)
  # length(dat$p[[1]]$model.form$terms[[2]]$inputs)

})

test_that("nodefactor interaction", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  riskg <- sample(rep(1:2, each = 50))
  race <- sample(rep(0:1, each = 50))
  nw <- set.vertex.attribute(nw, "riskg", riskg)
  nw <- set.vertex.attribute(nw, "race", race)

  est <- netest(nw = nw,
                formation = ~edges + nodefactor(c("riskg", "race")),
                target.stats = c(50, 25, 25, 25),
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
  dat$attr <- list()
  dat$attr$riskg <- riskg
  dat$attr$race <- race

  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p[[1]], p)

})


test_that("get_formula_term...", {

  form <- ~edges + concurrent(by = "riskg")
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_is(args, "list")
  expect_true(names(args) == "by")
  expect_true(args[[1]] == "riskg")

  # absdiff
  form <- ~edges + absdiff("riskg", 2)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "pow")

  form <- ~edges + absdiff("riskg", pow = 2)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "pow")

  form <- ~edges + absdiff("riskg")
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "pow")

  # absdiffby
  form <- ~edges + absdiffby("riskg", "by", "offset", values = c(0, 1))
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(args[[2]] == "by")
  expect_true(args[[3]] == "offset")
  expect_true(names(args)[4] == "values")

  form <- ~edges + absdiffby("riskg", "by", "offset", c(0, 1))
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(args[[2]] == "by")
  expect_true(args[[3]] == "offset")
  expect_true(names(args)[4] == "values")

  form <- ~edges + absdiffby("riskg", "by", "offset")
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(args[[2]] == "by")
  expect_true(args[[3]] == "offset")
  expect_true(names(args)[4] == "values")

  # nodefactor
  form <- ~edges + nodefactor("riskg", base = 1)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "base")

  form <- ~edges + nodefactor("riskg", 1)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "base")

  form <- ~edges + nodefactor("riskg")
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "base")

  # nodematch
  form <- ~edges + nodematch("riskg", diff = F)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "diff")

  form <- ~edges + nodematch("riskg", F)
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "diff")

  form <- ~edges + nodematch("riskg")
  args <- get_formula_term_args_in_formula_env(form, 2)
  expect_true(args[[1]] == "riskg")
  expect_true(names(args)[2] == "diff")

})
