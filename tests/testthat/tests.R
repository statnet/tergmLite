
test_that("concurrent", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + concurrent,
                target.stats = c(50, 25),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("concurrent_by", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 1, 0.5))

  est <- netest(nw = nw,
                formation = ~edges + concurrent(by = "riskg"),
                target.stats = c(50, 20, 10),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("degree, single", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + degree(1),
                target.stats = c(50, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("degree, multiple", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + degree(1:2),
                target.stats = c(50, 35, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("degree_by_attr", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", sample(rep(0:1, each = 50)))

  est <- netest(nw = nw,
                formation = ~edges + degree(1, by = "riskg"),
                target.stats = c(50, 15, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("degrange", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)


  ## from + to args

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4, to = 6),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("nodecov formula", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  risk <- runif(100)
  nw <- set.vertex.attribute(nw, "risk", risk)

  est <- netest(nw = nw,
                formation = ~edges + nodecov(~risk^2),
                target.stats = c(50, 65),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("nodecov function", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  risk <- runif(100)
  nw <- set.vertex.attribute(nw, "risk", risk)

  est <- netest(nw = nw,
                formation = ~edges + nodecov(function(x) exp(1 + log((x %v% "risk")^2))),
                target.stats = c(50, 200),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

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

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

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

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("nodemix levels", {

  library("EpiModel")
  nw <- network.initialize(200, directed = FALSE)
  race <- sample(rep(letters[1:4], each = 50))
  nw <- set.vertex.attribute(nw, "race", race)

  est <- netest(nw = nw,
                formation = ~edges + nodemix("race", levels = c("a", "b", "d"), levels2=-(2:3)),
                target.stats = c(200, 12, 25, 25, 12),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("triangle", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + triangle,
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("triangle_attr", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 1, 0.2))

  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg"),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("triangle_attrdiff", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 1, 0.5))

  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg", diff=TRUE),
                target.stats = c(50, 0, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("triangle_attrdifflevels", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)
  nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 2, 0.1))

  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg", diff=TRUE, levels=c(1,2)),
                target.stats = c(50, 0, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})


test_that("gwesp_true", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + gwesp(fixed=TRUE),
                target.stats = c(50, 2),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})

test_that("gwesp_truedecay", {

  library("EpiModel")
  nw <- network.initialize(100, directed = FALSE)

  est <- netest(nw = nw,
                formation = ~edges + gwesp(decay=0.8, fixed=TRUE),
                target.stats = c(50, 1.1),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, depend = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control)
  dat <- init_tergmLite(dat)

  p <- dat$p
  dat <- updateModelTermInputs(dat, network = 1)

  expect_identical(dat$p, p)

})
