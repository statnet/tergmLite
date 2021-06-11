
run_checks <- function(nw, est) {
  param <- param.net(inf.prob = 0.3)
  init <- init.net(i.num = 10)
  control <- control.net(type = "SI", nsteps = 100, nsims = 1, tergmLite = TRUE,
                         resimulate.network = TRUE)

  dat <- crosscheck.net(est, param, init, control)
  dat <- initialize.net(est, param, init, control, s = 1)

  attr(dat$el[[1]], "time") <- 0 # so tergmLite time matches tergm time, which defaults to 0
  dat <- updateModelTermInputs(dat, network = 1)

  ## classes and exact set of attributes differ, so we need to massage things 
  ## into a common form before comparison
  el_s <- structure(as.edgelist(matrix(unlist(dat$p[[1]]$state$el), ncol = 2), 
                      n = attr(dat$p[[1]]$state$el, "n"), 
                      directed = attr(dat$p[[1]]$state$el, "directed"), 
                      bipartite = attr(dat$p[[1]]$state$el, "bipartite"), 
                      loops = attr(dat$p[[1]]$state$el, "loops"),
                      vnames = NULL), time = 0)

  ## this checks that the edge sets and basic network attributes are the same
  expect_identical(dat$el[[1]], el_s)

  ## now do a tergm simulation with a networkLite, getting the ergm_state as output
  nw <- networkLite(dat$el[[1]], dat$attr)
  set.seed(0)
  es_t_n <- simulate(nw ~ Form(est$formation) + Diss(est$coef.diss$dissolution), coef = c(est$coef.form, est$coef.diss$coef.adj), output="ergm_state", dynamic=TRUE, control=dat$control$mcmc.control[[1]])

  ## now do a tergm simulation with a networkLite, getting the ergm_state as output
  nwL <- networkLite(dat$el[[1]], dat$attr)
  set.seed(0)
  es_t <- simulate(nwL ~ Form(est$formation) + Diss(est$coef.diss$dissolution), coef = c(est$coef.form, est$coef.diss$coef.adj), output="ergm_state", dynamic=TRUE, control=dat$control$mcmc.control[[1]])
  
  expect_equal(es_t_n, es_t)

  ## do an equivalent tergmLite simulation
  nwparam <- EpiModel::get_nwparam(dat, network = 1)
  set.seed(0)
  es_tL <- simulate_network(dat$p[[1]]$state, coef=c(nwparam$coef.form, nwparam$coef.diss$coef.adj), control=dat$control$mcmc.control[[1]], save.changes=TRUE)$state

  ## the two output states should be equal, up to some minor issues with attributes:
  
  ## force common set and ordering of network attributes here
  es_t$nw0$gal <- es_t$nw0$gal[names(es_tL$nw0$gal)]

  ## tergm state has some additional attributes that tergmLite state does not; remove them here
  attributes(es_t) <- attributes(es_t)[c("names", "class")]
  
  ## call fields will not match
  for(j in seq_along(es_t$model$terms)) {
    es_t$model$terms[[j]]$call <- NULL
  }

  for(j in seq_along(es_tL$model$terms)) {
    es_tL$model$terms[[j]]$call <- NULL
  }

  
  ## the two states should now be equal (not identical due to formulas)
  expect_equal(es_tL, es_t)
}

test_that("concurrent", {

  library("EpiModel")
  nw <- network_initialize(100)

  est <- netest(nw = nw,
                formation = ~edges + concurrent,
                target.stats = c(50, 25),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)
  
})

test_that("concurrent_by", {

  library("EpiModel")
  nw <- network_initialize(100)
  nw <- set_vertex_attribute(nw, "riskg", rbinom(100, 1, 0.5))

  est <- netest(nw = nw,
                formation = ~edges + concurrent(by = "riskg"),
                target.stats = c(50, 20, 10),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("degree, single", {

  library("EpiModel")
  nw <- network_initialize(100)

  est <- netest(nw = nw,
                formation = ~edges + degree(1),
                target.stats = c(50, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("degree, multiple", {

  library("EpiModel")
  nw <- network_initialize(100)

  est <- netest(nw = nw,
                formation = ~edges + degree(1:2),
                target.stats = c(50, 35, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("degree_by_attr", {

  library("EpiModel")
  nw <- network_initialize(100)
  nw <- set_vertex_attribute(nw, "riskg", sample(rep(0:1, each = 50)))

  est <- netest(nw = nw,
                formation = ~edges + degree(1, by = "riskg"),
                target.stats = c(50, 15, 20),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("degrangefrom", {

  library("EpiModel")
  nw <- network_initialize(100)

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("degrangefromto", {

  library("EpiModel")
  nw <- network_initialize(100)

  est <- netest(nw = nw,
                formation = ~edges + degrange(from = 4, to = 6),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("nodecov formula", {

  library("EpiModel")
  nw <- network_initialize(100)
  risk <- runif(100)
  nw <- set_vertex_attribute(nw, "risk", risk)

  est <- netest(nw = nw,
                formation = ~edges + nodecov(~risk^2),
                target.stats = c(50, 65),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("nodecov function", {

  library("EpiModel")
  nw <- network_initialize(100)
  risk <- runif(100)
  nw <- set_vertex_attribute(nw, "risk", risk)
 
  est <- netest(nw = nw,
                formation = ~edges + nodecov(function(x) exp(1 + log((x %v% "risk")^2))),
                target.stats = c(50, 200),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("nodefactor single", {

  library("EpiModel")
  nw <- network_initialize(100)
  riskg <- rep(1:4, each = 25)
  nw <- set_vertex_attribute(nw, "riskg", riskg)
 
  est <- netest(nw = nw,
                formation = ~edges + nodefactor("riskg"),
                target.stats = c(50, 25, 25, 25),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})


test_that("nodefactor interaction", {

  library("EpiModel")
  nw <- network_initialize(100)
  riskg <- sample(rep(1:2, each = 50))
  race <- sample(rep(0:1, each = 50))
  nw <- set_vertex_attribute(nw, "riskg", riskg)
  nw <- set_vertex_attribute(nw, "race", race)
 
  est <- netest(nw = nw,
                formation = ~edges + nodefactor(c("riskg", "race")),
                target.stats = c(50, 25, 25, 25),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("nodemix levels", {

  library("EpiModel")
  nw <- network_initialize(200)
  race <- sample(rep(letters[1:4], each = 50))
  nw <- set_vertex_attribute(nw, "race", race)
 
  est <- netest(nw = nw,
                formation = ~edges + nodemix("race", levels = c("a", "b", "d"), levels2=-(2:3)),
                target.stats = c(200, 12, 25, 25, 12),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("triangle", {

  library("EpiModel")
  nw <- network_initialize(100)
 
  est <- netest(nw = nw,
                formation = ~edges + triangle,
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("triangle_attr", {

  library("EpiModel")
  nw <- network_initialize(100)
  nw <- set_vertex_attribute(nw, "riskg", rbinom(100, 1, 0.2))
 
  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg"),
                target.stats = c(50, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)
  
})

test_that("triangle_attrdiff", {

  library("EpiModel")
  nw <- network_initialize(100)
  nw <- set_vertex_attribute(nw, "riskg", rbinom(100, 1, 0.5))

  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg", diff=TRUE),
                target.stats = c(50, 0, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})

test_that("triangle_attrdifflevels", {

  library("EpiModel")
  nw <- network_initialize(100)
  nw <- set_vertex_attribute(nw, "riskg", rbinom(100, 2, 0.1))

  est <- netest(nw = nw,
                formation = ~edges + triangle(attr = "riskg", diff=TRUE, levels=c(1,2)),
                target.stats = c(50, 0, 0),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)

})


test_that("gwesp_true", {

  library("EpiModel")
  nw <- network_initialize(100)
 
  est <- netest(nw = nw,
                formation = ~edges + gwesp(0, fixed=TRUE),
                target.stats = c(50, 2),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("gwesp_truedecay", {

  library("EpiModel")
  nw <- network_initialize(100)
 
  est <- netest(nw = nw,
                formation = ~edges + gwesp(decay=0.1, fixed=TRUE),
                target.stats = c(50, 10),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))
 
  run_checks(nw, est)

})

test_that("absdiff", {

  library("EpiModel")
  nw <- network_initialize(100)
  age <- sample(15:65, 100, TRUE)
  nw <- set_vertex_attribute(nw, "age", age)

  est <- netest(nw = nw,
                formation = ~edges + absdiff("age"),
                target.stats = c(50, 100),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)
  
})

test_that("absdiffby", {
  
  library("EpiModel")
  nw <- network_initialize(100)
  age <- sample(15:65, 100, TRUE)
  sex <- rbinom(100, 1, 0.5)
  nw <- set_vertex_attribute(nw, "age", age)
  nw <- set_vertex_attribute(nw, "sex", sex)
  
  est <- netest(nw = nw,
                formation = ~edges + absdiffby("age", "sex", 2),
                target.stats = c(50, 100),
                coef.diss = dissolution_coefs(~offset(edges), duration = 100))

  run_checks(nw, est)
  
})
