test_that("network and networkLite simulate identically in ergm", {
  require(ergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
  
  coef <- c(-4, 1, 1.5, 0.5, -1, 0.5)
  
  set.seed(0)
  nw_1 <- simulate(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  set.seed(0)
  nwL_1 <- simulate(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)
  
  set.seed(0)
  nw_2 <- simulate(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  set.seed(0)
  nwL_2 <- simulate(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), coef = coef, output = "network", dynamic = FALSE)
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)  
})

test_that("network and networkLite simulate identically in san", {
  require(ergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
    
  set.seed(0)
  nw_1 <- san(nw ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
  set.seed(0)
  nwL_1 <- san(nwL ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(1000, 500, 300, 200, 600, 1500))
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)
  
  set.seed(0)
  nw_2 <- san(nw_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
  set.seed(0)
  nwL_2 <- san(nwL_1 ~ edges + nodefactor("a") + nodecov(~b^2 + b), target.stats = c(800, 400, 200, 100, 600, 1200))
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)
})

test_that("network and networkLite simulate identically in tergm", {
  require(tergm)

  nw <- network.initialize(100, dir = FALSE)
  nw %v% "a" <- rep(letters[1:5], each = 20)
  nw %v% "b" <- runif(100)
  
  nwL <- as.networkLite(nw)
  
  coef <- c(-4, 1, 1.5, 0.5, -1, 0.5, 3)
  
  set.seed(0)
  nw_1 <- simulate(nw ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_1 <- simulate(nwL ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_1)), unclass(as.edgelist(nwL_1)), check.attributes = FALSE)
  expect_identical(nw_1 %n% "lasttoggle", nwL_1 %n% "lasttoggle")
  expect_identical(nw_1 %n% "time", nwL_1 %n% "time")
  
  set.seed(0)
  nw_2 <- simulate(nw_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_2 <- simulate(nwL_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_2)), unclass(as.edgelist(nwL_2)), check.attributes = FALSE)  
  expect_identical(nw_2 %n% "lasttoggle", nwL_2 %n% "lasttoggle")
  expect_identical(nw_2 %n% "time", nwL_2 %n% "time")

  set.seed(0)
  nw_3 <- simulate(nw_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  set.seed(0)
  nwL_3 <- simulate(nwL_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) + Diss(~edges), coef = coef, output = "final", dynamic = TRUE)
  
  expect_equal(unclass(as.edgelist(nw_3)), unclass(as.edgelist(nwL_3)), check.attributes = FALSE)  
  expect_identical(nw_3 %n% "lasttoggle", nwL_3 %n% "lasttoggle")
  expect_identical(nw_3 %n% "time", nwL_3 %n% "time")
})
