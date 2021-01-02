

test_that("simulate_network", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {
    # Set seed for reproducibility
    set.seed(1234)
    
    nw <- network_initialize(100)
    nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
    formation <- ~edges + nodefactor("group")
    target.stats <- c(15, 10)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.25)
    init <- init.net(i.num = 10, i.num.g2 = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # Full network structure after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    str(dat, max.level = 1)
    
    # Current network structure
    dat$el[[1]]
    
    # New network structure
    dat$el[[1]] <- simulate_network(state = dat$p[[1]]$state, 
                                    coef = c(dat$nwparam[[1]]$coef.form, 
                                             dat$nwparam[[1]]$coef.diss$coef.adj),
                                    control = dat$control$MCMC_control[[1]],
                                    save.changes = TRUE)$el
    dat$el[[1]]
    
    # Specific changes listed under changes list
    #    (new edges: to = 1; dissolved edges: to = 0):
    attributes(dat$el[[1]])$changes
  }
})


test_that("simulate_ergm", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {    
    # Set seed for reproducibility
    set.seed(1234)
    
    nw <- network_initialize(100)
    nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
    formation <- ~edges + nodefactor("group")
    target.stats <- c(15, 10)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.1)
    init <- init.net(i.num = 10, i.num.g2 = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # Full network structure after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    str(dat, max.level = 1)
    
    # Current network structure
    dat$el[[1]]
    
    # New network structure (all edges are new)
    dat$el[[1]] <- simulate_ergm(state = dat$p[[1]]$state,
                                 coef = dat$nwparam[[1]]$coef.form,
                                 control = dat$control$MCMC_control[[1]])$el
    dat$el[[1]]
  }
})


test_that("init_tergmLite", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    str(dat, max.level = 1)
    
    # networkLite representation used by tergmLite
    str(dat$p, max.level = 3)
    
    # Elements removed are nw (network class object)
    # Elements added are el (edgelist representation of network)...
    dat$el
    
    # ... and p (contains all relevant ERGM structural information for simulation)
    str(dat$p, max.level = 3)
  }
})


test_that("networkLite", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    
    # Conversion to networkLite class format
    nwl <- networkLite(dat$el[[1]], dat$attr)
    nwl
  }
})


test_that("updateModelTermInputs", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {    
    # Set seed for reproducibility
    set.seed(1234)
    
    nw <- network_initialize(100)
    nw <- set_vertex_attribute(nw, "group", rep(1:2, each = 50))
    formation <- ~edges + nodefactor("group")
    target.stats <- c(15, 10)
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3, inf.prob.g2 = 0.1)
    init <- init.net(i.num = 10, i.num.g2 = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # Full network structure after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    str(dat, max.level = 1)
    
    # Examine the network list structure for nodefactor term
    dat$p[[1]]$state$model$terms[[2]]
    
    # inputs vector corresponds to group attribute stored here
    dat$attr$group
    
    # As example of what could happen in EpiModel: randomly reshuffle group
    #   attribute values of 100 nodes
    dat$attr$group <- sample(dat$attr$group)
    dat$attr$group
    
    # Update network list structure
    dat <- updateModelTermInputs(dat)
    
    # Check that network list structure for nodefactor term has been updated
    dat$p[[1]]$state$model$terms[[2]]
  }
})


test_that("add_vertices", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # networkLite representation after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    
    # Check current network size
    attributes(dat$el[[1]])$n
    
    # Add 10 vertices
    dat$el[[1]] <- add_vertices(dat$el[[1]], 10)
    
    # Check new network size
    attributes(dat$el[[1]])$n
  }
})


test_that("delete_vertices", {

  library("EpiModel")
  if (packageVersion("EpiModel") >= "2.0.0") {
    set.seed(12345)
    nw <- network_initialize(100)
    formation <- ~edges
    target.stats <- 50
    coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
    x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
    
    param <- param.net(inf.prob = 0.3)
    init <- init.net(i.num = 10)
    control <- control.net(type = "SI", nsteps = 100, nsims = 5, tergmLite = TRUE)
    
    # Set seed for reproducibility
    set.seed(123456)
    
    # networkLite representation structure after initialization
    dat <- crosscheck.net(x, param, init, control)
    dat <- initialize.net(x, param, init, control)
    
    # Current edges
    head(dat$el[[1]], 20)
    
    # Remove nodes 1 and 2
    nodes.to.delete <- 1:2
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes.to.delete)
    
    # Newly permuted edges
    head(dat$el[[1]], 20)
  }
})
