test_that("manual and tergmLite dynamic simulations produce identical results", {
  require(tergm)
  require(EpiModel)

  age_vals <- 1:20
  race_vals <- c("A","B","C","D","E")
  sex_vals <- c("M","F")
  
  set.seed(0)
  nw <- network.initialize(1000, dir = FALSE)
  nw %v% "age" <- sample(age_vals, 1000, TRUE)
  nw %v% "race" <- sample(race_vals, 1000, TRUE)
  nw %v% "sex" <- sample(sex_vals, 1000, TRUE)
  
  formation <- ~edges + nodecov(~age) + nodefactor(~race)
  dissolution <- ~offset(edges)

  ## obtained from simulating previously chosen coefficients
  target_stats <- c(835.168, 18817.544, 378.888, 296.758, 320.370, 409.382)

  diss_coefs <- dissolution_coefs(dissolution, duration = 1 + exp(3))    
  
  ff <- ~Form(formation) + Diss(dissolution)

  pmat <- matrix(3 + runif(25), 5, 5)
  pmat <- pmat + t(pmat)
  
  constraints <- ~bd(maxout = 3) + blocks(~sex, levels2 = diag(TRUE, 2))

  ff_m <- ~edges + mean.age + degree(0:3) + degrange(4) + nodematch("sex")

  # set some arbitrary, non-default control to ensure it gets propagated correctly
  control <- control.simulate.network.tergm(MCMC.burnin.min = 54321, MCMC.burnin.max = 123456, MCMC.prop = ~Strat(~race, pmat=pmat) + discord + TNT)

  update_nw <- function(nw, nodes_to_remove, nodes_to_add, age_vals, race_vals) {
    el <- as.edgelist(nw)
    lt <- nw %n% "lasttoggle"
    time <- nw %n% "time"
    
    age <- nw %v% "age"
    race <- nw %v% "race"
    sex <- nw %v% "sex"
    
    new_indices <- seq_len(network.size(nw))
    new_indices[nodes_to_remove] <- 0
    new_indices <- new_indices - cumsum(new_indices == 0)
    
    nw <- network.initialize(network.size(nw) + nodes_to_add - length(nodes_to_remove), directed = FALSE)
    nw %v% "age" <- c(age[-nodes_to_remove], sample(age_vals, nodes_to_add, TRUE))
    nw %v% "race" <- c(race[-nodes_to_remove], sample(race_vals, nodes_to_add, TRUE))
    nw %v% "sex" <- c(sex[-nodes_to_remove], sample(sex_vals, nodes_to_add, TRUE))
    nw %n% "time" <- time
    
    ## need to update el and lt for removed nodes
    el <- el[!(el[,1] %in% nodes_to_remove) & !(el[,2] %in% nodes_to_remove),,drop=FALSE]
    lt <- lt[!(lt[,1] %in% nodes_to_remove) & !(lt[,2] %in% nodes_to_remove),,drop=FALSE]
    el <- matrix(new_indices[c(el)], ncol = 2)
    lt <- matrix(c(new_indices[c(lt[,1:2])],lt[,3]), ncol = 3)
    
    nw[el] <- 1L
    nw %n% "lasttoggle" <- lt
    
    nw
  }

  set.seed(0)
  nws <- nw

  nw_ergm <- ergm(nw ~ edges + nodecov(~age) + nodefactor(~race), target.stats = target_stats, constraints = constraints, eval.loglik = FALSE, control = list(init.method="MPLE"))
  nw_coef <- c(coef(nw_ergm)[1] - diss_coefs$coef.form.corr, coef(nw_ergm)[-1], diss_coefs$coef.adj)
  
  nw_el <- list()
  nw_lt <- list()
  nw_time <- list()
  
  nw_summstats <- NULL
  
  nw_summstats <- rbind(nw_summstats, c(summary(ff, basis=nws), summary(ff_m, basis=nws)))

  nws <- simulate(ff, basis = nws, coef = nw_coef, constraints = constraints, control = control, output = "final", dynamic = TRUE)

  nw_el[[length(nw_el) + 1]] <- as.edgelist(nws)
  nw_lt[[length(nw_lt) + 1]] <- nws %n% "lasttoggle"
  nw_time[[length(nw_time) + 1]] <- nws %n% "time"
  
  for(j in seq_len(10)) {
    nws %v% "age" <- nws %v% "age" + 1

    nodes_to_remove <- sample(seq_len(network.size(nws)), rpois(1, 25), FALSE)
    nodes_to_add <- rpois(1,30)
    
    nws <- update_nw(nws, nodes_to_remove, nodes_to_add, age_vals, race_vals)

    nw_summstats <- rbind(nw_summstats, c(summary(ff, basis=nws), summary(ff_m, basis=nws)))
        
    nws <- simulate(ff, basis = nws, coef = nw_coef, constraints = constraints, control = control, output = "final", dynamic = TRUE)

    nw_el[[length(nw_el) + 1]] <- as.edgelist(nws)
    nw_lt[[length(nw_lt) + 1]] <- nws %n% "lasttoggle"
    nw_time[[length(nw_time) + 1]] <- nws %n% "time"
  }


  ## now do the same thing via tergmLite, and test for identical el, lt, and time  
  update_dat <- function(dat, nodes_to_remove, nodes_to_add, age_vals, race_vals) {
    dat$attr$age <- c(dat$attr$age[-nodes_to_remove], sample(age_vals, nodes_to_add, TRUE))
    dat$attr$race <- c(dat$attr$race[-nodes_to_remove], sample(race_vals, nodes_to_add, TRUE))
    dat$attr$sex <- c(dat$attr$sex[-nodes_to_remove], sample(sex_vals, nodes_to_add, TRUE))
    
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes_to_remove)
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", nodes_to_remove)
    dat$el[[1]] <- add_vertices(dat$el[[1]], nodes_to_add)

    dat      
  }


  dat_el <- list()
  dat_lt <- list()
  dat_time <- list()

  dat_summstats <- NULL

  set.seed(0)
  nwL <- as.networkLite(nw)

  nwL_ergm <- ergm(nwL ~ edges + nodecov(~age) + nodefactor(~race), target.stats = target_stats, constraints = constraints, eval.loglik = FALSE, control = list(init.method="MPLE"))
  nwL_coef <- c(coef(nwL_ergm)[1] - diss_coefs$coef.form.corr, coef(nwL_ergm)[-1])

  
  dat <- list(nw = list(nw),
              attr = list(age = nw %v% "age", 
                          race = nw %v% "race",
                          sex = nw %v% "sex"),
              nwparam = list(list(formation = formation,
                                  coef.form = nwL_coef,
                                  coef.diss = diss_coefs,
                                  constraints = constraints)),
              control = list(track_duration = TRUE, MCMC_control = list(control)))

  dat <- init_tergmLite(dat)
  dat <- updateModelTermInputs(dat)

  nwL <- networkLite(dat$el[[1]], dat$attr)
  nwL %n% "time" <- dat$p[[1]]$state$nw0 %n% "time"
  nwL %n% "lasttoggle" <- dat$p[[1]]$state$nw0 %n% "lasttoggle"
  
  dat_summstats <- rbind(dat_summstats, c(summary(ff, basis=nwL), summary(ff_m, basis=nwL)))

  nwparam <- dat$nwparam[[1]]
  rv <- tergmLite::simulate_network(state = dat$p[[1]]$state,
                                    coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
                                    control = dat$control$MCMC_control[[1]],
                                    save.changes = TRUE)

  dat$el[[1]] <- rv$el
  dat$p[[1]]$state$nw0 %n% "time" <- rv$state$nw0 %n% "time"  
  dat$p[[1]]$state$nw0 %n% "lasttoggle" <- rv$state$nw0 %n% "lasttoggle"

  dat_el[[length(dat_el) + 1]] <- rv$el
  dat_lt[[length(dat_lt) + 1]] <- rv$state$nw0 %n% "lasttoggle"
  dat_time[[length(dat_time) + 1]] <- rv$state$nw0 %n% "time"

  for(j in seq_len(10)) {
    dat$attr$age <- dat$attr$age + 1

    nodes_to_remove <- sample(seq_len(attr(rv$el, "n")), rpois(1, 25), FALSE)
    nodes_to_add <- rpois(1,30)

    dat <- update_dat(dat, nodes_to_remove, nodes_to_add, age_vals, race_vals)
    
    dat <- updateModelTermInputs(dat)

    nwL <- networkLite(dat$el[[1]], dat$attr)
    nwL %n% "time" <- dat$p[[1]]$state$nw0 %n% "time"
    nwL %n% "lasttoggle" <- dat$p[[1]]$state$nw0 %n% "lasttoggle"

    dat_summstats <- rbind(dat_summstats, c(summary(ff, basis=nwL), summary(ff_m, basis=nwL)))
    
    nwparam <- dat$nwparam[[1]]
    rv <- tergmLite::simulate_network(state = dat$p[[1]]$state,
                                      coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
                                      control = dat$control$MCMC_control[[1]],
                                      save.changes = TRUE)
  
    dat$el[[1]] <- rv$el
    dat$p[[1]]$state$nw0 %n% "time" <- rv$state$nw0 %n% "time"  
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- rv$state$nw0 %n% "lasttoggle"
  
    dat_el[[length(dat_el) + 1]] <- rv$el
    dat_lt[[length(dat_lt) + 1]] <- rv$state$nw0 %n% "lasttoggle"
    dat_time[[length(dat_time) + 1]] <- rv$state$nw0 %n% "time"
  }

  expect_equal(lapply(nw_el, unclass), lapply(dat_el, unclass), check.attributes = FALSE)
  expect_equal(nw_lt, dat_lt, check.attributes = FALSE)
  expect_equal(nw_time, dat_time, check.attributes = FALSE)
  
  
  
  ## now do the same thing making use of EpiModel
  set.seed(0)
  nwL <- as.networkLite(nw)
  netest_ergm <- netest(nwL, formation, target.stats = target_stats, coef.diss = diss_coefs, constraints = constraints, set.control.ergm = list(init.method = "MPLE"))
  
  init_dat <- function(x, param, init, control, s) {
    nwL <- x$fit$network
    nwL[,] <- FALSE
    
    dat <- list(nw = list(nwL),
                attr = list(age = nwL %v% "age",
                            race = nwL %v% "race",
                            sex = nwL %v% "sex",
                            active = rep(TRUE, network.size(nwL))),
                nwparam = list(list(formation = x$formation,
                                    coef.form = x$coef.form,
                                    coef.diss = x$coef.diss,
                                    constraints = x$constraints)),
                param = list(groups = Inf), # hack
                epi = list(name = matrix(2,2,2)), # hack
                control = control,
                edgelist = list(),
                lasttoggle = list(),
                time = list())
  
    dat <- init_tergmLite(dat)
    
    dat
  }

  update_dat <- function(dat, at) {
    dat$time[[length(dat$time) + 1]] <- dat$p[[1]]$state$nw0 %n% "time"
    dat$lasttoggle[[length(dat$lasttoggle) + 1]] <- dat$p[[1]]$state$nw0 %n% "lasttoggle"
    dat$edgelist[[length(dat$edgelist) + 1]] <- dat$el[[1]]
    
    dat$attr$age <- dat$attr$age + 1
    
    nodes_to_remove <- sample(seq_len(attr(dat$el[[1]], "n")), rpois(1, 25), FALSE)
    nodes_to_add <- rpois(1,30)

    dat$attr$age <- c(dat$attr$age[-nodes_to_remove], sample(age_vals, nodes_to_add, TRUE))
    dat$attr$race <- c(dat$attr$race[-nodes_to_remove], sample(race_vals, nodes_to_add, TRUE))
    dat$attr$sex <- c(dat$attr$sex[-nodes_to_remove], sample(sex_vals, nodes_to_add, TRUE))
    dat$attr$active <- c(dat$attr$active[-nodes_to_remove], rep(TRUE, nodes_to_add))
    
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes_to_remove)
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", nodes_to_remove)
    dat$el[[1]] <- add_vertices(dat$el[[1]], nodes_to_add)

    dat      
  }

  netsim_control <- control.net(nsteps = 12, 
                                initialize.FUN = init_dat,
                                nwupdate.FUN = update_dat,
                                prevalence.FUN = NULL,
                                verbose.FUN = NULL,
                                resimulate.network = TRUE, 
                                tergmLite = TRUE,
                                skip.check = TRUE, 
                                track_duration = TRUE,
                                extract.summary.stats = TRUE, 
                                save.other = c("edgelist", "lasttoggle", "time"),
                                monitors = list(ff_m),
                                MCMC_control = list(control))

  sim <- netsim(netest_ergm, NULL, NULL, netsim_control)
  
  expect_equal(lapply(nw_el, unclass), lapply(sim$edgelist[[1]], unclass), check.attributes = FALSE)
  expect_equal(nw_lt, sim$lasttoggle[[1]], check.attributes = FALSE)
  expect_equal(nw_time, sim$time[[1]])  

  expect_identical(nw_summstats, dat_summstats)
  expect_identical(nw_summstats, sim$stats$summstats[[1]][[1]])
})
