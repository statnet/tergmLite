test_that("manual and tergmLite dynamic simulations produce identical results", {
  require(tergm)

  update_nw <- function(nw, nodes_to_remove, nodes_to_add, age_vals, race_vals) {
    el <- as.edgelist(nw)
    lt <- nw %n% "lasttoggle"
    time <- nw %n% "time"
    
    age <- nw %v% "age"
    race <- nw %v% "race"
    
    new_indices <- seq_len(network.size(nw))
    new_indices[nodes_to_remove] <- 0
    new_indices <- new_indices - cumsum(new_indices == 0)
    
    nw <- network.initialize(network.size(nw) + nodes_to_add - length(nodes_to_remove), directed = FALSE)
    nw %v% "age" <- c(age[-nodes_to_remove], sample(age_vals, nodes_to_add, TRUE))
    nw %v% "race" <- c(race[-nodes_to_remove], sample(race_vals, nodes_to_add, TRUE))
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

  nw_el <- list()
  nw_lt <- list()
  nw_time <- list()
  
  age_vals <- 1:20
  race_vals <- c("A","B","C","D","E")
  
  set.seed(0)
  nw <- network.initialize(1000, dir = FALSE)
  nw %v% "age" <- sample(age_vals, 1000, TRUE)
  nw %v% "race" <- sample(race_vals, 1000, TRUE)

  ff <- ~Form(~edges + nodecov(~age) + nodefactor(~race)) + Diss(~edges)
  
  coef <- c(-10, 0.05, 0.5, 0.3, 0.4, 0.9, 3)

  nw <- simulate(ff, basis = nw, coef = coef, output = "final", dynamic = TRUE)

  nw_el[[length(nw_el) + 1]] <- as.edgelist(nw)
  nw_lt[[length(nw_lt) + 1]] <- nw %n% "lasttoggle"
  nw_time[[length(nw_time) + 1]] <- nw %n% "time"
  
  for(j in seq_len(10)) {
    nw %v% "age" <- nw %v% "age" + 1

    nodes_to_remove <- sample(seq_len(network.size(nw)), rpois(1, 25), FALSE)
    nodes_to_add <- rpois(1,30)
    
    nw <- update_nw(nw, nodes_to_remove, nodes_to_add, age_vals, race_vals)
        
    nw <- simulate(ff, basis = nw, coef = coef, output = "final", dynamic = TRUE)

    nw_el[[length(nw_el) + 1]] <- as.edgelist(nw)
    nw_lt[[length(nw_lt) + 1]] <- nw %n% "lasttoggle"
    nw_time[[length(nw_time) + 1]] <- nw %n% "time"
  }


  ## now do the same thing via tergmLite, and test for identical el, lt, and time  
  update_dat <- function(dat, nodes_to_remove, nodes_to_add, age_vals, race_vals) {
    dat$attr$age <- c(dat$attr$age[-nodes_to_remove], sample(age_vals, nodes_to_add, TRUE))
    dat$attr$race <- c(dat$attr$race[-nodes_to_remove], sample(race_vals, nodes_to_add, TRUE))
    
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes_to_remove)
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", nodes_to_remove)
    dat$el[[1]] <- add_vertices(dat$el[[1]], nodes_to_add)

    dat      
  }


  dat_el <- list()
  dat_lt <- list()
  dat_time <- list()

  set.seed(0)
  nw <- network.initialize(1000, dir = FALSE)
  nw %v% "age" <- sample(age_vals, 1000, TRUE)
  nw %v% "race" <- sample(race_vals, 1000, TRUE)
  dat <- list(nw = list(nw),
              attr = list(age = nw %v% "age", 
                          race = nw %v% "race"),
              nwparam = list(list(formation = ~edges + nodecov(~age) + nodefactor(~race),
                                  coef.form = coef[-length(coef)],
                                  coef.diss = list(dissolution = ~offset(edges), duration = 1 + exp(3), coef.adj = 3),
                                  constraints = ~.)),
              control = list(track_duration = TRUE))

  dat <- init_tergmLite(dat)
  dat <- updateModelTermInputs(dat)

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
  require(EpiModel)

  set.seed(0)
  nw <- network.initialize(1000, dir = FALSE)
  nw %v% "age" <- sample(age_vals, 1000, TRUE)
  nw %v% "race" <- sample(race_vals, 1000, TRUE)

  ff <- ~edges + nodecov(~age) + nodefactor(~race)
  
  coef <- c(-10, 0.05, 0.5, 0.3, 0.4, 0.9, 3)
  
  x <- list(nw = nw, 
            formation = ff,
            coef.form = coef[-length(coef)],
            coef.diss = list(dissolution = ~edges, duration = 1 + exp(3), coef.adj = 3),
            constraints = ~.)
  
  init_dat <- function(x, param, init, control, s) {
    dat <- list(nw = list(x$nw),
                attr = list(age = x$nw %v% "age",
                            race = x$nw %v% "race",
                            active = rep(TRUE, network.size(x$nw))),
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
    dat$attr$active <- c(dat$attr$active[-nodes_to_remove], rep(TRUE, nodes_to_add))
    
    dat$el[[1]] <- delete_vertices(dat$el[[1]], nodes_to_remove)
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- delete_vertices(dat$p[[1]]$state$nw0 %n% "lasttoggle", nodes_to_remove)
    dat$el[[1]] <- add_vertices(dat$el[[1]], nodes_to_add)

    dat      
  }

  control <- control.net(nsteps = 12, 
                         initialize.FUN = init_dat,
                         nwupdate.FUN = update_dat,
                         prevalence.FUN = NULL,
                         verbose.FUN = NULL,
                         resimulate.network = TRUE, 
                         tergmLite = TRUE,
                         skip.check = TRUE, 
                         track_duration = TRUE,
            #             extract.summary.stats = TRUE, 
                         save.other = c("edgelist", "lasttoggle", "time"),
             #            monitors = list(~mean.age + concurrent),
                         MCMC_control = list(control.simulate.network.tergm()))

  sim <- netsim(x, NULL, NULL, control)
  
  expect_equal(lapply(nw_el, unclass), lapply(sim$edgelist[[1]], unclass), check.attributes = FALSE)
  expect_equal(nw_lt, sim$lasttoggle[[1]], check.attributes = FALSE)
  expect_equal(nw_time, sim$time[[1]])  

})
