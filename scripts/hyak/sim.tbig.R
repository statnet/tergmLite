
library(methods)
library(EpiModelHIV)
library(EpiModelHPC)
library(tergmLite)
# source("EMHmodules.R")

param <- param.hiv()
init <- init.hiv(i.prev.male = 0.05, i.prev.feml = 0.05)
control <- control.hiv(simno = 1,
                       nsteps = 52 * 50,
                       nsims = 16,
                       ncores = 16,
                       resim.int = 1,
                       initialize.FUN = initialize.hiv,
                       aging.FUN = aging.hiv,
                       cd4.FUN = cd4.hiv,
                       vl.FUN = vl.hiv,
                       dx.FUN = dx.hiv,
                       tx.FUN = tx.hiv,
                       deaths.FUN = deaths.hiv,
                       births.FUN = births.hiv,
                       edges_correct.FUN = edges_correct.hiv,
                       updatenwp.FUN = NULL,
                       resim_nets.FUN = simnet.hiv,
                       infection.FUN = infect.hiv,
                       get_prev.FUN = prevalence.hiv,
                       verbose.FUN = verbose.hiv,
                       save.network = FALSE,
                       save.other = NULL,
                       verbose = TRUE,
                       verbose.int = 25)

load("agemix.est.lim.rda")
sim <- netsim_par(est, param, init, control)
save(sim, file = "sim.tbig.rda")
