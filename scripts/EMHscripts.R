
## Age Mixing Scripts

library(EpiModelHIV)
source("scripts/EMHmodules.R")

# Estimation --------------------------------------------------------------

st <- make_nw.hiv(part.dur = 1025,
                  absdiff.offst = 5.38,
                  prop.conc.male = 0.178,
                  prop.conc.feml = 0.028,
                  ages.male = seq(18, 55, 7/365),
                  ages.feml = seq(18, 55, 7/365))

est <- netest(st$nw,
              formation = ~edges + concurrent(by = "male") +
                absdiffby("age", "male", 5.38) + offset(nodematch("male")),
              target.stats = st$stats[c(1:3, 6)],
              coef.diss = st$coef.diss,
              coef.form = -Inf,
              constraints = st$constraints,
              set.control.ergm = control.ergm(MCMLE.maxit = 500, MPLE.type = "penalized"),
              nonconv.error = TRUE)

save(est, file = "scripts/agemix.est.rda")

param <- param.hiv(aids.stage.mult = 2)
init <- init.hiv(i.prev.male = 0.01, i.prev.feml = 0.01)
control <- control.hiv(simno = 1,
                       nsteps = 50,
                       nsims = 1,
                       ncores = 1,
                       resim.int = 1,
                       initialize.FUN = new.initialize.hiv,
                       aging.FUN = aging.hiv,
                       cd4.FUN = cd4.hiv,
                       vl.FUN = vl.hiv,
                       dx.FUN = dx.hiv,
                       tx.FUN = tx.hiv,
                       deaths.FUN = new.deaths.hiv,
                       births.FUN = new.births.hiv,
                       edges_correct.FUN = edges_correct.hiv,
                       resim_nets.FUN = new.simnet.hiv,
                       infection.FUN = infect.hiv,
                       get_prev.FUN = prevalence.hiv,
                       verbose.FUN = verbose.hiv,
                       verbose = TRUE,
                       verbose.int = 1)

sim <- netsim(est, param, init, control)
