
## Age Mixing Scripts

library(EpiModelHIV)
library(tergmLite)
source("scripts/EMHmodules.R")

# Estimation --------------------------------------------------------------

st <- make_nw.hiv(n = 10000,
                  part.dur = 1025,
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

dx <- netdx(est, nsims = 5, nsteps = 300,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
dx <- netdx(est, nsims = 1000, dynamic = FALSE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6,
                                                     MCMC.interval = 1e4))
plot(dx)

# save(est, file = "scripts/agemix.est.rda")

param <- param.hiv(aids.stage.mult = 2)
init <- init.hiv(i.prev.male = 0.1, i.prev.feml = 0.1)
control <- control.hiv(simno = 1,
                       nsteps = 200,
                       nsims = 1,
                       ncores = 1,
                       resim.int = 1,
                       initialize.FUN = new.initialize.hiv,
                       aging.FUN = new.aging.hiv,
                       cd4.FUN = cd4.hiv,
                       vl.FUN = vl.hiv,
                       dx.FUN = dx.hiv,
                       tx.FUN = tx.hiv,
                       deaths.FUN = new.deaths.hiv,
                       births.FUN = new.births.hiv,
                       edges_correct.FUN = edges_correct.hiv,
                       updatepop.FUN = update_population,
                       updatenwp.FUN = update_nwp,
                       resim_nets.FUN = new.simnet.hiv,
                       infection.FUN = new.infect.hiv,
                       get_prev.FUN = prevalence.hiv,
                       verbose.FUN = verbose.hiv,
                       verbose = TRUE,
                       verbose.int = 1)

load("scripts/agemix.est.rda")
sim <- netsim(est, param, init, control)


at <- 1
dat <- new.initialize.hiv(est, param, init, control, s = 1)     # 1.3148

at <- 2
dat <- new.aging.hiv(dat, at)         # 0.0013
dat <- cd4.hiv(dat, at)               # 0.0051
dat <- vl.hiv(dat, at)                # 0.0034
dat <- dx.hiv(dat, at)                # 0.0001
dat <- tx.hiv(dat, at)                # 0.0025
dat <- new.deaths.hiv(dat, at)        # 0.0020
dat <- new.births.hiv(dat, at)        # 0.0016
dat <- edges_correct.hiv(dat, at)     # 0.0000
dat <- update_population(dat, at)     # 0.0129
dat <- update_nwp(dat, at)            # 0.0015
dat <- new.simnet.hiv(dat, at)        # 0.0248
dat <- new.infect.hiv(dat, at)        # 0.0015
dat <- prevalence.hiv(dat, at)        # 0.0144


fp <- profr(new.initialize.hiv(est, param, init, control, s = 1), interval = 0.005)
fp
ggplot(fp)

res <- microbenchmark(new.simnet.hiv(dat, at), times = 100)
summary(res, unit = "s")




