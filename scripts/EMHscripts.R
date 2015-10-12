
## Age Mixing Scripts

library(EpiModelHIV)
library(tergmLite)
source("scripts/EMHmodules.R")

# Estimation --------------------------------------------------------------

st <- make_nw.hiv(n = 10000,
                  part.dur = 2000,
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

est <- netest(st$nw,
              formation = ~edges + offset(nodematch("male")),
              target.stats = st$stats[1],
              coef.diss = st$coef.diss,
              coef.form = -Inf,
              constraints = st$constraints,
              set.control.ergm = control.ergm(MCMLE.maxit = 500, MPLE.type = "penalized"),
              nonconv.error = TRUE)


dx <- netdx(est, nsims = 5, nsteps = 500,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
dx <- netdx(est, nsims = 1000, dynamic = FALSE,
            set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6,
                                                     MCMC.interval = 1e4))
plot(dx, plots.joined = FALSE)

# save(est, file = "scripts/agemix.est.rda")


# Simulation --------------------------------------------------------------

param <- param.hiv()
init <- init.hiv(i.prev.male = 0.05, i.prev.feml = 0.05)
control <- control.hiv(simno = 1,
                       nsteps = 1000,
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
                       updatenwp.FUN = update_nwp_lim,
                       resim_nets.FUN = new.simnet.hiv,
                       infection.FUN = new.infect.hiv,
                       get_prev.FUN = prevalence.hiv,
                       verbose.FUN = verbose.hiv,
                       verbose = TRUE,
                       verbose.int = 1)

# load("scripts/agemix.est.rda")
# load("scripts/agemix.est.lim.rda")
sim <- netsim(est, param, init, control)


at <- 1
dat <- new.initialize.hiv(est, param, init, control, s = 1)

tlf <- function(dat, at = 2) {
  dat <- new.aging.hiv(dat, at)
  dat <- cd4.hiv(dat, at)
  dat <- vl.hiv(dat, at)
  dat <- dx.hiv(dat, at)
  dat <- tx.hiv(dat, at)
  dat <- new.deaths.hiv(dat, at)
  dat <- new.births.hiv(dat, at)
  dat <- edges_correct.hiv(dat, at)
  dat <- update_nwp_lim(dat, at)
  dat <- new.simnet.hiv(dat, at)
  dat <- new.infect.hiv(dat, at)
  dat <- prevalence.hiv(dat, at)
  return(dat)
}
dat <- tlf(dat, at = 6)
dat$el

fp <- profr(new.initialize.hiv(est, param, init, control, s = 1), interval = 0.005)
fp
ggplot(fp)

res <- microbenchmark(tlf(dat, at = 2), times = 100)
summary(res, unit = "s")
par(mar = c(3,3,1,1), mgp = c(2,1,0))
boxplot(res, unit = "s", outline = TRUE, log = FALSE)
autoplot(res)
