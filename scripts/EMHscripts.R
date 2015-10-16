
## Age Mixing Scripts

library(EpiModelHIV)
library(tergmLite)
library(microbenchmark)
library(profr)
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
                       nsteps = 25,
                       nsims = 2,
                       ncores = 2,
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
                       updatenwp.FUN = update_nwp,
                       resim_nets.FUN = new.simnet.hiv,
                       infection.FUN = new.infect.hiv,
                       get_prev.FUN = prevalence.hiv,
                       verbose.FUN = verbose.hiv,
                       save.nwstats = FALSE,
                       verbose = TRUE,
                       verbose.int = 1)

load("scripts/agemix.est.rda")
sim <- netsim_par(est, param, init, control)


at <- 1
dat <- new.initialize.hiv(est, param, init, control, s = 1)

tlf <- function(dat, at = 2) {
  dat <- new.aging.hiv(dat, at) # 0.0014
  dat <- cd4.hiv(dat, at) # 0.0014
  dat <- vl.hiv(dat, at) # 0.0029
  dat <- dx.hiv(dat, at) # 0.0006
  dat <- tx.hiv(dat, at) # 0.0025
  dat <- new.deaths.hiv(dat, at) # 0.0048
  dat <- new.births.hiv(dat, at) # 0.0016
  dat <- edges_correct.hiv(dat, at) # Nada
  dat <- update_nwp(dat, at) # 0.0009
  dat <- new.simnet.hiv(dat, at) # 0.0191
  dat <- new.infect.hiv(dat, at) # 0.0010
  dat <- prevalence.hiv(dat, at) # 0.0109
  return(dat)
}
# dat <- tlf(dat, at = 6)
# dat$el

# fp <- profr(tlf(dat, at = 2), interval = 0.005)
# fp
# ggplot(fp)

res <- microbenchmark(tlf(dat, at = 2), times = 100)
summary(res, unit = "s")

library(ggplot2)
autoplot(res)



# apples to apples benchmark ----------------------------------------------

control.old <- control.hiv(simno = 1,
                       nsteps = 1000,
                       nsims = 1,
                       ncores = 1,
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
                       verbose = TRUE,
                       verbose.int = 1)

at <- 1
dat2 <- initialize.hiv(est, param, init, control.old, s = 1)

tlf.old <- function(dat2, at = 2) {
  dat2 <- aging.hiv(dat2, at)
  dat2 <- cd4.hiv(dat2, at)
  dat2 <- vl.hiv(dat2, at)
  dat2 <- dx.hiv(dat2, at)
  dat2 <- tx.hiv(dat2, at)
  dat2 <- deaths.hiv(dat2, at)
  dat2 <- births.hiv(dat2, at)
  dat2 <- edges_correct.hiv(dat2, at)
  dat2 <- simnet.hiv(dat2, at)
  dat2 <- infect.hiv(dat2, at)
  dat2 <- prevalence.hiv(dat2, at)
  return(dat2)
}

res <- microbenchmark(tlf.old(dat2, at = 2), tlf(dat, at = 2), times = 100)
res <- microbenchmark(tlf(dat, at = 2))
summary(res, unit = "s")
summary(res, unit = "relative")

library(ggplot2)

pdf(file = "scripts/hyak/timing.pdf", height = 4, width = 8)
autoplot(res, log = FALSE)
dev.off()
autoplot(res, log = TRUE)
boxplot(res, outline = FALSE)

# validation of old versus new --------------------------------------------

load("scripts/hyak/sim.tbig.rda")
big <- sim

load("scripts/hyak/sim.tlite.rda")
lite <- sim

pdf(file = "scripts/hyak/validate.pdf", height = 4, width = 8)
par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(big, y = "i.num", ylim = c(0, 0.25), qnts = 1)
plot(lite, y = "i.num", qnts = 1, add = TRUE, mean.col = "firebrick",
     mean.lty = 2, qnts.col = "firebrick")
dev.off()

