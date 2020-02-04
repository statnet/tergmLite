
## Full scale EpiModelHIV example

# package install
install.packages("EpiModel", dep = TRUE)
remotes::install_github(c("statnet/EpiModel@2c131f0",
                          "statnet/EpiModelHPC@a64dbf2",
                          "statnet/tergmLite@v2.1.6",
                          "EpiModel/EpiABC@c32ecb6",
                          "EpiModel/ARTnetData@1d8ec6e",
                          "EpiModel/ARTnet@150c631"),
                        upgrade = FALSE)
remotes::install_github("EpiModel/EpiModelHIV-p", ref = "CombPrev", upgrade = FALSE)

## Packages ##
rm(list = ls())
suppressMessages(library("EpiModelHIV"))
suppressMessages(library("ARTnet"))

epistats <- build_epistats(city_name = "Atlanta")
netparams <- build_netparams(epistats = epistats, smooth.main.dur.55p = TRUE)
netstats <- build_netstats(epistats, netparams,
                           expect.mort = 0.000478213, network.size = 10000)

save(epistats, netparams, netstats, file = "inst/ARTnetTestParams.rda")

# 0. Initialize Network ---------------------------------------------------

num <- netstats$demog$num
nw <- network::network.initialize(num, directed = FALSE)

attr.names <- names(netstats$attr)
attr.values <- netstats$attr
nw <- network::set.vertex.attribute(nw, attr.names, attr.values)
nw_main <- nw_casl <- nw_inst <- nw


# 1. Main Model -----------------------------------------------------------

# Formula
model_main <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.casl", levels = -1) +
  concurrent +
  degrange(from = 3) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_main <- c(
  edges = netstats$main$edges,
  nodematch_age.grp = netstats$main$nodematch_age.grp,
  nodefactor_age.grp = netstats$main$nodefactor_age.grp[-1],
  nodematch_race = netstats$main$nodematch_race_diffF,
  nodefactor_race = netstats$main$nodefactor_race[-1],
  nodefactor_deg.casl = netstats$main$nodefactor_deg.casl[-1],
  concurrent = netstats$main$concurrent,
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
cbind(netstats_main)
netstats_main <- unname(netstats_main)

# Fit model
fit_main <- netest(nw_main,
                   formation = model_main,
                   target.stats = netstats_main,
                   coef.diss = netstats$main$diss.byage,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 3,
                                                   SAN.nsteps.times = 3),
                   verbose = FALSE)



# 2. Casual Model ---------------------------------------------------------

# Formula
model_casl <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = c(-1, -5)) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("deg.main", levels = -3) +
  concurrent +
  degrange(from = 4) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_casl <- c(
  edges = netstats$casl$edges,
  nodematch_age.grp = netstats$casl$nodematch_age.grp,
  nodefactor_age.grp = netstats$casl$nodefactor_age.grp[-c(1,5)],
  nodematch_race = netstats$casl$nodematch_race_diffF,
  nodefactor_race = netstats$casl$nodefactor_race[-1],
  nodefactor_deg.main = netstats$casl$nodefactor_deg.main[-3],
  concurrent = netstats$casl$concurrent,
  degrange = 0,
  nodematch_role.class = c(0, 0)
)
cbind(netstats_casl)
netstats_casl <- unname(netstats_casl)

# Fit model
fit_casl <- netest(nw_casl,
                   formation = model_casl,
                   target.stats = netstats_casl,
                   coef.diss = netstats$casl$diss.byage,
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 3,
                                                   SAN.nsteps.times = 3),
                   verbose = FALSE)


# 3. One-Off Model --------------------------------------------------------

# Formula
model_inst <- ~edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = -1) +
  nodematch("race", diff = FALSE) +
  nodefactor("race", levels = -1) +
  nodefactor("risk.grp", levels = -5) +
  nodefactor("deg.tot", levels = -1) +
  nodematch("role.class", diff = TRUE, levels = 1:2)

# Target Stats
netstats_inst <- c(
  edges = netstats$inst$edges,
  nodematch_age.grp = sum(netstats$inst$nodematch_age.grp),
  nodefactor_age.grp = netstats$inst$nodefactor_age.grp[-1],
  nodematch_race = netstats$inst$nodematch_race_diffF,
  nodefactor_race = netstats$inst$nodefactor_race[-1],
  nodefactor_risk.grp = netstats$inst$nodefactor_risk.grp[-5],
  nodefactor_deg.tot = netstats$inst$nodefactor_deg.tot[-1],
  nodematch_role.class = c(0, 0)
)
cbind(netstats_inst)
netstats_inst <- unname(netstats_inst)

# Fit model
fit_inst <- netest(nw_inst,
                   formation = model_inst,
                   target.stats = netstats_inst,
                   coef.diss = dissolution_coefs(~offset(edges), 1),
                   set.control.ergm = control.ergm(MCMLE.maxit = 500,
                                                   SAN.maxit = 3,
                                                   SAN.nsteps.times = 3),
                   verbose = FALSE)

out <- list(fit_main, fit_casl, fit_inst)


# 4. Run Diagnostics ------------------------------------------------------

## Main

fit_main <- out[[1]]

model_main_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = NULL) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = NULL) +
  nodefactor("deg.casl", levels = NULL) +
  degrange(from = 3) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:3)
dx_main <- netdx(fit_main, nsims = 10, ncores = 6, nsteps = 500,
                 nwstats.formula = model_main_dx, skip.dissolution = TRUE,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx_main)
plot(dx_main)

netstats$main


## Casual

fit_casl <- out[[2]]

model_casl_dx <- ~edges +
  nodematch("age.grp", diff = TRUE) +
  nodefactor("age.grp", levels = NULL) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = NULL) +
  nodefactor("deg.main", levels = NULL) +
  degrange(from = 4) +
  concurrent +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)
dx_casl <- netdx(fit_casl, nsims = 10, ncores = 6, nsteps = 500,
                 nwstats.formula = model_casl_dx, skip.dissolution = TRUE,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))
print(dx_casl, digits = 1)
plot(dx_casl)

netstats$casl


## One-Off

fit_inst <- out[[3]]

model_inst_dx <- ~edges +
  nodematch("age.grp", diff = FALSE) +
  nodefactor("age.grp", levels = NULL) +
  nodematch("race", diff = TRUE) +
  nodefactor("race", levels = NULL) +
  nodefactor("risk.grp", levels = NULL) +
  nodefactor("deg.tot", levels = NULL) +
  nodematch("role.class", diff = TRUE) +
  degree(0:4)
dx_inst <- netdx(fit_inst, nsims = 10000, dynamic = FALSE,
                 nwstats.formula = model_inst_dx,
                 set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e5))

print(dx_inst, digits = 2)

plot(dx_inst, sim.lines = TRUE, sim.lwd = 0.05)

netstats$inst


# 5. Run Epidemic Model ---------------------------------------------------


param <- param_msm(netstats = netstats,
                   epistats = epistats,
                   hiv.test.rate = c(0.00432, 0.00425, 0.00730),
                   hiv.test.late.prob = c(0, 0, 0),
                   tx.init.prob = c(0.1775, 0.190, 0.2521),
                   tt.part.supp = c(0.45, 0.40, 0.28),
                   tt.full.supp = c(0.55, 0.60, 0.72),
                   tt.dur.supp = c(0, 0, 0),
                   tx.halt.part.prob = c(0.009, 0.0084, 0.00768),
                   tx.halt.full.rr = c(0.45, 0.45, 0.45),
                   tx.halt.dur.rr = c(0.45, 0.45, 0.45),
                   tx.reinit.part.prob = c(0.0115, 0.0135, 0.0205),
                   tx.reinit.full.rr = c(1, 1, 1),
                   tx.reinit.dur.rr = c(1, 1, 1),
                   max.time.off.tx.full.int = 52 * 15,
                   max.time.on.tx.part.int = 52 * 10,
                   max.time.off.tx.part.int = 52 * 10,
                   aids.mr = 1/250,
                   trans.scale = c(2.64, 0.45, 0.285),
                   acts.scale = 1.00,
                   acts.aids.vl = 5.75,
                   prep.start = (52*60) + 1,
                   riskh.start = 52*59,
                   prep.start.prob = 0.66,
                   prep.require.lnt = TRUE,
                   prep.risk.reassess.method = "year")
init <- init_msm(prev.ugc = 0,
                 prev.rct = 0,
                 prev.rgc = 0,
                 prev.uct = 0)
control <- control_msm(simno = 1,
                       nsteps = 52*5,
                       nsims = 1,
                       ncores = 1,
                       save.nwstats = TRUE,
                       save.clin.hist = FALSE)

sim <- netsim(out, param, init, control)

df <- as.data.frame(sim, out = "mean")
names(df)

df$cc.test.int

par(mar = c(3,3,1,1), mgp = c(2,1,0))
plot(sim, y = "i.prev", mean.smooth = FALSE, ylim = c(0, 1))
