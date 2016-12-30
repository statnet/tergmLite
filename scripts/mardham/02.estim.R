
suppressMessages(library(EpiModelHIV))
rm(list = ls())

load("scripts/mardham/nwstats.rda")

# 1. Main Model -----------------------------------------------------------

# Initialize network
nw.main <- base_nw_msm(st)

# Assign degree
nw.main <- assign_degree(nw.main, deg.type = "pers", nwstats = st)

# Formulas
formation.m <- ~edges +
                nodemix("race", base = 1) +
                nodefactor("deg.pers", base = c(1, 4)) +
                absdiffnodemix("sqrt.age", "race") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.m <- netest(nw.main,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.m,
                coef.diss = st$coef.diss.m,
                constraints = ~bd(maxout = 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e10,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))


# 2. Casual Model ---------------------------------------------------------

# Initialize network
nw.pers <- nw.main

# Assign degree
nw.pers <- assign_degree(nw.pers, deg.type = "main", nwstats = st)

# Formulas
formation.p <- ~edges +
                nodemix("race", base = 1) +
                nodefactor("deg.main", base = c(1, 3)) +
                concurrent(by = "race") +
                absdiffnodemix("sqrt.age", "race") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.p <- netest(nw.pers,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = st$stats.p,
                coef.diss = st$coef.diss.p,
                constraints = ~bd(maxout = 2),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9,
                                                init.method = "zeros",
                                                MCMLE.maxit = 250))


# Fit inst model ----------------------------------------------------------

# Initialize network
nw.inst <- nw.main

# Assign degree
nw.inst <- set.vertex.attribute(nw.inst, "deg.main", nw.pers %v% "deg.main")
nw.inst <- set.vertex.attribute(nw.inst, "deg.pers", nw.main %v% "deg.pers")
table(nw.inst %v% "deg.main", nw.inst %v% "deg.pers")

# Formulas
formation.i <- ~edges +
                nodefactor(c("deg.main", "deg.pers")) +
                nodefactor(c("race", "riskg"), base = c(3, 8)) +
                nodematch("race") +
                absdiffnodemix("sqrt.age", "race") +
                offset(nodematch("role.class", diff = TRUE, keep = 1:2))

# Fit model
fit.i <- netest(nw.inst,
                formation = formation.i,
                target.stats = st$stats.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e9,
                                                MCMLE.maxit = 250))

# Save data
est <- list(fit.m, fit.p, fit.i)
save(est, file = "scripts/mardham/fit.rda")


# Diagnostics -------------------------------------------------------------

# dx <- netdx(fit.i, nsims = 10000, ncores = 1, dynamic = FALSE,
#             nwstats.formula = ~ edges + nodefactor(c("race", "riskg"), base = 0))
# dx
