
library(EpiModel)
library(tergmLite)
library(microbenchmark)
library(profr)

nw <- network.initialize(10000, directed = FALSE)
age <- sample(18:50, 10000, TRUE)
nw %v% "age" <- age

est <- netest(nw,
              formation = ~edges + concurrent,
              target.stats = c(3500, 1000),
              coef.diss = dissolution_coefs(~offset(edges), 1))
names(est)
dx <- netdx(est, nsims = 1000, dynamic = FALSE)
plot(dx)

sim <- simulate(est$fit)
sim

inst.formula <- update.formula(est$formation, sim ~ .)
sim <- simulate(inst.formula, coef = est$coef.form)
sim

res <- microbenchmark(simulate(inst.formula, coef = est$coef.form))
summary(res, unit = "ms")

sim <- simulate_ergm(inst.formula, coef = est$coef.form)

res <- microbenchmark(simulate_ergm(inst.formula, coef = est$coef.form))
summary(res, unit = "ms")

