
library(EpiModel)
library(tergmLite)
library(microbenchmark)
library(profr)

nw <- network.initialize(10000, directed = FALSE)
age <- sample(18:50, 10000, TRUE)
nw %v% "age" <- age

est <- netest(nw,
              formation = ~edges + concurrent + absdiff("age"),
              target.stats = c(3500, 1000, 3*3500),
              coef.diss = dissolution_coefs(~offset(edges), 1))
names(est)
dx <- netdx(est, nsims = 1000, dynamic = FALSE)
plot(dx)

sim <- simulate(est$fit)
sim

p <- ergm_prep(sim, formation = ~edges + concurrent + absdiff("age"),
               coef = est$coef.form, constraints = ~.)

el <- as.edgelist(sim)
attributes(el)$vnames <- NULL

new.el <- simulate_ergm(p, el, est$coef.form)


res <- microbenchmark(simulate_ergm(p, el, est$coef.form), times = 250)
summary(res, unit = "ms")

sim <- simulate_ergm(inst.formula, coef = est$coef.form)

res <- microbenchmark(simulate_ergm(inst.formula, coef = est$coef.form))
summary(res, unit = "ms")

