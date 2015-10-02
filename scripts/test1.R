
library(EpiModel)
library(microbenchmark)
library(profr)

nw <- network.initialize(10000, directed = FALSE)
age <- sample(18:50, 10000, TRUE)
nw %v% "age" <- age

est <- netest(nw,
              formation = ~edges + concurrent + absdiff("age"),
              target.stats = c(3500, 1000, 3500*3),
              coef.diss = dissolution_coefs(~offset(edges), 100))

dx <- netdx(est, nsims = 5, nsteps = 500)
dx

nw <- simulate(est$fit)


# Base function -----------------------------------------------------------

f <- function() {
  simulate(nw,
           formation = est$formation,
           dissolution = est$coef.diss$dissolution,
           coef.form = est$coef.form,
           coef.diss = est$coef.diss$coef.adj,
           time.start = 2,
           time.slices = 1,
           time.offset = 0)
}

# Start unpacking ---------------------------------------------------------

fn <- function() {
  simulate_network(nw,
           formation = est$formation,
           dissolution = est$coef.diss$dissolution,
           coef.form = est$coef.form,
           coef.diss = est$coef.diss$coef.adj,
           time.start = 2,
           time.slices = 1,
           time.offset = 0)
}
t2 <- fn()

fp <- profr(fn())
fp
ggplot(fp)

res <- microbenchmark(f(), fn(), times = 25)
summary(res, unit = "relative")

