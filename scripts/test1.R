
library(EpiModel)
library(tergmLite)
library(microbenchmark)
library(profr)

nw <- network.initialize(10000, directed = FALSE)
age <- sample(18:50, 10000, TRUE)
nw %v% "age" <- age

est <- netest(nw,
              formation = ~edges + concurrent + absdiff("age"),
              target.stats = c(3500, 1000, 3500*3),
              coef.diss = dissolution_coefs(~offset(edges), 100))

est <- netest(nw,
              formation = ~edges + concurrent,
              target.stats = c(3500, 1000),
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

fn2 <- function() {
  simulate_network(nw,
                   el = NULL,
                   formation = est$formation,
                   dissolution = est$coef.diss$dissolution,
                   coef.form = est$coef.form,
                   coef.diss = est$coef.diss$coef.adj,
                   time.start = 2,
                   time.slices = 1,
                   time.offset = 0,
                   output = "network")
}
t2 <- fn2()

fn3 <- function() {
  simulate_network(nw,
                   el = NULL,
                   formation = est$formation,
                   dissolution = est$coef.diss$dissolution,
                   coef.form = est$coef.form,
                   coef.diss = est$coef.diss$coef.adj,
                   time.start = 2,
                   time.slices = 1,
                   time.offset = 0,
                   output = "edgelist")
}
t2 <- fn3()

el <- as.edgelist(nw)
attributes(el)$vnames <- NULL
fn4 <- function() {
  simulate_network(nw = nw,
                   el = el,
                   formation = est$formation,
                   dissolution = est$coef.diss$dissolution,
                   coef.form = est$coef.form,
                   coef.diss = est$coef.diss$coef.adj,
                   time.start = 2,
                   time.slices = 1,
                   time.offset = 0,
                   output = "edgelist")
}
t2 <- fn4()


el <- as.edgelist(nw)
attributes(el)$vnames <- NULL
p <- ergm_prep(nw, est$formation, est$coef.diss$dissolution, est$coef.form,
               est$coef.diss$coef.adj, est$constraints, control = control.simulate.network())
fn5 <- function() {
  simulate_network(p = p,
                   el = el,
                   coef.form = est$coef.form,
                   coef.diss = est$coef.diss$coef.adj,
                   time.start = 2,
                   output = "edgelist")
}
t2 <- fn5()

fp <- profr(fn4(), interval = 0.005)
fp
ggplot(fp)

res <- microbenchmark(f(), fn5(), times = 50)
res <- microbenchmark(fn4(), times = 50)

summary(res, unit = "s")
summary(res, unit = "relative")

