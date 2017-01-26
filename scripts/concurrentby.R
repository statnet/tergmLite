
# testing for concurrent(by = "attr")

# for stergm_prep
# nw, formation, dissolution,
# coef.form, coef.diss, constraints,
# control = tergm::control.simulate.network()

library("EpiModel")
library("tergmLite")

nw <- network.initialize(100, directed = FALSE)
nw <- set.vertex.attribute(nw, "riskg", rbinom(100, 1, 0.5))
nw %v% "riskg"

est <- netest(nw = nw,
              formation = ~edges + concurrent(by = "riskg"),
              target.stats = c(50, 25, 15),
              coef.diss = dissolution_coefs(~offset(edges), duration = 100))

p <- stergm_prep(nw,
                 est$formation,
                 est$coef.diss$dissolution,
                 est$coef.form,
                 est$coef.diss$coef.adj,
                 est$constraints)

# str(p)
# str(p$model.form$terms)
# head(p$model.form$terms[[2]]$inputs, 20)
# head(nw %v% "riskg", 20)

dat <- list()
dat$p <- dat$el <- dat$attr <- dat$nwparam <- list()
dat$p[[1]] <- p
dat$el[[1]] <- as.edgelist(simulate(est$fit))
dat$attr$riskg <- nw %v% "riskg"
dat$nwparam[[1]] <- est

# debug(updateModelTermInputs)

dat <- updateModelTermInputs(dat, network = 1)

identical(dat$p[[1]], p)




# Test change ---------------------------------------------------------

deaths <- sample(1:100, 10)
dat$el[[1]] <- delete_vertices(dat$el[[1]], vid = deaths)
dat$attr <- deleteAttr(dat$attr, deaths)

sapply(dat$attr, length)

nw2 <- as.network(dat$el[[1]], directed = FALSE)
nw2 <- set.vertex.attribute(nw2, "riskg", dat$attr$riskg)

p2 <- stergm_prep(nw2,
                  est$formation,
                  est$coef.diss$dissolution,
                  est$coef.form,
                  est$coef.diss$coef.adj,
                  est$constraints)

dat <- updateModelTermInputs(dat, network = 1)

dat$p[[1]][[1]]$formula <- NULL
dat$p[[1]][[2]]$formula <- NULL
p2[[1]]$formula <- NULL
p2[[2]]$formula <- NULL

identical(dat$p[[1]], p2)
