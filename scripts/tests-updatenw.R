
# tests to network update functions

n <- 10000
edges <- 3500
el <- t(replicate(edges, sample(n, 2)))

delete_vertices(el, 8)

deaths <- sample(1:n, 100)

library(microbenchmark)
library(network)

nw <- network.initialize(n, directed = FALSE)
nw <- add.edges(nw, tail = el[, 1], head = el[, 2])

res <- microbenchmark(delete_vertices(el, deaths), delete.vertices(nw, deaths))
summary(res, unit = "s")


test_delete_vertices <- function(n, edges, dfrac) {

  require(network)
  el <- t(replicate(edges, sample(n, 2)))
  attributes(el)$n <- n
  nw <- as.network(el, directed = FALSE, matrix.type = "edgelist")
  el <- as.edgelist(nw)

  dth <- sample(1:n, round(n * dfrac))
  new.nw <- nw
  new.nw <- delete.vertices(new.nw, dth)
  new.el <- delete_vertices(el, dth)
  tel <- as.edgelist(new.nw)

  if (!(all(new.el[, 1] == tel[, 1])) |
      !(all(new.el[, 2] == tel[, 2])) |
      !(attributes(new.el)$n == attributes(tel)$n)) {
    # stop("test problem when n = ", n, " and edges = ", edges, " and dth vec = ", dth)
    browser()
  }

}

sample.n <- 100
max.n <- 1e4
grid <- data.frame(n = sample(1:max.n, sample.n))
grid$edges <- round((grid$n/2) * runif(sample.n, min = 0.1, max = 1.5))
grid$dfrac <- runif(sample.n, min = 0, max = 0.5)

for (i in 1:nrow(grid)) {
  test_delete_vertices(grid$n[i], grid$edges[i], grid$dfrac[i])
}
