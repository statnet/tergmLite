
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


test_delete_vertices <- function(n, edges) {

  require(network)
  el <- t(replicate(edges, sample(n, 2)))
  attributes(el)$n <- n
  nw <- as.network(el, directed = FALSE, matrix.type = "edgelist")
  el <- as.edgelist(nw)

  dth <- sample(1:n, round(n * 0.05))
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

max.n <- 1e4
grid <- data.frame(n = seq(10, max.n, 20))
grid <- data.frame(n = sample(1:max.n, 100))
grid$edges <- round((grid$n/2) * runif(nrow(grid), min = 0.1, max = 1.5))

for (i in 1:nrow(grid)) {
  test_delete_vertices(grid$n[i], grid$edges[i])
}
test_delete_vertices(30, 13)


el <- structure(c(1L, 2L, 2L, 3L, 5L, 10L, 10L, 13L, 17L, 8L, 14L,
            16L, 21L, 10L, 18L, 19L, 14L, 23L), .Dim = c(9L, 2L), n = 30, vnames = 1:30, directed = FALSE, bipartite = FALSE, loops = FALSE, inverted = FALSE, class = c("edgelist", "matrix"))
dth <- c(4L, 25L)

new.el <- delete_vertices(el, dth)
