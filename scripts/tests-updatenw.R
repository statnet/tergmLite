
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
