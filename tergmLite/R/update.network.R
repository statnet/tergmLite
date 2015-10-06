
## functions to add and delete "vertices"

# adding vertices just means changing n in the attributes of the edgelist

# deleting vertices

#' @export
delete_vertices <- function(el, vid) {

  el.rows.to.del <- which(el[, 1] %in% vid | el[, 2] %in% vid)
  new.el <- el[-el.rows.to.del, ]

  elv <- as.vector(new.el)
  shifted.elv <- vapply(1:length(elv),
                        function(x) elv[x] - sum(elv[x] > vid), FUN.VALUE = integer(1))

  # vid <- sort(vid)
  # for (i in 1:length(vid)) {
  #   new.el[new.el > vid[i]] <- new.el[new.el > vid[i]] - 1
  # }

  return(matrix(shifted.elv, ncol = 2))
  # return(new.el)
}
