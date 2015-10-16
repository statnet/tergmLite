
## functions to add and delete "vertices"

# adding vertices just means changing n in the attributes of the edgelist

# deleting vertices

#' @export
delete_vertices <- function(el, vid) {

  new.el <- el
  if (length(vid) > 0) {
    el.rows.to.del <- which(el[, 1] %in% vid | el[, 2] %in% vid)
    if (length(el.rows.to.del) > 0) {
      new.el <- el[-el.rows.to.del, , drop = FALSE]
    }
    if (nrow(new.el) > 0) {
      elv <- as.vector(new.el)
      # shifted.elv <- vapply(1:length(elv),
      #                       function(x) elv[x] - sum(elv[x] > vid), FUN.VALUE = integer(1))
      shifted.elv <- shiftVec(elv, vid)
      new.el <- matrix(shifted.elv, ncol = 2)
    }
    attributes(new.el)$n <- attributes(el)$n - length(vid)
  }

  return(new.el)
}
