

#' @title Fast Version of network::add.vertices for Edgelist-formated Network
#'
#' @description This function performs a simple operation of updating the
#'              edgelist attribute \code{n} that tracks the total network
#'              size implicit in an edgelist representation of the network.
#'
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size
#' @param nv A integer equal to the number of nodes to add to the network
#'           size at the given time step
#' @export
#'
add_vertices <- function(el, nv) {
  attributes(el)$n <- attributes(el)$n + nv
  return(el)
}



#' @title Fast Version of network::delete.vertices for Edgelist-formated Network
#'
#' @description Given a current two-column matrix of edges and a vector of IDs
#'              to delete from the matrix, this function first removes any rows
#'              of the edgelist in which the IDs are present and then permutes
#'              downward the index of IDs on the edgelist that were numerically
#'              larger than the IDs deleted.
#'
#' @param el A two-column matrix of current edges (edgelist) with an attribute
#'           variable \code{n} containing the total current network size
#' @param vid A vector of IDs to delete from the edgelist
#'
#' @export
#'
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
