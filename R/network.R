
#' @title Initialize Network Object
#'
#' @description Initialize an undirected network object for use in EpiModel
#'              workflows.
#'
#' @param n Network size.
#' @param ... Additional arguments to pass to \code{network.initialize}.
#'
#' @details
#' This function is used in \code{EpiModel} workflows to initialize an empty
#' network object with the directed network attribute hard set to \code{FALSE}.
#'
#' @return
#' Returns an object of class \code{network}.
#'
#' @export
#'
#' @examples
#' nw <- network_initialize(100)
#' nw
#'
network_initialize <- function(n, ...) {
  nw <- network.initialize(n, directed = FALSE, ...)
  return(nw)
}
