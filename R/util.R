
deinf <- function (x, replace = 1/.Machine$double.eps) {
  ifelse(is.nan(x) | abs(x) < replace, x, sign(x) * replace)
}

#' @rdname networkLite
#' @title networkLite constructor
#' @description Constructor function for networkLite object.
#' @details This function takes an edge list \code{el} with network attributes attached, and
#'   a list of vertex attributes \code{attr}, and returns a networkLite object, which is a list
#'   with named fields \code{el}, \code{attr}, and \code{gal}, with each of the first two corresponding to
#'   the argument of the same name, and \code{gal} being the list of network attributes
#'   (copied from \code{attributes(el)}) for compatibility with some \code{network} accessors.
#'   Missing attributes \code{directed}, \code{bipartite}, \code{loops}, \code{hyper}, and \code{multiple}
#'   are defaulted to \code{FALSE}.  The network size attribute \code{n} must not be missing.
#' @param el an edgelist-formatted network representation, including network attributes
#' @param attr a list of named vertex attributes for the network represented by \code{el}
#' @return a networkLite object with edge list \code{el}, vertex attributes \code{attr}, and
#'   network attributes \code{gal}
#' @export
networkLite <- function(el, attr) {
  x <- list(el = el, attr = attr, gal = attributes(el))

  # network size attribute is required
  if(is.null(x$gal$n)) stop("networkLite constructor requires network size attribute.")

  # other common attributes default to FALSE
  if(is.null(x$gal$directed))  x$gal$directed  <- FALSE
  if(is.null(x$gal$bipartite)) x$gal$bipartite <- FALSE
  if(is.null(x$gal$loops))     x$gal$loops     <- FALSE
  if(is.null(x$gal$hyper))     x$gal$hyper     <- FALSE
  if(is.null(x$gal$multiple))  x$gal$multiple  <- FALSE

  class(x) <- "networkLite"
  x
}

#' @name networkLitemethods
#' @title networkLite methods
#' @description S3 methods for networkLite class, for generics defined in network package.
#' @details Allows use of networkLite objects in \code{ergm_model}.
#' @param nw a networkLite object
#' @param x a networkLite object
#' @param attrname the name of a vertex attribute in \code{x}
#' @param ... any additional arguments

#' @importFrom network as.network
#' @rdname networkLitemethods
#' @export
as.network.networkLite <- function(nw, ...) {
  class(nw) <- c("networkLite", "network")
  nw
}

#' @importFrom network get.vertex.attribute
#' @rdname networkLitemethods
#' @export
get.vertex.attribute.networkLite <- function(x,attrname,...) {
  x$attr[[attrname]]
}

#' @importFrom network list.vertex.attributes
#' @rdname networkLitemethods
#' @export
list.vertex.attributes.networkLite <- function(x) {
  sort(unique(names(x$attr)))
}
