
#' @title networkLite Constructor Utility
#'
#' @description Constructor function for a networkLite object.
#'
#' @param el an edgelist-formatted network representation, including network
#'        attributes.
#' @param attr a list of named vertex attributes for the network represented
#'        by \code{el}.
#'
#' @details
#' This function takes an edge list \code{el} with network attributes attached,
#' and a list of vertex attributes \code{attr}, and returns a networkLite object,
#' which is a list with named fields \code{el}, \code{attr}, and \code{gal}, with
#' each of the first two corresponding to the argument of the same name, and
#' \code{gal} being the list of network attributes (copied from \code{attributes(el)})
#' for compatibility with some \code{network} accessors. Missing attributes
#' \code{directed}, \code{bipartite}, \code{loops}, \code{hyper}, and \code{multiple}
#' are defaulted to \code{FALSE}. The network size attribute \code{n} must not
#' be missing.
#'
#' This new data structure is then used within the \code{\link{updateModelTermInputs}}
#' function for updating the structural information on the network used for ERGM
#' simulation.
#'
#' @return
#' A networkLite object with edge list \code{el}, vertex attributes \code{attr},
#' and network attributes \code{gal}.
#'
#' @rdname networkLite
#' @export
#'
#' @examples
#' library("EpiModel")
#' nw <- network.initialize(n = 100, directed = FALSE)
#' formation <- ~edges
#' target.stats <- 50
#' coef.diss <- dissolution_coefs(dissolution = ~offset(edges), duration = 20)
#' x <- netest(nw, formation, target.stats, coef.diss, verbose = FALSE)
#'
#' param <- param.net(inf.prob = 0.3)
#' init <- init.net(i.num = 10)
#' control <- control.net(type = "SI", nsteps = 100, nsims = 5, depend = TRUE)
#'
#' # Full network structure after initialization
#' dat <- crosscheck.net(x, param, init, control)
#' dat <- initialize.net(x, param, init, control)
#'
#' # networkLite representation used by tergmLite
#' dat <- init_tergmLite(dat)
#'
#' # Conversion to networkLite class format
#' nwl <- networkLite(dat$el[[1]], dat$attr)
#' nwl
#'
networkLite <- function(el, attr) {
  x <- list(el = el, attr = attr, gal = attributes(el))

  # network size attribute is required
  if (is.null(x$gal$n)) {
    stop("networkLite constructor requires network size attribute.")
  }

  # other common attributes default to FALSE
  if (is.null(x$gal$directed))  {
    x$gal$directed <- FALSE
  }
  if (is.null(x$gal$bipartite)) {
    x$gal$bipartite <- FALSE
  }
  if (is.null(x$gal$loops)) {
    x$gal$loops <- FALSE
  }
  if (is.null(x$gal$hyper)) {
    x$gal$hyper <- FALSE
  }
  if (is.null(x$gal$multiple)) {
    x$gal$multiple <- FALSE
  }

  class(x) <- "networkLite"
  return(x)
}

#' @name networkLitemethods
#' @title networkLite Methods
#'
#' @description S3 methods for networkLite class, for generics defined in
#'              network package.
#'
#' @param nw a \code{networkLite} object.
#' @param x a \code{networkLite} object.
#' @param attrname the name of a vertex attribute in \code{x}.
#' @param ... any additional arguments.
#'
#' @details Allows use of networkLite objects in \code{ergm_model}.
#'
#' @importFrom network as.network
#'
#' @rdname networkLitemethods
#' @export
#'
as.network.networkLite <- function(nw, ...) {
  class(nw) <- c("networkLite", "network")
  nw
}

#' @importFrom network get.vertex.attribute
#' @rdname networkLitemethods
#' @export
#'
get.vertex.attribute.networkLite <- function(x,attrname,...) {
  x$attr[[attrname]]
}

#' @importFrom network list.vertex.attributes
#' @rdname networkLitemethods
#' @export
#'
list.vertex.attributes.networkLite <- function(x) {
  sort(unique(names(x$attr)))
}
