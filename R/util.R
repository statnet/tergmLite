
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
#' be missing.  Attributes \code{class}, \code{dim}, and \code{vnames} (if present)
#' are not copied from \code{el} to the networkLite.
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
#' nw <- network_initialize(100)
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
networkLite <- function(el, attr = NULL) {
  x <- list(el = el, 
            attr = attr, 
            gal = attributes(el)[setdiff(names(attributes(el)), c("class", "dim", "vnames"))])

  storage.mode(x$el) <- "integer" # some "double" edgelists have been coming through...

  # network size attribute is required
  if (is.null(x$gal$n)) {
    stop("networkLite constructor requires network size attribute.")
  }

  # other common attributes default to FALSE
  if (is.null(x$gal$directed)) {
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

  class(x) <- c("networkLite", "network")
  return(x)
}

#' @name networkLitemethods
#' @title networkLite Methods
#'
#' @description S3 methods for networkLite class, for generics defined in
#'              network package.
#'
#' @param x a \code{networkLite} object.
#' @param attrname the name of a vertex attribute in \code{x}.
#' @param ... any additional arguments.
#'
#' @details Allows use of networkLite objects in \code{ergm_model}.
#'
#' @rdname networkLitemethods
#' @export
#'
get.vertex.attribute.networkLite <- function(x,attrname,...) {
  x$attr[[attrname]]
}

#' @rdname networkLitemethods
#' @export
#'
list.vertex.attributes.networkLite <- function(x) {
  sort(unique(names(x$attr)))
}

#' @rdname networkLitemethods
#' @export
#'
network.edgecount.networkLite <- function(x, ...) {
  NROW(x$el)
}

#' @importFrom tibble as_tibble
#' @rdname networkLitemethods
#' @param output Type of edgelist to output.
#' @export
#'
as.edgelist.networkLite <- function(x, output=c("matrix","tibble"), ...) {
  output <- match.arg(output)

  if(output == "matrix") {
    m <- x$el
  } else {
    m <- as_tibble(list(.tail = x$el[,1], .head = x$el[,2]))
    class(m) <- c("tibble_edgelist", "edgelist", class(m))
  }

  attr(m, "n") <- as.integer(network.size(x))
  attr(m, "directed") <- as.logical(is.directed(x))
  bip <- if(is.bipartite(x)) x %n% "bipartite" else FALSE
  attr(m, "bipartite") <- if(is.numeric(bip)) as.integer(bip) else bip
  attr(m, "loops") <- as.logical(has.loops(x))
  attr(m, "vnames") <- network.vertex.names(x)
  
  return(m)
}

#' @rdname networkLitemethods
#' @param i,j Nodal indices (must be missing for networkLite method)
#' @param value Value to set edges to (must be FALSE for networkLite method)
#' @export
"[<-.networkLite" <- function(x, i, j, value) {
  if(missing(i) && missing(j) && isTRUE(all(value == FALSE))) {
    x$el <- structure(x$el[NULL,,drop=FALSE], class = class(x$el), n = x$gal$n)
    return(x)
  } else {
    stop("networkLite `[<-` operator only supports removing all edges at this time")
  }
}

#' @rdname networkLitemethods
#' @export
print.networkLite <- function(x, ...) {
  cat("networkLite with properties:\n")
  cat("  Network size:", x$gal$n, "\n")
  cat("  Edge count:", NROW(x$el), "\n")
  cat("  Network attributes:", sort(unique(names(x$gal))), "\n")
  cat("  Vertex attributes:", sort(unique(names(x$attr))), "\n")
  invisible(x)
}
