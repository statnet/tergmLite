% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{networkLitemethods}
\alias{networkLitemethods}
\alias{get.vertex.attribute.networkLite}
\alias{set.vertex.attribute.networkLite}
\alias{list.vertex.attributes.networkLite}
\alias{get.network.attribute.networkLite}
\alias{set.network.attribute.networkLite}
\alias{list.network.attributes.networkLite}
\alias{network.edgecount.networkLite}
\alias{as.edgelist.networkLite}
\alias{mixingmatrix.networkLite}
\alias{[<-.networkLite}
\alias{print.networkLite}
\alias{network.naedgecount.networkLite}
\alias{add.edges.networkLite}
\alias{as.networkLite}
\alias{as.networkLite.network}
\alias{as.networkLite.networkLite}
\title{networkLite Methods}
\usage{
\method{get.vertex.attribute}{networkLite}(x, attrname, ...)

\method{set.vertex.attribute}{networkLite}(x, attrname, value, v = seq_len(network.size(x)), ...)

\method{list.vertex.attributes}{networkLite}(x, ...)

\method{get.network.attribute}{networkLite}(x, attrname, ...)

\method{set.network.attribute}{networkLite}(x, attrname, value, ...)

\method{list.network.attributes}{networkLite}(x, ...)

\method{network.edgecount}{networkLite}(x, ...)

\method{as.edgelist}{networkLite}(x, output = c("matrix", "tibble"), ...)

\method{mixingmatrix}{networkLite}(object, attr, ...)

\method{[}{networkLite}(x, i, j) <- value

\method{print}{networkLite}(x, ...)

\method{network.naedgecount}{networkLite}(x, ...)

\method{add.edges}{networkLite}(
  x,
  tail,
  head,
  names.eval = NULL,
  vals.eval = NULL,
  ...,
  check.unique = FALSE
)

as.networkLite(x, ...)

\method{as.networkLite}{network}(x, ...)

\method{as.networkLite}{networkLite}(x, ...)
}
\arguments{
\item{x}{a \code{networkLite} object.}

\item{attrname}{the name of an attribute in \code{x}.}

\item{...}{any additional arguments.}

\item{value}{Value to set edges to (must be FALSE for networkLite method)}

\item{v}{indices at which to set vertex attribute values.}

\item{output}{Type of edgelist to output.}

\item{object}{a \code{networkLite} object}

\item{attr}{specification of a vertex attribute in \code{object} as
described in \code{\link[ergm]{nodal_attributes}}}

\item{i, j}{Nodal indices (must be missing for networkLite method)}

\item{tail}{vector of tails of edges to add to the networkLite}

\item{head}{vector of heads of edges to add to the networkLite}

\item{names.eval}{currently unsupported by add.edges.networkLite}

\item{vals.eval}{currently unsupported by add.edges.networkLite}

\item{check.unique}{should a check to ensure uniqueness of edges
in the final edgelist be performed?}
}
\description{
S3 methods for networkLite class, for generics defined in
             network package.
}
\details{
Allows use of networkLite objects in \code{ergm_model}.
}
