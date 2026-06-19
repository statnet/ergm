#' @name dyad_attr
#'
#' @title Get and set dyad attribute matrices robust to transformations
#'
#' @param x an object capable of having vertex and network attributes.
#' @param attrname name of an attribute; uses the "namespace" of network attributes.
#' @param value matrix of appropriate dimension.
#' @param ... additional arguments to methods.
NULL

#' @describeIn dyad_attr Set a dyad attribute.
#' @export
set.dyad.attribute <- function(x, attrname, value, ...) UseMethod("set.dyad.attribute")

#' @rdname dyad_attr
#' @export
set.dyad.attribute.network <- function(x, attrname, value, ...) {
  i <- seq_len(network.size(x))
  if (is.bipartite(x)) {
    if ((nrow(value) != b1.size(x) || ncol(value) != b2.size(x)))
      stop("dyad attribute is not a mode-1-size by mode-2-size rectangular matrix")
    i[i > (b <- x%v%"bipartite")] <- i - b
  } else {
    if (any(dim(value) != network.size(x)))
      stop("dyad attribute is not an n by n matrix")
  }

  x %v% paste0(".vid_orig.", attrname) <- i
  x %n% attrname <- value
  modify_in_place(x)
}

#' @describeIn dyad_attr Alias for `set.dyad.attribute(x, attrname, value)`.
#' @export
`%d%<-` <- function(x, attrname, value) UseMethod("set.dyad.attribute")

#' @describeIn dyad_attr Get a dyad attribute, adjusting for any transformations such as taking of a subgraph.
#' @export
get.dyad.attribute <- function(x, attrname) UseMethod("get.dyad.attribute")

#' @rdname dyad_attr
#' @export
get.dyad.attribute.network <- function(x, attrname) {
  a <- x %n% attrname
  va <- paste0(".vid_orig.", attrname)
  if (va %in% list.vertex.attributes(x)) {
    xi <- x %v% va
    # Indexing will automatically pad any vertices that don't have an
    # index with rows or columns of NA.
    if (is.bipartite(x)) {
      xi1 <- xi[seq_len(x %n% "bipartite")]
      xi2 <- xi[-seq_len(x %n% "bipartite")]
      if (!identical(xi1, seq_len(nrow(a))) ||
          !identical(xi2, seq_len(ncol(a)))) a <- a[xi1, xi2]
    } else {
      if (!identical(xi, seq_len(nrow(a)))) a <- a[xi, xi]
    }
  }
  a
}

#' @describeIn dyad_attr Alias for `get.dyad.attribute(x, attrname)`.
#' @export
`%d%` <- function(x, attrname) UseMethod("get.dyad.attribute")
