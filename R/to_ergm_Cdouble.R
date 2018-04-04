#' Methods to serialize objects into numeric vectors for passing to the C side.
#'
#' @param x object to be serialized.
#' @param ... arguments for methods.
#' @return A [`double`]-precision vector.
#'
#' For serialization of edge lists, this usually takes the form of a
#' \eqn{2 e + 1}-vector, whose first element is the number of edges
#' and whose subsequent elements are tails and heads of the edges,
#' sorted in tail-major order, with an optional third attribute
#' column.
#' 
#' @export
to_ergm_Cdouble <- function(x, ...){
  UseMethod("to_ergm_Cdouble")
}
