#  File R/to_ergm_Cdouble.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Methods to serialize objects into numeric vectors for passing to the C side.
#'
#' These methods return a vector of [`double`]s. For edge lists, this
#' usually takes the form of a \eqn{2 e + 1}- or \eqn{3 e + 1}-vector,
#' containing the number of edges followed a column-major
#' serialization of the edgelist matrix.
#'
#' @param x object to be serialized.
#' @param ... arguments for methods.
#' 
#' @export
to_ergm_Cdouble <- function(x, ...){
  UseMethod("to_ergm_Cdouble")
}
