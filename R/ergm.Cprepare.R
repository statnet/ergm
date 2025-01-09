#  File R/ergm.Cprepare.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' @describeIn to_ergm_Cdouble
#'
#' Method for [`network`] objects.
#'
#' @param attrname name of an edge attribute.
#' 
#' @export
to_ergm_Cdouble.network <- function(x, attrname=NULL, ...){
  xm <- as.edgelist(x, attrname=attrname)
  c(nrow(xm),c(na.omit(xm)))
}

#' @describeIn to_ergm_Cdouble
#'
#' Method for [`ergm_state`] objects, extracting their edgelists.
#'
#' @export
to_ergm_Cdouble.ergm_state <- to_ergm_Cdouble.network


#' @describeIn to_ergm_Cdouble
#'
#' Method for [`matrix`] objects, assumed to be edgelists.
#'
#' @param prototype A network whose relevant attributes (size,
#'   directedness, bipartitedness, and presence of loops) are imposed
#'   on the output edgelist if \code{x} is already an edgelist. (For
#'   example, if the prototype is undirected, `to_ergm_Cdouble`
#'   will ensure that \eqn{t < h}.)
#' @keywords internal
#' @export
to_ergm_Cdouble.matrix <- function(x, prototype=NULL, ...){
  x <- if(!is.null(prototype)) as.edgelist(x, n=network.size(prototype), directed=is.directed(prototype),
                                           bipartite=if(is.bipartite(prototype)) prototype%n%"bipartite" else 0,
                                           loops=has.loops(prototype))
       else x[order(x[,1],x[,2]),,drop=FALSE]
  c(nrow(x),c(na.omit(x)))
}
