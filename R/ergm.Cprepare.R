#  File R/ergm.Cprepare.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

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

#' @noRd
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

#' Storing last toggle information in a network
#' 
#' An informal extension to \code{\link{network}} objects allowing
#' some limited temporal information to be stored.
#' WARNING: THIS DOCUMENTATION IS PROVIDED AS A COURTESY, AND THE API
#' DESCRIBED IS SUBJECT TO CHANGE WITHOUT NOTICE, DOWN TO COMPLETE
#' REMOVAL. NOT ALL FUNCTIONS THAT COULD SUPPORT IT DO. USE AT YOUR
#' OWN RISK.
#' 
#' While \code{\link[networkDynamic]{networkDynamic}} provides a flexible,
#' consistent method for storing dynamic networks, the \code{C} routines of
#' \code{\link[=ergm-package]{ergm}} and
#' \code{\link[tergm:tergm-package]{tergm}} required a simpler and more
#' lightweight representation.
#' 
#' Though this is an API intended for internal use, some functions,
#' like \code{\link[tergm]{stergm}} (for EGMME),
#' \code{\link[tergm:simulate.stergm]{simulate}}, and
#' \code{\link[=summary.formula]{summary}} can be passed networks with
#' this information using the following \code{\link{network}} (i.e.,
#' \code{\link{\%n\%}}) attributes:
#' \describe{
#'
#' \item{`"time"`}{the
#' time stamp associated with the network}
#'
#' \item{`"lasttoggle"`}{an integer vector analogous to the one returned by [to_ergm_Cdouble.network()] with an attribute: number of elements, a list of tails, a list of heads, and time the edge was last toggled.}
#' }
#'
#' On the C side, it is represented by a hash table.
#' 
#' Again, this API is subject to change without notice.
#'
#' @aliases lasttoggle last.toggle last-toggle
#' @keywords internal
#' @name lasttoggle
NULL
