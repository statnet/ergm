#  File R/is.durational.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
###############################################################################
# These functions are used to detect whether a ERGM formula/model/etcs are 
# durational dependent or not, based on the (T)ERGM term used.
# To make an (T)ERGM term durational dependent, simply add an "duration" object
# to terms in InitErgmTerm.duration.R
###############################################################################


#' Testing for durational dependent models
#' 
#' These functions test whether an ERGM model or formula is durational
#' dependent or not. If the formula or model does not include any terms that
#' need information about the duration of existing ties, the ergm proceass can
#' use more efficient internal data structures.
#' 
#' @param object An \code{\link{ergm}} object or an ERGM formula, or some
#' characters, e.g., object="all" for monitoring purpose.
#' @param \dots Unused at this time.
#' @return \code{TRUE} if the ERGM terms in the formula or model are durational
#' dependent ; \code{FALSE} otherwise.
#' @keywords model
#' @export
is.durational<-function(object,...) UseMethod("is.durational")

#' @rdname is.durational
#' @description The method for `NULL` always returns `FALSE` by
#'   convention.
#' @export
is.durational.NULL <- function(object, ...) FALSE # By convention.

#' @rdname is.durational
#' @description The method for `character` always returns `FALSE` by
#'   convention.
#' @export
is.durational.character <- function(object,...) FALSE # for mon="all"

#' @describeIn ergm_model Test if the model has duration-dependent terms, which call for [lasttoggle] data structures.
#' @export
is.durational.ergm_model <- function(object, ...){
	any(object$duration)
}

#' @rdname is.durational
#' @template response
#' @param basis See [ergm()].
#' @export
is.durational.formula<-function(object,response=NULL,basis=NULL,...){
	# If basis is not null, replace network in formula by basis.
	# In either case, let nw be network object from formula.
	if(is.null(nw <- basis)) {
		nw <- ergm.getnetwork(object)
	}
	
	nw <- as.network(nw)
	if(!is.network(nw)){
		stop("A network object on the LHS of the formula or via",
				" the 'basis' argument must be given")
	}
	
	# work around when durational dependent terms do not has role="target"
#	if(	deparse(substitute(object))=="monitor")
	m<-ergm_model(object, nw, response=response, role=NULL, ...)
	is.durational(m)
}

