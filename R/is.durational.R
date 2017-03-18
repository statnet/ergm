#  File R/is.durational.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###############################################################################
# These functions are used to detect whether a ERGM formula/model/etcs are 
# durational dependent or not, based on the (T)ERGM term used.
# To make an (T)ERGM term durational dependent, simply add an "duration" object
# to terms in InitErgmTerm.duration.R
###############################################################################
is.durational<-function(object,...) UseMethod("is.durational")

is.durational.NULL <- function(object, ...) FALSE # By convention.
is.durational.character <- function(object,...) FALSE # for mon="all"
is.durational.ergm.model <- function(object, ...){
	any(object$duration)
}

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
	
	# New formula (no longer use 'object'):
	form <- ergm.update.formula(object, nw ~ ., from.new="nw")
	# work around when durational dependent terms do not has role="target"
#	if(	deparse(substitute(object))=="monitor")
	m<-ergm.getmodel(form, nw, response=response, role=NULL)
	is.durational(m)
}

