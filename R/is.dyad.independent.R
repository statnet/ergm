#  File R/is.dyad.independent.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################

is.dyad.independent<-function(object,...) UseMethod("is.dyad.independent")

is.dyad.independent.NULL <- function(object, ...) TRUE # By convention.

is.dyad.independent.ergm.model <- function(object, ...){
  ! any(sapply(object$terms, function(term) is.null(term$dependence) || term$dependence))
}

is.dyad.independent.formula<-function(object,response=NULL,basis=NULL,...){
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
  
  m<-ergm.getmodel(form, nw, response=response)
  is.dyad.independent(m)
}

is.dyad.independent.conlist <- function(object, object.obs=NULL, ...){
  dind <- TRUE
  for(con in names(object)){
    if(con=="bd" && isTRUE(all.equal(unlist(object[[con]]),FALSE,check.attributes=FALSE))) next
    if(is.null(object[[con]]$free.dyads)) dind <- FALSE
  }
  dind && if(!is.null(object.obs)) is.dyad.independent(object.obs) else TRUE
}

is.dyad.independent.ergm<-function(object,...){
  with(object,
       is.dyad.independent(formula,object$response,network)
       && is.dyad.independent(object$constrained, object$constrained.obs))
}
