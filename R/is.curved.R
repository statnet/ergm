#  File R/is.curved.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################

is.curved<-function(object,...) UseMethod("is.curved")

is.curved.NULL <- function(object, ...) FALSE # By convention.

is.curved.ergm.model <- function(object, ...){
  any(object$etamap$canonical==0)
}

is.curved.formula<-function(object,response=NULL,basis=NULL,...){
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
  is.curved(m)
}

is.curved.ergm<-function(object,...){
  any(object$etamap$canonical==0)
}
