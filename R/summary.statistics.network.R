#  File ergm/R/summary.statistics.network.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#############################################################################
# Each of the <summary.statistics.X> functions and <summary.formula> checks
# that the implicit formula is correctly given as 'nw ~ term(s)' and returns
# the global statistics of the network specified by the formula
#############################################################################

summary.statistics <- function(object, ..., basis=NULL) {
  UseMethod("summary.statistics")
}


summary.formula <- function(object, ...){
  current.warn <- options()$warn
  on.exit(options(warn=current.warn))
  options(warn=0)
  if(length(object)!=3 || object[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  lhs <- eval(object[[2]], envir = environment(object))
  UseMethod("summary.statistics",object=lhs)
}



summary.statistics.formula <- function(object, ..., basis=NULL) {
  summary.statistics.network(object, ..., basis=basis)
}



summary.statistics.ergm <- function(object, ..., basis=NULL)
{
  summary.statistics.network(object$formula, ..., basis=basis)
}

summary.statistics.network.list <- function(object, ..., basis=NULL){
  if(!is.null(basis)){
    if(inherits(basis,'network.list'))
      object[[2]] <- basis
    else stop('basis, if specified, should be the same type as the LHS of the formula (network.list, in this case).')
  }
  nwl <- eval(object[[2]], envir=environment(object))
  out<-lapply(nwl, function(nw) summary.statistics.network(object, ..., basis=nw))
  do.call(rbind,out)
}

summary.statistics.default <-
summary.statistics.matrix <- 
summary.statistics.network <- function(object,...,basis=NULL) {
  current.warn <- options()$warn
  on.exit(options(warn=current.warn))
  options(warn=0)
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(object)
    formula[[2]] <- as.name("basis") # This seems irrelevant; network name
                                     # not needed by ergm.getmodel
  }else{
    formula <- object
    nw <- ergm.getnetwork(formula)
  }
  m <- ergm.getmodel(formula, nw, role="target",...)
  gs <- ergm.getglobalstats(nw, m)
  gs
}




