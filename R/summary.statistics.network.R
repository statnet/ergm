#  File R/summary.statistics.network.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
#==========================================================================
# This file contains the following 5 functions for computing summary stats
#      <summary.statistics>           <summary.statistics.formula>
#      <summary.formula>              <summary.statistics.ergm>
#      <summary.statisitcs.default>   <summary.statistics.network>
#      <summary.statisitics.matrix>
#==========================================================================




#############################################################################
# Each of the <summary.statistics.X> functions and <summary.formula> checks
# that the implicit formula is correctly given as 'nw ~ term(s)' and returns
# the global statistics of the network specified by the formula
#
# --PARAMETERS--
#   object:  a formula, matrix, ergm, or network, as appropriate
#   basis :  optionally, the network from the formula; if a network
#            is passed to 'basis', it is assumed that 'object' is the
#            formula
#
# --RETURNED--
#   gs: the vector of global stats, as returned by <ergm.getglobalstats>
#############################################################################

summary.statistics <- function(object, ..., basis=NULL) {
  UseMethod("summary.statistics")
}


summary.formula <- function(object, ...){
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

summary.statistics.network.list <- function(object, response=NULL, ..., basis=NULL){
  if(!is.null(basis)){
    if(inherits(basis,'network.list'))
      object[[2]] <- basis
    else stop('basis, if specified, should be the same type as the LHS of the formula (network.list, in this case).')
  }
  nwl <- eval(object[[2]], envir=environment(object))
  out<-lapply(nwl, function(nw) summary.statistics.network(object, response=response, ..., basis=nw))
  do.call(rbind,out)
}

summary.statistics.default <-
summary.statistics.matrix <- 
summary.statistics.network <- function(object, response=NULL,...,basis=NULL) {
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(object)
    formula[[2]] <- as.name("basis") # This seems irrelevant; network name
                                     # not needed by ergm.getmodel
  }else{
    formula <- object
    nw <- ergm.getnetwork(formula)
  }
  m <- ergm.getmodel(formula, nw, response=response, role="target",...)
  gs <- ergm.getglobalstats(nw, m, response=response)
  gs
}




