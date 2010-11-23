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
#   drop  :  whether to drop degenerate terms when the model is
#            constructed from the formula (T or F); default=FALSE
#   basis :  optionally, the network from the formula; if a network
#            is passed to 'basis', it is assumed that 'object' is the
#            formula
#
# --RETURNED--
#   gs: the vector of global stats, as returned by <ergm.getglobalstats>
#############################################################################

summary.statistics <- function(object, ..., drop=FALSE, basis=NULL) {
  UseMethod("summary.statistics")
}


summary.formula <- function(object, ...){
  current.warn <- options()$warn
  options(warn=0)
  trms<-terms(object)
  if (trms[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  if(length(trms)<3){return(object)}
  parent <- sys.parent()
  rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  while(inherits(rhs,"try-error") & parent > 1){                      
    parent <- parent - 1
    rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  }
  options(warn=current.warn)
  UseMethod("summary.statistics",object=rhs)
}



summary.statistics.formula <- function(object, ..., drop=FALSE, basis=NULL) {
  summary.statistics.network(object, ..., drop=drop, basis=basis)
}



summary.statistics.ergm <- function(object, ..., drop=FALSE, basis=NULL)
{
  summary.statistics.network(object$formula, ..., drop=drop, basis=basis)
}



summary.statistics.default <-
summary.statistics.matrix <- 
summary.statistics.network <- function(object,...,drop=FALSE, basis=NULL) {
  current.warn <- options()$warn
  options(warn=0)
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(object)
    formula[[2]] <- as.name("basis") # This seems irrelevant; network name
                                     # not needed by ergm.getmodel
  }else{
    formula <- object
    trms <- terms(formula)
    if(length(trms)>2){
      parent <- sys.parent()
      nw <- try(eval(trms[[2]],parent), silent = TRUE)
      while(inherits(nw,"try-error") & parent > 1){
        parent <- parent - 1
        nw <- try(eval(trms[[2]],parent), silent = TRUE)
      }
      if (inherits(nw, "try-error")) {
        stop(trms[[2]], " is not a network or network.series object")
      }
      if(class(nw) =="network.series")
        nw <- nw$networks[[1]]
      nw <- as.network(nw, ...)
    }else{
      stop("Must specify a network object")
    }
  }
  m <- ergm.getmodel(formula, nw, drop=drop)
  gs <- ergm.getglobalstats(nw, m)
  options(warn=current.warn)
  gs
}




