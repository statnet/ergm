#  File R/summary.statistics.network.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
#==========================================================================
# This file contains the following 5 functions for computing summary stats
#      <summary_statistics>           <summary_statistics.formula>
#      <summary.formula>              <summary_statistics.ergm>
#      <summary.statisitcs.default>   <summary_statistics.network>
#      <summary.statisitics.matrix>
#==========================================================================




#############################################################################
# Each of the <summary_statistics.X> functions and <summary.formula> checks
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



#' Calculation of network or graph statistics
#' 
#' Used to calculate the specified statistics for an observed network if its
#' argument is a formula for an \code{\link{ergm}}.  See
#' \code{\link{ergm-terms}} for more information on the statistics that may be
#' specified.
#' 
#' If \code{object} is of class \code{\link{formula}}, then
#' \code{\link[base]{summary}} may be used in lieu of \code{summary_statistics}
#' because \code{summary.formula} calls the \code{summary_statistics} function.
#' %If neither of those are given as the
#' %first argument then a \code{ergm.Cprepare} object is expected. This last
#' %option is meant for internal use.
#' The function actually cumulates the change statistics when
#' removing edges from the observed network one by one until the empty network
#' results.  Since each model term has a prespecified value (zero by default)
#' for the corresponding statistic(s) on an empty network, these change
#' statistics give the absolute statistics on the original network.
#' 
#' \code{summary.formula} for networks understands the \code{\link{lasttoggle}}
#' "API".
#' 
#' %More information can be found by looking at the documentation of
#' %\code{\link{ergm}}.
#' 
#' @aliases summary
#' @param object Either an \code{\link{formula}} object (see above) or an
#' \code{\link{ergm}} model object.  In the latter case,
#' \code{summary_statistics} is called for the \code{object$formula} object.
#' In the former case, \code{object} is of the form \code{y ~ <model terms>},
#' where \code{y} is a \code{\link[network]{network}} object or a matrix that
#' can be coerced to a \code{\link[network]{network}} object.  For the details
#' on the possible \code{<model terms>}, see \code{\link{ergm-terms}}.  To
#' create a \code{\link[network]{network}} object in , use the \code{network()}
#' function, then add nodal attributes to it using the \code{\%v\%} operator if
#' necessary.
#' @param response Name of the edge attribute whose value is to be modeled.
#' Defaults to \code{NULL} for simple presence or absence, modeled via binary
#' ERGM terms. Passing anything but \code{NULL} uses valued ERGM terms.
#' @param basis An optional \code{\link[network]{network}} object relative to
#' which the global statistics should be calculated.
#' @param \dots further arguments passed to or used by methods.
#' @return A vector of statistics measured on the network.
#' @seealso [ergm()], [network()], [ergm-terms]
#' @keywords models
#' @examples
#' 
#' #
#' # Lets look at the Florentine marriage data
#' #
#' data(florentine)
#' #
#' # test the summary_statistics function
#' #
#' summary(flomarriage ~ edges + kstar(2))
#' m <- as.matrix(flomarriage)
#' summary(m ~ edges)  # twice as large as it should be
#' summary(m ~ edges, directed=FALSE) # Now it's correct
#' 
#' @export summary_statistics
summary_statistics <- function(object, ..., basis=NULL) {
  UseMethod("summary_statistics")
}

#' @rdname summary_statistics
#' @usage \method{summary}{formula}(object, \dots)
#' @export
summary.formula <- function(object, ...){
  if(length(object)!=3 || object[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  lhs <- eval(object[[2]], envir = environment(object))
  UseMethod("summary_statistics",object=lhs)
}


## #' @describeIn summary_statistics an [ergm()] [`formula`] method.
## #' @export
## summary_statistics.formula <- function(object, ..., basis=NULL) {
##   summary_statistics.network(object, ..., basis=basis)
## }


#' @describeIn summary_statistics an [`ergm`] fit method, extracting its model from the fit.
#' @export
summary_statistics.ergm <- function(object, ..., basis=NULL)
{
  summary_statistics.network(object$formula, ..., basis=basis)
}

#' @describeIn summary_statistics a method for a [`network.list`] on the LHS of the formula.
#' @export
summary_statistics.network.list <- function(object, response=NULL, ..., basis=NULL){
  if(!is.null(basis)){
    if(inherits(basis,'network.list'))
      object[[2]] <- basis
    else stop('basis, if specified, should be the same type as the LHS of the formula (network.list, in this case).')
  }
  nwl <- eval(object[[2]], envir=environment(object))
  out<-lapply(nwl, function(nw) summary_statistics.network(object, response=response, ..., basis=nw))
  do.call(rbind,out)
}

#' @describeIn summary_statistics a method for a [`network`] on the LHS of the formula.
#' @export
summary_statistics.network <- function(object, response=NULL,...,basis=NULL) {
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

#' @describeIn summary_statistics a method for a [`matrix`] on the LHS of the formula.
#' @export
summary_statistics.matrix <- summary_statistics.network
#' @describeIn summary_statistics a fallback method.
#' @export
summary_statistics.default <- summary_statistics.network



