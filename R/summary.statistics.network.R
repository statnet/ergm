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
#      <summary_formula>           <summary_formula.formula>
#      <summary.formula>              <summary_formula.ergm>
#      <summary.statisitcs.default>   <summary_formula.network>
#      <summary.statisitics.matrix>
#==========================================================================




#############################################################################
# Each of the <summary_formula.X> functions and <summary.formula> checks
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


#' Calculation of network or graph statistics or other attributes
#' specified on a formula
#' 
#' Most generally, this function computes those summaries of the
#' object on the LHS of the formula that are specified by its RHS.  In
#' particular, if given a network as its LHS and
#' \code{\link{ergm-terms}} on its RHS, it computes the sufficient
#' statistics associated with those terms.
#' 
#' 
#' In practice, [summary.formula()] is a thin wrapper around the
#' [summary_formula()] generic, which dispatches methods based on the
#' class of the LHS of the formula.
#'
#' \code{summary.formula} for networks understands the
#' \code{\link{lasttoggle}} "API".
#' 
#' @param object A formula having as its LHS a
#'   \code{\link[network]{network}} object or a matrix that can be
#'   coerced to a \code{\link[network]{network}} object, a
#'   [`network.list`] or other methods. (See
#'   `methods('summary_formula') for the possible LHS methods.
#' @param \dots further arguments passed to or used by methods.
#' @return A vector of statistics specified in RHS of the formula.
#' @seealso [ergm()], [network()], [ergm-terms]
#' @keywords models
#' @examples
#' 
#' #
#' # Lets look at the Florentine marriage data
#' #
#' data(florentine)
#' #
#' # test the summary_formula function
#' #
#' summary(flomarriage ~ edges + kstar(2))
#' m <- as.matrix(flomarriage)
#' summary(m ~ edges)  # twice as large as it should be
#' summary(m ~ edges, directed=FALSE) # Now it's correct
#'
#' @export
summary.formula <- function(object, ...){
  summary_formula(object, ...)
}

#' Dispatching a summary function based on the class of the LHS of a
#' formula.
#' 
#' The generic [summary_formula()] (note the underscore) expects a
#' formula argument and will attempt to identify the class of the LHS
#' of the formula and dispatch to the appropriate `summary_formula`
#' method.
#' 
#' @param object A two-sided formula.
#' @param basis Optional object of the same class as the LHS of the formula, substituted in place of the LHS.
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
#' # test the summary_formula function
#' #
#' summary(flomarriage ~ edges + kstar(2))
#' m <- as.matrix(flomarriage)
#' summary(m ~ edges)  # twice as large as it should be
#' summary(m ~ edges, directed=FALSE) # Now it's correct
#' 
#' @export
summary_formula <- function(object, ..., basis=NULL) {
  if(length(object)!=3 || object[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  lhs <- eval(object[[2]], envir = environment(object))
  UseMethod("summary_formula",object=lhs)
}



## #' @describeIn summary_formula an [ergm()] [`formula`] method.
## #' @export
## summary_formula.formula <- function(object, ..., basis=NULL) {
##   UseMethod("summary_formula")
## }


#' @describeIn summary_formula an [`ergm`] fit method, extracting its model from the fit.
#' @export
summary_formula.ergm <- function(object, ..., basis=NULL)
{
  summary_formula.network(object$formula, ..., basis=basis)
}

#' @describeIn summary_formula a method for a [`network.list`] on the LHS of the formula.
#' @param response Name of the edge attribute whose value is to be modeled.
#' Defaults to \code{NULL} for simple presence or absence, modeled via binary
#' ERGM terms. Passing anything but \code{NULL} uses valued ERGM terms.
#' @export
summary_formula.network.list <- function(object, response=NULL, ..., basis=NULL){
  if(!is.null(basis)){
    if(inherits(basis,'network.list'))
      object[[2]] <- basis
    else stop('basis, if specified, should be the same type as the LHS of the formula (network.list, in this case).')
  }
  nwl <- eval(object[[2]], envir=environment(object))
  out<-lapply(nwl, function(nw) summary_formula.network(object, response=response, ..., basis=nw))
  do.call(rbind,out)
}

#' @describeIn summary_formula a method for a [`network`] on the LHS of the formula.
#' @export
summary_formula.network <- function(object, response=NULL,...,basis=NULL) {
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(object)
    formula[[2]] <- as.name("basis") # This seems irrelevant; network name
                                     # not needed by ergm_model
  }else{
    formula <- object
    nw <- ergm.getnetwork(formula)
  }
  m <- ergm_model(formula, nw, response=response, role="target",...)
  gs <- ergm.getglobalstats(nw, m, response=response)
  gs
}

#' @describeIn summary_formula a method for a [`matrix`] on the LHS of the formula.
#' @export
summary_formula.matrix <- summary_formula.network
#' @describeIn summary_formula a fallback method.
#' @export
summary_formula.default <- summary_formula.network



