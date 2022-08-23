#  File R/summary.statistics.network.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#==========================================================================
# This file contains the following 5 functions for computing summary stats
#      <summary_formula>           <summary_formula.formula>
#      <summary.formula>              <summary_formula.ergm>
#      <summary.statisitcs.default>   <summary_formula.network>
#      <summary.statisitics.matrix>
#==========================================================================


#' Calculation of network or graph statistics or other attributes
#' specified on a formula
#' 
#' Most generally, this function computes those summaries of the
#' object on the LHS of the formula that are specified by its RHS.  In
#' particular, if given a network as its LHS and
#' \code{\link{ergmTerm}} on its RHS, it computes the sufficient
#' statistics associated with those terms.
#' 
#' 
#' In practice, [summary.formula()] is a thin wrapper around the
#' [summary_formula()] generic, which dispatches methods based on the
#' class of the LHS of the formula.
#'
#' @aliases summary
#' @param object A formula having as its LHS a
#'   \code{\link[network]{network}} object or a matrix that can be
#'   coerced to a \code{\link[network]{network}} object, a
#'   [`network.list`], or other types to be summarized using a formula. (See
#'   `methods('summary_formula') for the possible LHS types.
#' @param \dots further arguments passed to or used by methods.
#' @return A vector of statistics specified in RHS of the formula.
#' @seealso [ergm()], [network()], [`ergmTerm`]
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
#' @seealso [ergm()], [network()], [`ergmTerm`]
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
#' @keywords internal
#' @export
summary_formula <- function(object, ..., basis=NULL) {
  if(!is.null(basis)){
    lhs <- basis
  }else{
    if(length(object)!=3 || object[[1]]!="~")
      stop ("Formula must be of form 'y ~ model'.")
    lhs <- eval_lhs.formula(object)
  }

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
#' @template response
#' @export
summary_formula.network.list <- function(object, response=NULL, ..., basis=eval_lhs.formula(object)){
  out<-lapply(basis, function(nw) summary_formula.network(object, response=response, ..., basis=nw))
  do.call(rbind,out)
}

#' @describeIn summary_formula a method for a [`network`] on the LHS of the formula.
#' @seealso [summary.ergm_model()]
#' @export
summary_formula.network <- function(object, response=NULL,...,basis=ergm.getnetwork(object)){
  ergm_preprocess_response(basis,response)
  m <- ergm_model(object, basis, ...)
  summary(m, basis)
}

#' @describeIn summary_formula a method for the semi-internal [`ergm_state`] on the LHS of the formula.
#' @export
summary_formula.ergm_state <- function(object,...,basis=NULL) {
  if(is(basis,"ergm_state")){
    s <- basis
  }else{
    s <- eval_lhs.formula(object)
  }
  summary(s)
}

#' @describeIn summary_formula a method for a [`matrix`] on the LHS of the formula.
#' @export
summary_formula.matrix <- summary_formula.network
#' @describeIn summary_formula a fallback method.
#' @export
summary_formula.default <- summary_formula.network
