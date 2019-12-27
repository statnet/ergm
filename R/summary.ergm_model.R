#  File R/summary.ergm_model.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

#' Evaluate network summary statistics from an initialized ergm model
#' 
#' Returns a vector of the model's statistics for a given network or
#' an empty network. This is a low-level function that should not be
#' used by end-users, but may be useful to developers.
#'
#' @param object an [`ergm_model`] object.
#' @param nw a [`network`] whose statistics are to be evaluated.
#' @template response
#' @template dotdotdot
#' 
#' @seealso [summary_formula()]
#' @keywords internal
#' @export
summary.ergm_model <- function(object, nw, response=NULL,...){
  m <- object
  if(nparam(m,canonical=TRUE)==0) return(numeric(0)) # Escape if the model has 0 statistics.
    
  NVL(response) <- nw %ergmlhs% "response"
  state <- ergm_state(nw, response=response, model=m)
  summary(state)
}

#' @describeIn ergm_state a very low-level function that calculates summary statistics associated with an [`ergm_state`] object.
#' @export
summary.ergm_state <- function(object, ...){
  state <- object

  gs <-
    if(!is.valued(state))
      .Call("network_stats_wrapper",
            state,
            PACKAGE="ergm")
    else
      .Call("wt_network_stats_wrapper",
            state,
            PACKAGE="ergm")
  
  names(gs) <- param_names(state,canonical=TRUE)

  gs
}
