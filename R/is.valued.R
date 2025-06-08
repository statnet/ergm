#  File R/is.valued.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Function to check whether an ERGM fit or some aspect of it is valued
#' @param object the object to be tested.
#' @param ... additional arguments for methods, currently unused.
#' @export
is.valued <- function(object, ...) UseMethod("is.valued")

#' @describeIn is.valued a method for [`ergm_state`] objects.
#' @export
is.valued.ergm_state <- function(object, ...){
  is.valued(object$el)
}

#' @describeIn is.valued a method for [`edgelist`] objects.
#' @export
is.valued.edgelist <- function(object, ...){
  ncol(object)>2
}

#' @describeIn is.valued a method for [`ergm`] objects.
#' @export
is.valued.ergm <- function(object, ...) object$info$valued

#' @describeIn is.valued a method for [`network`] objects that tests whether the network has been instrumented with a valued [`%ergmlhs%`] `"response"` specification, typically by [ergm_preprocess_response()]. Note that it is *not* a test for whether a network has edge attributes. This method is primarily for internal use.
#' @export
is.valued.network <- function(object, ...){
  NVL(attr(object %ergmlhs% "response", "valued"), FALSE)
}
