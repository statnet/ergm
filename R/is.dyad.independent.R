#  File R/is.dyad.independent.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################



#' Testing for dyad-independence
#' 
#' These functions test whether an ERGM fit, a formula, or some other
#' object represents a dyad-independent model.
#' 
#' Dyad independence is determined by checking if all of the
#' constituent parts of the object (formula, ergm terms, constraints,
#' etc.) are flagged as dyad-independent.
#' 
#' @param object The object to be tested for dyadic independence.
#' @param \dots Unused at this time.
#' @return \code{TRUE} if the model implied by the object is
#'   dyad-independent; \code{FALSE} otherwise.
#' @keywords model
#' @export is.dyad.independent
is.dyad.independent<-function(object,...) UseMethod("is.dyad.independent")

#' @rdname is.dyad.independent
#' @description The method for `NULL` always returns `TRUE` by
#'   convention.
#' @export
is.dyad.independent.NULL <- function(object, ...) TRUE # By convention.

#' @describeIn ergm_model Tests whether the model is dyad-independent.
#' @param ignore_aux A flag to specify whether a dyad-dependent
#'   auxiliary should make the model dyad-dependent or should be
#'   ignored.
#' @export
is.dyad.independent.ergm_model <- function(object, ..., ignore_aux=TRUE){
  ## NB: Auxiliaries (i.e., terms without statistics) do not affect dyadic dependence.
  ! any(sapply(object$terms, function(term) (!ignore_aux || length(term$coef.names)) && (is.null(term$dependence) || term$dependence)))
}

#' @rdname is.dyad.independent
#' @template response
#' @param basis See [ergm()].
#' @export
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

  ergm_preprocess_response(nw, response)
  m<-ergm_model(object, nw, ...)
  is.dyad.independent(m, ...)
}

#' @rdname is.dyad.independent
#' @param object.obs For the [`ergm_conlist`] method, the observed data constraint.
#' @export
is.dyad.independent.ergm_conlist <- function(object, object.obs=NULL, ...){
  dind <- TRUE
  for(con in object){
    if(con$dependence) dind <- FALSE
  }
  dind && NVL3(object.obs, is.dyad.independent(.), TRUE)
}

#' @rdname is.dyad.independent
#' @param how one of `"overall"` (the default), `"terms"`, or
#'   "`space`", to specify which aspect of the ERGM is to be tested
#'   for dyadic independence.
#' @export
is.dyad.independent.ergm<-function(object, how = c("overall", "terms", "space"), ...){
  how <- match.arg(how)
  if(how %in% c("overall", "terms")) terms_dind <- NVL(object$info$terms_dind, is.dyad.independent(object$formula, basis=object$network,...))
  if(how %in% c("overall", "space")) space_dind <- NVL(object$info$space_dind, is.dyad.independent(object$constrained, object$constrained.obs))

  switch(how,
         overall = terms_dind && space_dind,
         terms = terms_dind,
         space = space_dind)
}
