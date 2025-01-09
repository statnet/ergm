#  File R/nonidentifiability.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' A heuristic check for nonidentifiability of an ERGM by examining linear dependence in the model statistics.
#'
#' @param x a matrix of MPLE-style covariates, sampled statistics.
#' @param theta current model parameters (particularly for curved models); should be `NULL` if `x` is already estimating functions.
#' @param model an [`ergm_model`].
#' @param tol tolerance before declaring a statistic linearly dependent.
#' @param type the type of matrix that is `x`.
#' @param action,nonident_action,nonvar_action what to do in case of non-identifiable and nonvarying statistics; `action` controls both but can be overridden.
#' @noRd
check_nonidentifiability <- function(x, theta, model, tol=1e-10, type=c("covariates","statistics"), action=c("warning","message","error"), nonident_action=action, nonvar_action=action){
  type <- match.arg(type)
  nonident_action <- match.arg(nonident_action)
  nonvar_action <- match.arg(nonvar_action)
  if(!missing(theta) && !is.null(theta)) x <- ergm.estfun(x,theta,model)
  v <- switch(type,
              covariates = crossprod(x),
              statistics = cov(x)
              )
  nonvarying <- diag(v)==0
  nonvarying_names <- param_names(model, canonical=FALSE, offset=FALSE)[nonvarying]

  msg <- structure(vector("list",3), names=c("message","warning","error"))
  if(any(nonvarying)) msg[[nonvar_action]] <- paste0("Model statistics ", paste.and(sQuote(nonvarying_names)), " are not varying. This may indicate that the observed data occupies an extreme point in the sample space or that the estimation has reached a dead-end configuration.")

  v <- v[!nonvarying,!nonvarying]

  q <- qr(v, tol=tol)
  redundant <- q$pivot[-seq_len(q$rank)]
  redundant_names <- param_names(model, canonical=FALSE, offset=FALSE)[redundant]
  if(length(redundant)) msg[[nonident_action]] <- paste(c(msg[[nonident_action]], paste0("Model statistics ", paste.and(sQuote(redundant_names)), " are linear combinations of some set of preceding statistics at the current stage of the estimation. This may indicate that the model is nonidentifiable.")), collapse=" ")

  for(action in names(msg))
    if(!is.null(m <- msg[[action]]))
      switch(action,
             error = stop(m, call.=FALSE),
             warning = warning(m, immediate.=TRUE, call.=FALSE), # Warn immediately, so the user gets the warning before the MCMC starts.
             message = message(m)
           )

  invisible(list(qr=q, redundant=redundant, redundant_names=redundant_names, nonvarying=nonvarying, nonvarying_names=nonvarying_names))
}
