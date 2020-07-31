#' A heuristic check for nonidentifiability of an ERGM by examining linear dependence in the model statistics.
#'
#' @param x a matrix of MPLE-style covariates, sampled statistics, or a symmetric square matrix.
#' @param theta current model parameters (particularly for curved models); should be `NULL` if `x` is already estimating functions.
#' @param model an [`ergm_model`].
#' @param tol tolerance before declaring a statistic linearly dependent.
#' @param type the type of matrix that is `x`.
#' @param action what to do in case of non-identifiability.
#' @noRd

check_nonidentifiability <- function(x, theta, model, tol=1e-10, type=c("covariates","statistics"), action=c("warning","message","error")){
  type <- match.arg(type)
  action <- match.arg(action)
  if(!is.null(theta)) x <- ergm.estfun(x,theta,model)
  v <- switch(type,
              covariates = crossprod(x),
              statistics = cov(x)
              )
  q <- qr(v, tol=tol)
  redundant <- q$pivot[-seq_len(q$rank)]
  redundant_names <- param_names(model, canonical=FALSE, offset=FALSE)[redundant]
  if(length(redundant)){
    msg <- paste0("Model statistics ", paste.and(sQuote(redundant_names)), " are linear combinations of some set of preceding statistics at the current stage of the estimation. This may indicate that the model is nonidentifiable.")
    switch(action,
           error = stop(msg, call.=FALSE),
           warning = warning(msg, immediate.=TRUE, call.=FALSE), # Warn immediately, so the user gets the warning before the MCMC starts.
           message = message(msg)
           )
  }
  invisible(list(qr=q, redundant=redundant, redundant_names=redundant_names))
}
