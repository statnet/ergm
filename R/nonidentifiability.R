#  File R/nonidentifiability.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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

  lindep <- ergm_lindep(v, tol = tol)

  if (nrow(lindep))
    msg[[nonident_action]] <-
      paste(c(msg[[nonident_action]],
              paste0("The following linear dependence has been detected among the model statistics:\n",
                     paste("  ", format(lindep), collapse = "\n"),
                     "\nThis may indicate that the model is nonidentifiable.")),
            collapse = " ")

  for(action in names(msg))
    if(!is.null(m <- msg[[action]]))
      switch(action,
             error = stop(m, call.=FALSE),
             warning = warning(m, immediate.=TRUE, call.=FALSE), # Warn immediately, so the user gets the warning before the MCMC starts.
             message = message(m)
           )

  invisible(list(lindep = lindep, nonvarying = nonvarying, nonvarying_names = nonvarying_names))
}


ergm_lindep <- function(m, tol) {
  n <- min(dim(as.matrix(m)))
  qr <- qr(m, tol = tol)
  r <- qr.R(qr)
  rindep <- r[seq_len(qr$rank), seq_len(qr$rank)]
  rdep <- r[seq_len(qr$rank), -seq_len(qr$rank), drop = FALSE]

  coefs <- matrix(NA, n - qr$rank, n, dimnames = list(NULL, colnames(m)))

  coefs[, qr$pivot[seq_len(qr$rank)]] <- t(backsolve(rindep, rdep))
  coefs[, qr$pivot[-seq_len(qr$rank)]] <- -diag(ncol(rdep))

  structure(coefs, class = "ergm_lindep")
}

#' @noRd
#' @importFrom MASS fractions
#' @export
format.ergm_lindep <- function(x, ...) {
  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }

  format_term <- function(coeff, var_name) {
    frac <- as.character(fractions(abs(coeff)))
    var_name <- quote_var_name(var_name)
    if (frac == "1") var_name
    else if (frac == "0") character(0)
    else paste0(frac, " * ", var_name)
  }

  apply(x, 1, function(row) {
    lhs <- c()
    rhs <- c()
    for (i in seq_along(row)) {
      coeff <- row[i]
      var_name <- colnames(x)[i]
      term <- format_term(coeff, var_name)
      if (coeff > 0) {
        lhs <- c(lhs, term)
      } else {
        rhs <- c(rhs, term)
      }
    }

    # Shorter set goes on the LHS.
    if (length(lhs) > length(rhs)) {
      tmp <- lhs
      lhs <- rhs
      rhs <- tmp
    }

    lhs_str <- if (length(lhs) > 0) paste(lhs, collapse = " + ") else "CONSTANT"
    rhs_str <- if (length(rhs) > 0) paste(rhs, collapse = " + ") else "CONSTANT"
    paste(lhs_str, "=", rhs_str)
  })
}

#' @noRd
#' @export
print.ergm_lindep <- function(x, ...) {
  cat(strwrppst("Detected linear dependence in ERGM sufficient statistics or estimating functions:"), "\n")
  cat(strwrppst(format(x, ...), exdent = 2, parsep = "\n"))
}

#' @describeIn ergm Extract a matrix of detected linear dependence
#'   among the model's sufficient statistics or estimating functions
#'   (if curved). Each row, if any, contains coefficients for a linear
#'   combination of the statistics that results in a constant. These
#'   are pretty-printed as a series of equations.
#'
#' @export
alias.ergm <- function(object, ...) object$lindep
