#  File R/anova.ergm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
################################################################################
# The <anova.ergm> function computes an analysis of variance table for a
# single model fit
#
# --PARAMETERS--
#   object:  an ergm object
#   ...   :  additional ergm objects. If this argument is provided,
#            the <anova.ergmlist> function is used instead
#
#
# --RETURNED--
#   an anova object with the analysis of variance table for the given ergm
#
#################################################################################



#' ANOVA for ERGM Fits
#' 
#' Compute an analysis of variance table for one or more ERGM fits.
#' 
#' Specifying a single object gives a sequential analysis of variance table for
#' that fit.  That is, the reductions in the residual sum of squares as each
#' term of the formula is added in turn are given in the rows of a table, plus
#' the residual sum of squares.
#' 
#' The table will contain F statistics (and P values) comparing the mean square
#' for the row to the residual mean square.
#' 
#' If more than one object is specified, the table has a row for the residual
#' degrees of freedom and sum of squares for each model.  For all but the first
#' model, the change in degrees of freedom and sum of squares is also given.
#' (This only make statistical sense if the models are nested.)  It is
#' conventional to list the models from smallest to largest, but this is up to
#' the user.
#' 
#' Optionally the table can include test statistics.  Normally the F statistic
#' is most appropriate, which compares the mean square for a row to the
#' residual sum of squares for the largest model considered.  If \code{scale}
#' is specified chi-squared tests can be used. Mallows' \eqn{C_p}{Cp} statistic
#' is the residual sum of squares plus twice the estimate of
#' \eqn{\sigma^2}{sigma^2} times the residual degrees of freedom.
#' 
#' If any of the objects do not have estimated log-likelihoods, produces an
#' error, unless \code{eval.loglik=TRUE}.
#' 
#' @aliases anova.ergm anova.ergmlist
#' @param object,... objects of class \code{\link{ergm}}, usually, a result of a
#' call to \code{\link{ergm}}.
#' @param eval.loglik a logical specifying whether the log-likelihood will be
#' evaluated if missing.
#' @return An object of class \code{"anova"} inheriting from class
#' \code{"data.frame"}.
#' @section Warning: The comparison between two or more models will only be
#' valid if they are fitted to the same dataset. This may be a problem if there
#' are missing values and 's default of \code{na.action = na.omit} is used, and
#' \code{\link{anova.ergmlist}} will detect this with an error.
#' @seealso The model fitting function \code{\link{ergm}}, \code{\link{anova}},
#' \code{\link{logLik.ergm}} for adding the log-likelihood to an existing
#' \code{\link[=ergm.object]{ergm}} object.
#' @keywords regression models
#' @examples
#' 
#' data(molecule)
#' molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#' fit0 <- ergm(molecule ~ edges)
#' anova(fit0)
#' fit1 <- ergm(molecule ~ edges + nodefactor("atomic type"))
#' anova(fit1)
#' 
#' fit2 <- ergm(molecule ~ edges + nodefactor("atomic type") +  gwesp(0.5,
#'   fixed=TRUE), eval.loglik=TRUE) # Note the eval.loglik argument.
#' anova(fit0, fit1)
#' anova(fit0, fit1, fit2)
#' 
#' @export
anova.ergm <- function (object, ..., eval.loglik=FALSE) 
{
  if (length(list(object, ...)) > 1) 
    return(anova.ergmlist(object, ...,eval.loglik=eval.loglik))
  
  logl <- try(logLik(object,eval.loglik=eval.loglik), silent=TRUE)
  if(inherits(logl,"try-error"))
    stop(NO_LOGLIK_MESSAGE)

  n<- nobs(logl)
  df <- nparam(object, offset = FALSE)
  Rdf <- n - df
  logl.null <- if(is.null(object$null.lik)) 0 else object$null.lik

  df <- c(0, df)
  Rdf <- c(n, Rdf)
  logl <- c(logl.null, logl)
  pv <- pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)
  table <- data.frame(c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- paste(deparse(formula(object)), collapse = "\n")
  colnames(table) <- c("Df", "Deviance", "Resid. Df", "Resid. Dev", 
                       "Pr(>|Chisq|)")
    rownames(table) <- c("NULL", "Model 1:")
  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1), ": ", variables, sep = "", 
                   collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
