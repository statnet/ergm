#  File R/anova.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

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
#' If any of the objects do not have estimated log-likelihoods, produces an
#' error, unless \code{eval.loglik=TRUE}.
#' 
#' @aliases anova.ergm anova.ergmlist
#' @param object,... objects of [`ergm`], usually, a result of a
#' call to [ergm()].
#' @param eval.loglik a logical specifying whether the log-likelihood will be
#' evaluated if missing.
#' @return An object of class \code{"anova"} inheriting from class
#' \code{"data.frame"}.
#' @section Warning: The comparison between two or more models will only be
#' valid if they are fitted to the same dataset. This may be a problem if there
#' are missing values and 's default of \code{na.action = na.omit} is used, and
#' [anova.ergmlist()] will detect this with an error.
#' @seealso The model fitting function [ergm()], [anova()],
#' [logLik.ergm()] for adding the log-likelihood to an existing
#' [`ergm`] object.
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
    return(anova.ergmlist(object, ..., eval.loglik=eval.loglik))
  
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
  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1), ": ", variables, sep = "", 
                   collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}

#' @rdname anova.ergm
#' @export
anova.ergmlist <- function(object, ..., eval.loglik = FALSE) {
  objects <- list(object=object, ...)
  if(!all(sapply(objects[-1], is.ergm))) stop("All arguments to ", sQuote("anova.ergm()"), " other than ", sQuote("eval.loglik="), " must be ", sQuote("ergm"), " fits.", call.=FALSE)

  responses <- sapply(objects, function(x) deparse1(x$formula[[2]]))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    objects <- objects[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from ", "model 1")
  }
  nmodels <- length(objects)
  if (nmodels == 1) 
    return(anova.ergm(object))
  logl <- df <- Rdf <- rep(0, nmodels)
  logl.null <- if(is.null(objects[[1]][["null.lik"]])) 0 else objects[[1]][["null.lik"]]
  for (i in 1:nmodels) {
    logli <- logLik(objects[[i]], eval.loglik = eval.loglik)
    n <- nobs(logli)
    df[i] <- nparam(objects[[i]], offset = FALSE) 
    Rdf[i] <- n - df[i]
    logl[i] <- logli
  }
    df <- c(0, df)
    Rdf <- c(n, Rdf)
    logl <- c(logl.null, logl)
  pv <- pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)

  table <- data.frame(c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- lapply(objects, function(x) paste(deparse(formula(x)), 
                                                 collapse = "\n"))
  colnames(table) <- c("Df","Deviance", "Resid. Df",
                              "Resid. Dev", "Pr(>|Chisq|)")
    rownames(table) <- c("NULL", 1:nmodels)

  title <- "Analysis of Deviance Table\n"
  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
