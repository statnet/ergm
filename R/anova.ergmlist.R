#  File R/anova.ergmlist.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
#################################################################################
# The <anova.ergmlist> function computes an analysis of variance table for one
# or more linear model fits with the same response.
#
# --PARAMETERS--
#   object:  an ergm object
#   ...   :  additional ergm objects. If these have a different response than
#            that of object, these are ignored. If this argument is not provided,
#            the <anova.ergm> function is used instead
#
# --RETURNED--
#   an anova object with the analysis of variance table for the considered ergms
#
#################################################################################

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
