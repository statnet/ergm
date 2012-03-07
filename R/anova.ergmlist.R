#  File ergm/R/anova.ergmlist.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#################################################################################
# The <anova.ergmlist> function computes an analysis of variance table for one
# or more linear model fits with the same response.
#################################################################################
anova.ergmlist <- function (object, ..., eval.loglik=FALSE, scale = 0, test = "F") 
{
  objects <- list(object, ...)
  responses <- as.character(lapply(objects, function(x) deparse(x$formula[[2]])))
  sameresp <- responses == responses[1]
  if (!all(sameresp)) {
    objects <- objects[sameresp]
    warning("Models with response ", deparse(responses[!sameresp]), 
            " removed because response differs from ", "model 1")
  }
  nmodels <- length(objects)
  if (nmodels == 1) 
    return(anova.ergm(object))
  n <- network.size(object$newnetwork)
  logl <- df <- Rdf <- rep(0, nmodels)
  for (i in 1:nmodels) {
    nodes<- network.size(objects[[i]]$newnetwork)
    n <- network.dyadcount(objects[[i]]$newnetwork)
    df[i] <- length(objects[[i]]$coef) 
    Rdf[i] <- n - df[i]
    logl[i] <- logLik(objects[[i]],eval.loglik=eval.loglik)
  }
  k <- nmodels
# k <- 1 + length(objects[[i]]$glm$coef)
#
# if (k >= 2) {
#    k <- k+1
#    if(length(object$glm$coef) > 3)
#      varlist <- attr(object$terms, "variables")
#    x <- if (n <- match("x", names(object$glm), 0))
#      object$glm[[n]]
#    else
#      model.matrix(object$glm)
#    varseq <- attr(x, "assign")
#    nvars <- max(0, varseq)
#    resdev <- resdf <- NULL
#    if(nvars>1)
#      for(i in 1:(nvars - 1))
#      {
#        fit <- glm.fit(x = x[, varseq <= i, drop = FALSE], 
#                     y = object$glm$y, weights = object$prior.weights, 
#                     start = object$glm$start, offset = object$glm$offset, 
#                     family = object$glm$family, control = object$glm$control)
#        resdev <- c(resdev, fit$deviance)
#        resdf <- c(resdf, fit$df.residual)
#      }
#
#    df <- c(0, object$glm$df.null - object$glm$df.residual, df)
#    Rdf <- c(object$glm$df.null, resdf,object$glm$df.residual, Rdf)
#    df <- n - Rdf
#    if(length(resdev>0))
#      logl <- c(-object$glm$null.deviance/2, -resdev/2,-object$glm$deviance/2, logl)
#    else logl <- c(-object$glm$null.deviance/2, -object$glm$deviance/2, logl)
#  } else {
    df <- c(0, df)
#   Rdf <- c(object$glm$df.null, Rdf)
#   logl <- c(-object$glm$null.deviance/2, logl)
    Rdf <- c(n, Rdf)
    logl <- c(-n*log(2), logl)
#  }
  pv <- pchisq(abs(2 * diff(logl)), abs(diff(df)), lower.tail = FALSE)

  table <- data.frame(c(NA, -diff(Rdf)), c(NA, diff(2 * logl)), 
                      Rdf, -2 * logl, c(NA, pv))
  variables <- lapply(objects, function(x) paste(deparse(formula(x)), 
                                                 collapse = "\n"))
  colnames(table) <- c("Df","Deviance", "Resid. Df",
                              "Resid. Dev", "Pr(>|Chisq|)")
  if (k > 2) 
    rownames(table) <- c("NULL", object$glm.names,1:nmodels)
  else
    rownames(table) <- c("NULL", 1:nmodels)

  title <- "Analysis of Variance Table\n"
  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  structure(table, heading = c(title, topnote), class = c("anova", 
                                                  "data.frame"))
}
