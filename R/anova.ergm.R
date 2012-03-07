#  File ergm/R/anova.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# The <anova.ergm> function computes an analysis of variance table for a
# single model fit
#################################################################################
anova.ergm <- function (object, ..., eval.loglik=FALSE) 
{
  if (length(list(object, ...)) > 1) 
    return(anova.ergmlist(object, ...,eval.loglik=eval.loglik))
  
  logl <- try(logLik(object,eval.loglik=eval.loglik), silent=TRUE)
  if(inherits(logl,"try-error"))
    stop(nologLik.message(deparse(substitute(object))))

  nodes<- network.size(object$newnetwork)
  n<- network.dyadcount(object$network)
  df <- length(object$coef)
  Rdf <- n - df

  k <- 1 + (length(object$mplefit$glm$coef) >= 2)
  df <- c(0, df)
  Rdf <- c(n, Rdf)
  logl <- c(-n*log(2), logl)
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
