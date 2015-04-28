#  File R/anova.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
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

anova.ergm <- function (object, ..., eval.loglik=FALSE) 
{
  if (length(list(object, ...)) > 1) 
    return(anova.ergmlist(object, ...,eval.loglik=eval.loglik))
  
  logl <- try(logLik(object,eval.loglik=eval.loglik), silent=TRUE)
  if(inherits(logl,"try-error"))
    stop(nologLik.message(deparse(substitute(object))))

  nodes<- network.size(object$newnetwork)
  n<- nobs(logl)
  df <- length(object$coef)
  Rdf <- n - df

  k <- 1 + (length(object$mplefit$glm$coef) >= 2)
  df <- c(0, df)
  Rdf <- c(n, Rdf)
  logl <- c(0, logl)
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
