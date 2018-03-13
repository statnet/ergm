#  File tests/target_offset.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
library(ergm)

data(florentine)

# Based on http://tolstoy.newcastle.edu.au/R/help/04/06/0217.html by Luke Tierney.
# TODO: Rewrite to use testthat.
withWarnings <- function(expr){
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}

out <- withWarnings({ergm(flomarriage~offset(edges)+edges+degree(1)+offset(degree(0)),target.stats=summary(flomarriage~edges+degree(1)),
            offset.coef=c(0,-0.25),control=control.ergm(init=c(0,-1.47,0.462,-0.25)))})

if(length(out$warnings)!=1) stop("Unexpected number of warnings.")
summary(out$value)
mcmc.diagnostics(out$value)

set.seed(10)

out <- withWarnings({ergm(flomarriage~offset(edges)+edges+gwdegree(fix=FALSE)+degree(0)+offset(degree(1)),target.stats=summary(flomarriage~edges+gwdegree(fix=FALSE)+degree(0)), offset.coef=c(0,-0.25))})

if(length(out$warnings)!=1) stop("Unexpected number of warnings.")
summary(out$value)
mcmc.diagnostics(out$value)

