#  File R/InitErgmTerm.extra.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
# These are InitErgm functions that were never included in the public release.
#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch

InitErgmTerm.concurrentties<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,bipartite=NULL,
                      varnames = c("by", "levels"),
                      vartypes = c("character", "character,numeric,logical"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  byarg <- a$by
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "concurrentties")
    u <- NVL(a$levels, sort(unique(nodecov)))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrentties() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
   
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    #  See comment in d_concurrent_by_attr function
    coef.names <- paste("concurrentties",".", byarg, u, sep="")
    name <- "concurrent_ties_by_attr"
    inputs <- c(ui, nodecov)
  }else{
    coef.names <- "concurrentties"
    name <- "concurrent_ties"
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0)
}



