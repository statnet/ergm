# These are InitErgm functions that were never included in the public release.
#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch

InitErgmTerm.concurrentties<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
                      varnames = c("byarg"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  if(!is.null(a$byarg)) {
    nodecov <- get.node.attr(nw, a$byarg, "concurrentties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrentties() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
  out <- list(name="concurrent_ties",                      #name: required
              coef.names = "concurrentties",               #coef.names: required
              minval = 0
              ) 
  if(!is.null(a$byarg)) {
    #  No covariates here, so input component 1 is arbitrary
    out$coef.names <- paste("concurrentties",".", a$byarg, u, sep="")
    out$inputs <- c(0, length(u), length(u)+length(nodecov), ui, nodecov)
    out$name="concurrent_ties_by_attr"
    # See comment in d_concurrent_ties_by_attr function
  }else{
    #  No covariates here, so input component 1 is arbitrary
    out$inputs <- c(0, 1, 0)
    out$coef.names <- "concurrentties"
  }
  out
}



