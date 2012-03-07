#  File ergm/R/InitErgmTerm.extra.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
# These are InitErgm functions that were never included in the public release.

#########################################################
InitErgmTerm.cyclicalties<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
 ### Process the arguments
  ### Construct the list to return
  out <- list(name="cyclicalties",                      #name: required
              coef.names = "cyclicalties"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "cyclicalties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties() has only one value", call.=FALSE)
    if (!a$diff) {
      out$coef.names <- paste("cyclicalties", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("cyclicalties", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "cyclicalties"
    out$inputs <- c(0,1,0)
  }
  out$minval <- 0
  out
}

#########################################################
InitErgmTerm.cyclicalties2<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
 ### Process the arguments
  ### Construct the list to return
  out <- list(name="cyclicalties2",                      #name: required
              coef.names = "cyclicalties2"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "cyclicalties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties2() has only one value", call.=FALSE)

    if (!a$diff) {
      out$coef.names <- paste("cyclicalties2", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("cyclicalties2", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "cyclicalties2"
    out$inputs <- c(0,1,0)
  }
  out$minval <- 0
  out
}

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
  out <- list(name="concurrentties",                      #name: required
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

#########################################################
InitErgmTerm.transitiveties2<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  out <- list(name="transitiveties2",                      #name: required
              coef.names = "transitiveties2"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "transitiveties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to transitiveties2() has only one value", call.=FALSE)

    if (!a$diff) {
      out$coef.names <- paste("transitiveties2", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("transitiveties2", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "transitiveties2"
    out$inputs <- c(0,1,0)
  }
  out$minval <- 0
  out
}

InitErgmTerm.concurrentties2<-function(nw, m, arglist, ...) {
  ergm.checkdirected("concurrentties2", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("concurrentties2", arglist,
                      varnames = c("byarg"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  out <- list(name="cyclicalties",                      #name: required
              coef.names = "cyclicalties",              #coef.names: required
              minval = 0
              ) 
  if(!is.null(a$byarg)) {
    nodecov <- get.node.attr(nw, a$byarg, "concurrentties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrentties2() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)

  }
  
  out <- list(name="concurrentties2",                      #name: required
              coef.names = "concurrentties2",              #coef.names: required
              minval = 0
              ) 
  if(!is.null(a$byarg)) {
    #  No covariates here, so input component 1 is arbitrary
    out$coef.names <- paste("concurrentties2",".", a$byarg, u, sep="")
    out$inputs <- c(0, length(u), length(u)+length(nodecov), ui, nodecov)
    out$name="concurrent_ties2_by_attr"
    # See comment in d_concurrent_ties_by_attr function
  }else{
    #  No covariates here, so input component 1 is arbitrary
    out$inputs <- c(0, 1, 0)
    out$coef.names <- "concurrentties2"
  }
  out
}
