#  File R/InitErgmTerm.transitiveties.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
# This new InitErgmTerm function still needs to be tested:

#################################################################################
InitErgmTerm.transitiveties<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname", "diff", "levels"),
                      vartypes = c("character", "logical", "character,numeric,logical"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE))
  if (a$diff) stop("diff=TRUE is not currently implemented in transitiveties")
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "transitiveties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      warning ("Attribute given to transitiveties() has only one value", call.=FALSE)
    if (!diff) {
      coef.names <- paste("transitiveties",attrname,sep=".")
      inputs <- c(nodecov)
     } else { 
       coef.names <- paste("transitiveties",attrname, u, sep=".")
       inputs <- c(ui, nodecov)
       attr(inputs, "ParamsBeforeCov") <- length(ui)
     }
  }else{
    coef.names <- "transitiveties"
    inputs <- NULL
  }
  list(name="transitiveties", coef.names=coef.names, inputs=inputs, minval=0)
}

#################################################################################
InitErgmTerm.cyclicalties<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname", "diff", "levels"),
                      vartypes = c("character", "logical", "character,numeric,logical"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE))
  if (a$diff) stop("diff=TRUE is not currently implemented in cyclicalties")
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "cyclicalties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      warning ("Attribute given to cyclicalties() has only one value", call.=FALSE)
    if (!diff) {
      coef.names <- paste("cyclicalties",attrname,sep=".")
      inputs <- c(nodecov)
     } else { 
       coef.names <- paste("cyclicalties",attrname, u, sep=".")
       inputs <- c(ui, nodecov)
       attr(inputs, "ParamsBeforeCov") <- length(ui)
     }
  }else{
    coef.names <- "cyclicalties"
    inputs <- NULL
  }
  list(name="cyclicalties", coef.names=coef.names, inputs=inputs, minval=0)
}

