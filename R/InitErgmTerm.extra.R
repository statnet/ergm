#  File R/InitErgmTerm.extra.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
# These are InitErgm functions that were never included in the public release.
#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch

#' @templateVar name concurrentties
#' @title Concurrent tie count
#' @description This term adds one network statistic to the model, equal to the number of
#'   ties incident on each actor beyond the first. 
#'   This term can only be used with undirected networks.
#'
#' @usage
#' # binary: concurrentties(by=NULL, levels=NULL)
#' @param by a vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.);
#'   it functions just like the `by` argument of the `degree` term
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.concurrentties<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,bipartite=NULL,
                        varnames = c("by", "levels"),
                        vartypes = c("character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,bipartite=NULL,
                        varnames = c("by", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- a$levels
  }
  byarg <- a$by
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
	attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    if (length(u)==1)
      ergm_Init_abort ("Attribute given to concurrentties() has only one value")
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
   
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    #  See comment in d_concurrent_by_attr function
    coef.names <- paste("concurrentties",".", attrname, u, sep="")
    name <- "concurrent_ties_by_attr"
    inputs <- c(ui, nodecov)
  }else{
    coef.names <- "concurrentties"
    name <- "concurrent_ties"
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0)
}



