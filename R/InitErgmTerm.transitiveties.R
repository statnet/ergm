#  File R/InitErgmTerm.transitiveties.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
# This new InitErgmTerm function still needs to be tested:

#################################################################################

#' @templateVar name transitiveties
#' @title Transitive ties
#' @description This term adds one statistic, equal to the number of ties
#'   \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
#'   \eqn{i} to \eqn{j} . (Related to the `ttriple` term.)
#'
#' @usage
#' # binary: transitiveties(attr=NULL, levels=NULL)
#'
#' @param attr quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If set, all three nodes involved ( \eqn{i} , \eqn{j} , and the node on the two-path) must match
#'   on this attribute in order for \eqn{i\rightarrow j}{i-->j} to be counted.
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept triad-related
#' @concept categorical nodal attribute
InitErgmTerm.transitiveties<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "diff", "levels"),
                        vartypes = c("character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, TRUE),
                        required = c(FALSE, FALSE, FALSE))
	attrarg <- a$attrname
	levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))  
	attrarg <- a$attr
	levels <- a$levels
  }
  if (a$diff) ergm_Init_abort("diff=TRUE is not currently implemented in transitiveties")

  diff <- a$diff
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
	attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      ergm_Init_warn ("Attribute given to transitiveties() has only one value")
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

#' @templateVar name cyclicalties
#' @title Cyclical ties
#' @description This term adds one statistic, equal to the number of ties
#'   \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
#'   \eqn{j} to \eqn{i} . (Related to the `ttriple` term.)
#'
#' @usage
#' # binary: cyclicalties(attr=NULL, levels=NULL)
#' @param attr quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If set, all three nodes involved ( \eqn{i} , \eqn{j} , and the node on the two-path) must match
#'   on this attribute in order for \eqn{i\rightarrow j}{i-->j} to be counted.
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.cyclicalties<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "diff", "levels"),
                        vartypes = c("character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, TRUE),
                        required = c(FALSE, FALSE, FALSE))
	attrarg <- a$attrname
	levels <- if(!is.null(a$levels)) I(a$levels) else NULL						
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, TRUE),
                        required = c(FALSE, FALSE, FALSE))
	attrarg <- a$attr
	levels <- a$levels  
  }
  if (a$diff) ergm_Init_abort("diff=TRUE is not currently implemented in cyclicalties")

  diff <- a$diff
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
	attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      ergm_Init_warn ("Attribute given to cyclicalties() has only one value")
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

