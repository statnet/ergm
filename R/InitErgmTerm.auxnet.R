#  File R/InitErgmTerm.test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
InitErgmTerm..sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("mode"),
                      vartypes = c("character"),
                      defaultvalues = list("integer"),
                      required = c(FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "_isociomatrix")
  
  list(name=name,
       coef.names=c(), dependence=FALSE)
}

InitErgmTerm..discord.sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "mode"),
                      vartypes = c("network", "character"),
                      defaultvalues = list(nw, "integer"),
                      required = c(FALSE, FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "_discord_isociomatrix")
  
  list(name=name,
       coef.names=c(), dependence=FALSE,
       auxiliaries = ~ .discord.net(a$x, implementation="Network"))
}

InitErgmTerm..discord.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "implementation"),
                      vartypes = c("network,matrix", "character"),
                      defaultvalues = list(NULL, "DyadSet"),
                      required = c(TRUE, FALSE))

  impl <- match.arg(a$implementation, c("Network","DyadSet"))
  
  list(name=paste0("_discord_net_",impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..intersect.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "assume_all_toggles_in_list", "implementation"),
                      vartypes = c("network,matrix", "logical", "character"),
                      defaultvalues = list(NULL, FALSE, "DyadSet"),
                      required = c(TRUE, FALSE, FALSE))
  
  impl <- match.arg(a$implementation, c("Network","DyadSet"))

  list(name=paste0(if(a$assume_all_toggles_in_list) "_intersect_net_toggles_in_list_" else "_intersect_net_", impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..union.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "implementation"),
                      vartypes = c("network,matrix", "character"),
                      defaultvalues = list(NULL, "DyadSet"),
                      required = c(TRUE, FALSE))

  impl <- match.arg(a$implementation, c("Network","DyadSet"))
  
  list(name=paste0("_union_net_",impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

## Exports a network that tracks the current network but excludes ties outside of the blocks specified by attrname.

InitErgmTerm..blockdiag.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov,u)

  list(name="_blockdiag_net", coef.names = c(), inputs=nodecov, dependence=FALSE)
}

InitErgmTerm..undir.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = c("rule"),
                      vartypes = c("character"),
                      defaultvalues = list("weak"),
                      required = c(FALSE))
  RULES <- c("weak","strong","upper","lower")
  rule <- match.arg(a$rule, RULES)
  ruleID <- which(RULES==rule)

  list(name="_undir_net", coef.names = c(), inputs=ruleID,
       dependence=rule%in%c("weak","strong")) # Just discarding half the network does not induce dependence.
}

InitErgmTerm..subgraph.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrs", "auxnw"),
                      vartypes = c(ERGM_VATTR_SPEC, "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  # If LHS is itself an auxiliary.
  if(!is.null(a$auxnw))
  
  spec <- a$attrs
  if(!is(spec, "formula")) spec <- call("~",spec) # Embed into formula.
  if(length(spec)==2) spec <- call("~", spec[[2]], spec[[2]]) # Duplicate if one-sided.

  environment(spec) <- environment(a$attrs)

  # Obtain subformulas for tail and head.
  tailspec <- spec[-3]
  headspec <- spec[-2]

  # Obtain the boolean indicators or numeric indices. If the network
  # is bipartite in the first place, expect bipartite indices.
  bip <- nw %n% "bipartite"

  tailsel <- ergm_get_vattr(tailspec, nw, accept="numeric", bip=if(bip) "b1" else "n")
  tailname <- attr(tailsel, "name")
  
  headsel <- ergm_get_vattr(headspec, nw, accept="numeric", bip=if(bip) "b2" else "n")
  headname <- attr(headsel, "name")

  # Convert to numeric selectors.
  if(is.logical(tailsel)) tailsel <- sort(which(tailsel))
  if(is.logical(headsel)) headsel <- sort(which(headsel))

  # TODO: Check if 1-node graphs cause crashes.
  if(length(tailset)==0 || length(headset)==0) stop("Empty subgraph selected.")

  # If bipartite, remap to the whole network's indexes.
  if(bip) headsel <- headsel + bip
  
  if(bip){
    if(ult(tailsel)>bip || headset[1]<=bip)
      stop("Invalid vertex subsets selected for a bipartite graph.")
  }else{
    if(!identical(tailset,headset)){ # Rectangular selection: output bipartite.
      if(!(ult(tailset)<headset[1] || ult(headset)<tailset[1]))) stop("Vertex subsets constructing a bipartite subgraph must have disjoint ranges.")

      type <- "bipartite"
    }
    
  }
  
  list(name="_blockdiag_net", coef.names = c(), inputs=nodecov, dependence=FALSE)
}
