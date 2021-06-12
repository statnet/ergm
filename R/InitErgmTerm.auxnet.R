#  File R/InitErgmTerm.auxnet.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
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

  x <- a$x

  list(name=name,
       coef.names=c(), dependence=FALSE,
       auxiliaries = trim_env(~ .discord.net(x, implementation="Network"), "x"))
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
       iinputs=to_ergm_Cdouble(a$x, prototype=nw),
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
       iinputs=to_ergm_Cdouble(a$x, prototype=nw),
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
       iinputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

## Exports a network that tracks the current network but excludes ties outside of the blocks specified by attrname.

InitErgmTerm..blockdiag.net <- function(nw, arglist, ...){
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

  list(name="_blockdiag_net", coef.names = c(), iinputs=nodecov, dependence=FALSE)
}

InitErgmTerm..undir.net <- function(nw, arglist, ...){
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

InitErgmTerm..subgraph.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("tailsel", "headsel"),
                      vartypes = c("numeric", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))

  tailsel <- a$tailsel
  headsel <- a$headsel

  type <-
    if(!identical(tailsel, headsel)) "bip"
    else if(is.directed(nw)) "dir" else "undir"

  tailmap <- numeric(network.size(nw))
  tailmap[tailsel] <- seq_along(tailsel)
  if(type=="bip"){
    headmap <- numeric(network.size(nw))
    headmap[headsel] <- length(tailsel) + seq_along(headsel)
  }else headmap <- numeric(0)

  TYPES <- c("dir", "undir", "bip")

  list(name=paste0("_subgraph_net"), coef.names = c(), iinputs=c(match(type, TYPES), length(tailsel), if(type=="bip") length(headsel), tailmap, headmap), dependence=FALSE)
}
