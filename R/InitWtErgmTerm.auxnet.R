#  File R/InitWtErgmTerm.auxnet.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
## InitWtErgmTerm..sociomatrix<-function(nw, arglist, ...) {
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("mode"),
##                       vartypes = c("character"),
##                       defaultvalues = list("integer"),
##                       required = c(FALSE))

##   mode <- match.arg(a$mode, c("integer"))
##   name <- switch(mode,
##                  integer = "_isociomatrix")

##   list(name=name,
##        coef.names=c(), dependence=FALSE)
## }

## InitWtErgmTerm..discord.sociomatrix<-function(nw, arglist, ...) {
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("x", "mode"),
##                       vartypes = c("network", "character"),
##                       defaultvalues = list(nw, "integer"),
##                       required = c(FALSE, FALSE))

##   mode <- match.arg(a$mode, c("integer"))
##   name <- switch(mode,
##                  integer = "_discord_isociomatrix")

##   x <- a$x

##   list(name=name,
##        coef.names=c(), dependence=FALSE,
##        auxiliaries = trim_env(~ .discord.net(x, implementation="Network"), "x"))
## }

## InitWtErgmTerm..discord.net<-function(nw, arglist, ...) {
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("x", "implementation"),
##                       vartypes = c("network,matrix", "character"),
##                       defaultvalues = list(NULL, "DyadSet"),
##                       required = c(TRUE, FALSE))

##   impl <- match.arg(a$implementation, c("Network","DyadSet"))

##   list(name=paste0("_discord_net_",impl),
##        coef.names=c(),
##        iinputs=to_ergm_Cdouble(a$x, prototype=nw),
##        dependence=FALSE)
## }

## InitWtErgmTerm..intersect.net<-function(nw, arglist, ...) {
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("x", "assume_all_toggles_in_list", "implementation"),
##                       vartypes = c("network,matrix", "logical", "character"),
##                       defaultvalues = list(NULL, FALSE, "DyadSet"),
##                       required = c(TRUE, FALSE, FALSE))

##   impl <- match.arg(a$implementation, c("Network","DyadSet"))

##   list(name=paste0(if(a$assume_all_toggles_in_list) "_intersect_net_toggles_in_list_" else "_intersect_net_", impl),
##        coef.names=c(),
##        iinputs=to_ergm_Cdouble(a$x, prototype=nw),
##        dependence=FALSE)
## }

## InitWtErgmTerm..union.net<-function(nw, arglist, ...) {
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("x", "implementation"),
##                       vartypes = c("network,matrix", "character"),
##                       defaultvalues = list(NULL, "DyadSet"),
##                       required = c(TRUE, FALSE))

##   impl <- match.arg(a$implementation, c("Network","DyadSet"))

##   list(name=paste0("_union_net_",impl),
##        coef.names=c(),
##        iinputs=to_ergm_Cdouble(a$x, prototype=nw),
##        dependence=FALSE)
## }

## ## Exports a network that tracks the current network but excludes ties outside of the blocks specified by attrname.

## InitWtErgmTerm..blockdiag.net <- function(nw, arglist, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c("attrname"),
##                       vartypes = c("character"),
##                       defaultvalues = list(NULL),
##                       required = c(TRUE))
##   nodecov <-
##     if(length(a$attrname)==1)
##       get.node.attr(nw, a$attrname)
##     else{
##       do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
##     }
##   u <- sort(unique(nodecov))
##   nodecov <- match(nodecov,u)

##   list(name="_blockdiag_net", coef.names = c(), iinputs=nodecov, dependence=FALSE)
## }

InitWtErgmTerm..undir.net <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm..undir.net
  term <- f(...)
  term$name <- "_wtundir_net"
  term
}

InitWtErgmTerm..subgraph.net <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm..subgraph.net
  term <- f(...)
  term$name <- "_wtsubgraph_net"
  term
}
