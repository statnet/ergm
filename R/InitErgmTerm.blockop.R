## Creates a submodel that ignores any edges not within the
## blocks.
##
## FIXME: Handle curved terms.

InitErgmTerm.OnMatch <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrname"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clist_as_num(Clist)

  gs <- ergm.emptynwstats.model(m)

  list(name="OnMatch", coef.names = paste0('OnMatch(',m$coef.names,a$attrname,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs, auxiliaries=~.blockdiag.net(a$attrname))
}

## Exports a network that tracks the current network but excludes ties outside of the blocks specified by attrname.

InitErgmTerm..blockdiag.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  ### Process the arguments
    ### Process the arguments
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
