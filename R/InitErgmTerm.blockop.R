## Creates a submodel that ignores any edges not within the
## blocks.
##
## FIXME: Handle curved terms.

InitErgmTerm.within.block <- function(nw, arglist, response=NULL, ...){
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

  list(name="within_block", coef.names = paste0('within.block(',m$coef.names,a$attrname,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs, auxiliaries=~.within.block(attrname))
}

## Creates a submodel that tracks the given formula.

InitErgmTerm..within.block <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)

  list(name="_within_block", coef.names = c(), inputs=nodecov, dependence=FALSE)
}
