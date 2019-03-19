## Creates a submodel that ignores any edges not within the
## blocks.

InitErgmTerm.NodematchFilter <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrname"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  f <- a$formula
  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  gs <- summary(m)

  c(list(name="OnAuxnet", coef.names = paste0('NodematchFilter(',m$coef.names,a$attrname,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs, auxiliaries=~.blockdiag.net(a$attrname)),
    passthrough.curved.ergm_model(m, function(x) paste0('NodematchFilter(',x,a$attrname,')')))
}

