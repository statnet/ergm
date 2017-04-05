
InitWtErgmTerm.b <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form<-match.arg(a$form,c("sum","nonzero"))

  f <- a$formula
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw,...)

  if(!is.dyad.independent(m) && form=="sum") stop("Only dyad-independent binary terms can be imported with form 'sum'.")
  
  Clist <- ergm.Cprepare(nw, m)
  inputs <- pack.Clistasnum(Clist)
  
  gs <- rep(0, Clist$nstats)

  if(form=="nonzero")
    gs <- ergm.emptynwstats(m)

  list(name=paste("import_binary_term",form,sep="_"), coef.names = paste0(form,'(',m$coef.names,')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs, auxiliaries=if(form=="nonzero") ~.binary.nonzero.net)
}

InitWtErgmTerm..binary.nonzero.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_binary_nonzero_net", depenence=FALSE)
}
