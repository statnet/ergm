InitErgmTerm.degcrossprod<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 
  ### Construct the list to return
  list(name="degcrossprod",                            #name: required
       coef.names = "degcrossprod",                    #coef.names: required
       emptynwstats=0,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
