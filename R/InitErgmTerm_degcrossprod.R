InitErgmTerm.degcrossprod<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 
  ### Construct the list to return
  list(name="degcrossprod",                            #name: required
       coef.names = "degcrossprod",                    #coef.names: required
       inputs=2*summary(nw ~ edges),
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
InitErgmTerm.degcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(mean(alldeg*alldeg)-mean(alldeg)^2)
  ### Construct the list to return
  list(name="degcor",                            #name: required
       coef.names = "degcor",                    #coef.names: required
       inputs=sigma2*2*summary(nw ~ edges),
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
