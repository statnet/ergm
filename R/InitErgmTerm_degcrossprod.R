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
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="degcor",                            #name: required
       coef.names = "degcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
InitErgmTerm.adegcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="adegcor",                            #name: required
       coef.names = "adegcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
InitErgmTerm.rdegcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="rdegcor",                            #name: required
       coef.names = "rdegcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
InitErgmTerm.pdegcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE) 

  el=as.matrix.network.edgelist(nw)
  deg1<-summary(nw ~ sender(base=0))[el[,1]]
  deg2<-summary(nw ~ receiver(base=0))[el[,2]]
  deg12na<-is.na(deg1)|is.na(deg2)
  deg1<-deg1[!deg12na]
  deg2<-deg2[!deg12na]
  sigma1<-(sum(deg1*deg1)-length(deg1)*(mean(deg1)^2))
  sigma2<-(sum(deg2*deg2)-length(deg2)*(mean(deg2)^2))
  sigma1 <- 0
  sigma2 <- 0
  ### Construct the list to return
  list(name="pdegcor",                            #name: required
       coef.names = "pdegcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
