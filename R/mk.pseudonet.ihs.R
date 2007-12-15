mk.pseudonet<-function(meanstats,f,y,ntoggles=length(meanstats)+1,covS=length(meanstats)^2*8,verbose=FALSE){
  if(verbose) cat("Constructing a fake network with correct (or close) meanstats:\n")
  oldwarn<-options()$warn
  on.exit(options(warn=oldwarn))
  options(warn=-1)
  
  n<-network.size(y)
  dir<-is.directed(y)
  bi<-is.bipartite(y)
  y<-network.copy(y)
  max.edges<-n*(n-1)/(2-dir)
    
  rbern.net<-function(y,p){
    y<-network.copy(y)
    y<-network.update(y,as.matrix.network.edgelist(as.network(n,density=p,directed=dir,bipartite=bi)),matrix.type="edgelist")
    y
  }
  summ.net<-function(y){
    summary(as.formula(paste("y~",f[3])))
  }
  
  if(verbose) cat("n=",n,", dir=",dir,", max.edges=",max.edges,"\n",sep="")
  
  all.edges<-as.matrix.network.edgelist(rbern.net(y,1))
  all.edges<-all.edges[order(runif(dim(all.edges)[1])),,drop=FALSE]
  
  if(verbose) cat("Estimating the covariance matrix of network statistics... ")

  dens.stats<-t(sapply(1:covS,function(i)
                       summ.net(rbern.net(y,(i-1)/(covS-1)))))
  
  wt<-robust.inverse(cov(dens.stats))
  if(verbose) cat("Finished.\n")
  
  if(verbose){
    cat("Weight matrix:\n")
    print(wt)
  }
  
  start.density<-(which.min(mahalanobis(dens.stats,meanstats,wt,inverted=TRUE))-1)/(covS-1)

  if(verbose) cat("Starting density:",start.density,"\n")
  
  decider<-function(target,cur,prop,wt)
    mahalanobis(cur,target,wt,inverted=TRUE)-mahalanobis(prop,target,wt,inverted=TRUE)>=sqrt(.Machine$double.eps)*runif(1,-1,1)

  y<-rbern.net(y,start.density)
  ms<-summ.net(y)
  y.prop<-network.copy(y)
               
  t<-0
  i<-floor(runif(ntoggles,1,max.edges+1))
  balance<-0
  last.acc<-0
  while(t<max.edges && (t<=n || t/2>=(t-last.acc))){
    if(verbose>=2) {print(summary(y))}
    
    if(all(ms==meanstats)){
      if(verbose) cat("\nNetwork with right statistics found.\n")
      return(y)
    }
    
    nt<-floor(runif(1,1,ntoggles+1))
    
    e<-all.edges[unique(i[1:nt]),,drop=FALSE]

    if(dim(e)[1]==1) y.e<-y[e[1],e[2]]
    else y.e<-y[e]

    if(exp(sum((-1)^(1-y.e))*balance/length(y.e))>runif(1)){
      if(verbose) cat(t,"/",max.edges,": i=",paste(i,collapse=","),"
      bal=",balance," ",sep="")
      if(verbose) cat("nt=",length(y.e),": ",y.e,"->",1-y.e,sep="")
      t<-t+1

      if(length(y.e)==1) y.prop[e[1],e[2]]<-1-y.e
      else y.prop[e]<-1-y.e
      balance<-balance+if(length(y.e)==1) y.prop[e[1],e[2]]-y.e else
      sum(y.prop[e]-y.e)
      ms.prop<-summ.net(y.prop)
      
      if(verbose) cat(":",ms.prop-ms)
      
      if(decider(meanstats,ms,ms.prop,wt)){
        y<-network.copy(y.prop)
        ms<-ms.prop
        last.acc<-t
        if(verbose) cat(": Acc")
      }else{
        if(verbose) cat(": Rej")
        y.prop<-network.copy(y)
        ms.prop<-ms
      }

      if(verbose)cat(": stats=",ms,"\n")
    }
    i<-(i+1:ntoggles)%%max.edges+1
  }
  return(y)
}
