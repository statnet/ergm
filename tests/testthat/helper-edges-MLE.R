logit<-function(p) log(p/(1-p))

MLE.tools <- new.env()

MLE.tools$mk.missnet<-function(n,d,m,directed=TRUE,bipartite=FALSE){
  y<-network.initialize(n, directed=directed, bipartite=bipartite)
  y<-simulate(y~edges, coef=logit(d), control=control.simulate(MCMC.burnin=2*n^2))
  if(m>0){
    y.miss<-simulate(y~edges, coef=logit(m))
    y[as.edgelist(y.miss)]<-NA
  }
  y
}

MLE.tools$edges.theta<-function(y){
  e<-network.edgecount(y)
  d<-network.dyadcount(y)
  logit(e/d)
}

MLE.tools$edges.llk<-function(y, theta=NULL, e=NULL){
  e<-NVL(e, network.edgecount(y))
  d<-network.dyadcount(y)
  NVL(theta) <- logit(e/d)
  e*theta - d*log1p(exp(theta))
}
