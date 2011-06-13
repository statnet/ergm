library(ergm)

theta0err<--1 # Perturbation in the initial values
maxit<-20 # Maximum number of iterations
tolerance<-0.01 # Result must be within 1% of truth.

n<-20 # Number of nodes
b<-3 # Bipartite split

d<-.1 # Density
m<-.1 # Missingness rate

logit<-function(p) log(p/(1-p))

cat("n=",n,", density=",d,", missing=",m,"\n",sep="")
mk.missnet<-function(n,d,m,directed=TRUE,bipartite=0){
  y<-network.initialize(n,directed=directed,bipartite=bipartite)
  y<-simulate(y~edges, theta0=logit(d),burnin=2*n^2)
  if(m>0){
    y.miss<-simulate(y~edges, theta0=logit(m))
    y[as.matrix(y.miss,matrix.type="edgelist")]<-NA
  }
  y
}

correct.edges.theta<-function(y){
  e<-summary(y~edges)
  d<-network.dyadcount(y)
  m<-network.naedgecount(y)

  logit(e/d)
}


run.miss.test<-function(y){
  truth<-correct.edges.theta(y)
  cat("Correct estimate =",truth,"\n")
  
  mplefit<-ergm(y~edges)
  mpleOK<-all.equal(truth,coef(mplefit),check.attributes=FALSE,tolerance=tolerance)
  cat("MPLE estimate =",coef(mplefit),if(isTRUE(mpleOK)) "OK" else mpleOK,"\n")

  mcmcfit<-ergm(y~edges,control=control.ergm(force.mcmc=TRUE),theta0=truth+theta0err,interval=ceiling(n^(3/2)),maxit=maxit)
  mcmcOK<-all.equal(truth,coef(mcmcfit),check.attributes=FALSE,tolerance=tolerance)
  cat("MCMCMLE estimate =",coef(mcmcfit),if(isTRUE(mcmcOK)) "OK" else mcmcOK,"\n")
  
  return(isTRUE(mpleOK) && isTRUE(mcmcOK))
}

# Directed
cat("\n\nDirected Network\n")
y<-mk.missnet(n,d,m,TRUE,0)
stopifnot(run.miss.test(y))

# Undirected
cat("\n\nUndirected Network\n")
y<-mk.missnet(n,d,m,FALSE,0)
stopifnot(run.miss.test(y))  

# Bipartite Undirected
cat("\n\nBipartite Undirected Network\n")
y<-mk.missnet(n,d,m,FALSE,b)
stopifnot(run.miss.test(y))

# Add the curved+missing test here for now

n <- 50
y <- network.initialize(n,directed=FALSE) # Create an empty network
y <- simulate(y~edges, theta0=logit(0.12),burnin=2*n^2)
y.miss <- simulate(y~edges, theta0=logit(0.1))
y[as.matrix(y.miss,matrix.type="edgelist")] <- NA

cat("Network statistics:\n")
print(summary(y~edges+gwesp(0.5)))
truth<-correct.edges.theta(y)
cat("Correct estimate =",truth,"\n")

mcmcfit<-ergm(y~edges+gwesp(0.5),maxit=5)
summary(mcmcfit)
stopifnot(abs(coef(mcmcfit)[1]-truth)/sqrt(mcmcfit$covar[1])<2)
