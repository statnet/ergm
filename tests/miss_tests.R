library(ergm)

n<-100
d<-.1
m<-0
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

  logit(e/(d-m))
}

# Directed
y<-mk.missnet(n,d,m,TRUE,0)
truth<-correct.edges.theta(y)
cat("Directed test. Correct coefficient=",truth,".\n",sep="")

mplefit<-ergm(y~edges)
print(mplefit)
mcmcfit<-ergm(y~edges,control=control.ergm(force.mcmc=TRUE),theta0=truth+sample(c(-1,1),1),interval=ceiling(n^(3/2)),maxit=20)
print(mcmcfit)

# Undirected
y<-mk.missnet(n,d,m,FALSE,0)
truth<-correct.edges.theta(y)
cat("Undirected test. Correct coefficient=",truth,".\n",sep="")

mplefit<-ergm(y~edges)
print(mplefit)
mcmcfit<-ergm(y~edges,control=control.ergm(force.mcmc=TRUE),theta0=truth+sample(c(-1,1),1),interval=ceiling(n^(3/2)),maxit=20)
print(mcmcfit)

# Bipartite Undirected
y<-mk.missnet(n,d,m,FALSE,3)
truth<-correct.edges.theta(y)
cat("Bipartite undirected test. Correct coefficient=",truth,".\n",sep="")

mplefit<-ergm(y~edges)
print(mplefit)
mcmcfit<-ergm(y~edges,control=control.ergm(force.mcmc=TRUE),theta0=truth+sample(c(-1,1),1),interval=ceiling(n^(3/2)),maxit=20)
print(mcmcfit)

