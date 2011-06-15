library(ergm)

tolerance<-0.01
n<-10
theta<--1.5

logit<-function(p) log(p/(1-p))

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(network.dyadcount(y0)-network.edgecount(y0-is.na(y1))-network.naedgecount(y1)))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

y0<-network.initialize(n)
y0<-simulate(y0~edges,theta0=theta,burnin=n^2*2)

cat("Complete data:\n")

y1<-simulate(y0~edges,theta0=theta,burnin=n^2*2)

fit.form<-ergm((y0|y1)~edges,constraints=~atleast(y0),maxit=10)
stopifnot(all.equal(form.mle(y0,y1),coef(fit.form),tolerance=tolerance,check.attributes=FALSE))

fit.diss<-ergm((y0&y1)~edges,constraints=~atmost(y0),maxit=10)
stopifnot(all.equal(diss.mle(y0,y1),coef(fit.diss),tolerance=tolerance,check.attributes=FALSE))

cat("Missing data:\n")

y1m<-network.copy(y1)
y1m[as.matrix(simulate(y0~edges,theta0=theta,burnin=n^2*2),matrix.type="edgelist")]<-NA

fit.form<-ergm((y0|y1m)~edges,constraints=~atleast(y0),maxit=10)
stopifnot(all.equal(form.mle(y0,y1m),coef(fit.form),tolerance=tolerance,check.attributes=FALSE))

fit.diss<-ergm((y0&y1m)~edges,constraints=~atmost(y0),maxit=10)
stopifnot(all.equal(diss.mle(y0,y1m),coef(fit.diss),tolerance=tolerance,check.attributes=FALSE))

