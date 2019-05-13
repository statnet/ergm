#  File tests/termTests.valued.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

library(ergm)

tst <- function(truth, test){
  stopifnot(isTRUE(all.equal(truth,test,check.attributes=FALSE)))
}

n <- 129
n1 <- 100
n2 <- n-n1
nz <- 60
pint <- .5

f <- rep(letters[1:3], length.out=n)
f1 <- f[1:n1]
f2 <- f[-(1:n1)]
q <- rep(3:1, length.out=n)
q1 <- q[1:n1]
q2 <- q[-(1:n1)]

# a bipartite nw
set.seed(143)
b1 <- floor(runif(nz, 1,n1+1))
b2 <- floor(runif(nz, n1+1, n+1))
exbip.el <- unique(cbind(b1,b2))
v <- runif(nrow(exbip.el), -4, 4)
v <- ifelse(rbinom(length(v),1,pint), round(v), v)
exbip.el <- cbind(exbip.el,v)
attr(exbip.el, "n") <- n
bipnw <- as.network(exbip.el, matrix.type="edgelist", bipartite=n1, directed=FALSE, ignore.eval=FALSE, names.eval="w")
bipnw %v% "f" <- f
bipnw %v% "q" <- q
bipnw %n% "e" <- bipe <- matrix(rnorm(n1*n2), n1,n2)
bipm <- as.matrix(bipnw, a="w")
bipvt <- c(0, sort(v)[6], sort(v)[length(v)-6], runif(1, -4, 4))

# a directed nw
set.seed(143)
t <- floor(runif(nz, 1, n+1))
h <- floor(runif(nz, 1, n+1))
exdir.el <- unique(cbind(t,h))
v <- runif(nrow(exdir.el), -4, 4)
v <- ifelse(rbinom(length(v),1,pint), round(v), v)
exdir.el <- cbind(exdir.el,v)
attr(exdir.el, "n") <- n
dirnw <- as.network(exdir.el, matrix.type="edgelist", directed=TRUE, ignore.eval=FALSE, names.eval="w")
dirnw %v% "f" <- f
dirnw %v% "q" <- q
dirnw %n% "e" <- dire <- matrix(rnorm(n*n), n,n)
dirm <- as.matrix(dirnw, a="w")
diag(dirm) <- NA
dirvt <- c(0, sort(v)[6], sort(v)[length(v)-6], runif(1, -4, 4))

# an undirected nw
set.seed(143)
t <- floor(runif(nz, 1, n+1))
h <- floor(runif(nz, 1, n+1))
exund.el <- t(apply(unique(cbind(t,h)),1,range))
v <- runif(nrow(exund.el), -4, 4)
v <- ifelse(rbinom(length(v),1,pint), round(v), v)
exund.el <- cbind(exund.el,v)
attr(exund.el, "n") <- n
undnw <- as.network(exund.el, matrix.type="edgelist", directed=FALSE, ignore.eval=FALSE, names.eval="w")
undnw %v% "f" <- f
undnw %v% "q" <- q
undnw %n% "e" <- matrix(rnorm(n*n), n,n)
undnw %n% "e" <- unde <- undnw %n% "e" + t(undnw %n% "e")
undm <- as.matrix(undnw, a="w")
diag(undm) <- NA
undvt <- c(0, sort(v)[6], sort(v)[length(v)-6], runif(1, -4, 4))


# absdiff
tst(sum(abs(outer(q,q,"-"))*dirm,na.rm=TRUE), summary(dirnw ~ absdiff("q"), response="w"))
tst(sum(abs(outer(q,q,"-")^2)*dirm,na.rm=TRUE), summary(dirnw ~ absdiff(~q,pow=2), response="w"))
tst(sum(abs(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), summary(dirnw ~ absdiff(function(x) x %v% "q", form="nonzero"), response="w"))
tst(sum(abs(outer(q,q,"-")^2)*(dirm!=0),na.rm=TRUE), summary(dirnw ~ absdiff("q",pow=2, form="nonzero"), response="w"))

tst(sum(abs(outer(q,q,"-"))*undm,na.rm=TRUE)/2, summary(undnw ~ absdiff("q"), response="w"))
tst(sum(abs(outer(q,q,"-")^2)*undm,na.rm=TRUE)/2, summary(undnw ~ absdiff(~q,pow=2), response="w"))
tst(sum(abs(outer(q,q,"-"))*(undm!=0),na.rm=TRUE)/2, summary(undnw ~ absdiff(function(x) x %v% "q", form="nonzero"), response="w"))
tst(sum(abs(outer(q,q,"-")^2)*(undm!=0),na.rm=TRUE)/2, summary(undnw ~ absdiff("q",pow=2, form="nonzero"), response="w"))

tst(sum(abs(outer(q1,q2,"-"))*bipm,na.rm=TRUE), summary(bipnw ~ absdiff("q"), response="w"))
tst(sum(abs(outer(q1,q2,"-")^2)*bipm,na.rm=TRUE), summary(bipnw ~ absdiff(~q,pow=2), response="w"))
tst(sum(abs(outer(q1,q2,"-"))*(bipm!=0),na.rm=TRUE), summary(bipnw ~ absdiff(function(x) x %v% "q", form="nonzero"), response="w"))
tst(sum(abs(outer(q1,q2,"-")^2)*(bipm!=0),na.rm=TRUE), summary(bipnw ~ absdiff("q",pow=2, form="nonzero"), response="w"))

# absdiffcat
diffs <- sort(unique(c(abs(outer(q,q,"-")))))
diffs <- diffs[diffs!=0]
for(base in c(0, seq_along(diffs))){
  keep <- if(all(base==0)) seq_along(diffs) else seq_along(diffs)[-base]
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*dirm,na.rm=TRUE)), summary(dirnw ~ absdiffcat("q",levels=keep), response="w"))
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*(dirm!=0),na.rm=TRUE)), summary(dirnw ~ absdiffcat(~q,base=base, form="nonzero"), response="w"))
  
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*undm,na.rm=TRUE))/2, summary(undnw ~ absdiffcat("q",levels=keep), response="w"))
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*(undm!=0),na.rm=TRUE))/2, summary(undnw ~ absdiffcat(function(x) x %v% "q",levels=keep, form="nonzero"), response="w"))
  
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q1,q2,"-"))==x)*bipm,na.rm=TRUE)), summary(bipnw ~ absdiffcat("q", levels=keep), response="w"))
  tst(sapply(diffs[keep], function(x) sum((abs(outer(q1,q2,"-"))==x)*(bipm!=0),na.rm=TRUE)), summary(bipnw ~ absdiffcat(~q, base=base, form="nonzero"), response="w"))
}

# atleast
for(v in dirvt) tst(sum(dirm >= v,na.rm=TRUE), summary(dirnw ~ atleast(v), response="w"))
for(v in undvt) tst(sum(undm >= v,na.rm=TRUE)/2, summary(undnw ~ atleast(v), response="w"))
for(v in bipvt) tst(sum(bipm >= v,na.rm=TRUE), summary(bipnw ~ atleast(v), response="w"))

# atmost
for(v in dirvt) tst(sum(dirm <= v,na.rm=TRUE), summary(dirnw ~ atmost(v), response="w"))
for(v in undvt) tst(sum(undm <= v,na.rm=TRUE)/2, summary(undnw ~ atmost(v), response="w"))
for(v in bipvt) tst(sum(bipm <= v,na.rm=TRUE), summary(bipnw ~ atmost(v), response="w"))

# b1cov
tst(sum(q1*bipm,na.rm=TRUE), summary(bipnw ~ b1cov("q"), response="w"))
tst(c(sum(q1*bipm,na.rm=TRUE),sum(q1^2*bipm,na.rm=TRUE)), summary(bipnw ~ b1cov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q1*(bipm!=0),na.rm=TRUE), summary(bipnw ~ b1cov(~q, form="nonzero"), response="w"))
tst(c(sum(q1*(bipm!=0),na.rm=TRUE),sum(q1^2*(bipm!=0),na.rm=TRUE)), summary(bipnw ~ b1cov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

# b1factor
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(sapply(sort(unique(f1))[keep], function(x) sum((f1==x)*bipm,na.rm=TRUE)), summary(bipnw ~ b1factor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f1))[keep], function(x) sum((f1==x)*(bipm!=0),na.rm=TRUE)), summary(bipnw ~ b1factor(~f, base=base, form="nonzero"), response="w"))
}

# b1sociality
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(apply(bipm, 1, sum)[keep], summary(bipnw ~ b1sociality(nodes=keep), response="w"))
  tst(apply(bipm!=0, 1, sum)[keep], summary(bipnw ~ b1sociality(nodes=keep, form="nonzero"), response="w"))
}


# b2cov
tst(sum(q2*t(bipm),na.rm=TRUE), summary(bipnw ~ b2cov("q"), response="w"))
tst(c(sum(q2*t(bipm),na.rm=TRUE),sum(q2^2*t(bipm),na.rm=TRUE)), summary(bipnw ~ b2cov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q2*t(bipm!=0),na.rm=TRUE), summary(bipnw ~ b2cov(function(x) x %v% "q", form="nonzero"), response="w"))
tst(c(sum(q2*t(bipm!=0),na.rm=TRUE),sum(q2^2*t(bipm!=0),na.rm=TRUE)), summary(bipnw ~ b2cov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

# b2factor
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(sapply(sort(unique(f2))[keep], function(x) sum((f2==x)*t(bipm),na.rm=TRUE)), summary(bipnw ~ b2factor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f2))[keep], function(x) sum((f2==x)*t(bipm!=0),na.rm=TRUE)), summary(bipnw ~ b2factor(~f, base=base, form="nonzero"), response="w"))
}

# b2sociality
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(apply(bipm, 2, sum)[keep], summary(bipnw ~ b2sociality(nodes=keep), response="w"))
  tst(apply(bipm!=0, 2, sum)[keep], summary(bipnw ~ b2sociality(nodes=keep, form="nonzero"), response="w"))
}

# edgecov
tst(sum(dire*dirm,na.rm=TRUE), summary(dirnw ~ edgecov("e"), response="w"))
tst(sum(dire*(dirm!=0),na.rm=TRUE), summary(dirnw ~ edgecov("e", form="nonzero"), response="w"))

tst(sum(unde*undm,na.rm=TRUE)/2, summary(undnw ~ edgecov("e"), response="w"))
tst(sum(unde*(undm!=0),na.rm=TRUE)/2, summary(undnw ~ edgecov("e", form="nonzero"), response="w"))

tst(sum(bipe*bipm,na.rm=TRUE), summary(bipnw ~ edgecov("e"), response="w"))
tst(sum(bipe*(bipm!=0),na.rm=TRUE), summary(bipnw ~ edgecov("e", form="nonzero"), response="w"))

# edges
tst(sum((dirm!=0),na.rm=TRUE), summary(dirnw ~ edges, response="w"))
tst(sum((undm!=0),na.rm=TRUE)/2, summary(undnw ~ edges, response="w"))
tst(sum((bipm!=0),na.rm=TRUE), summary(bipnw ~ edges, response="w"))

# nonzero
tst(sum((dirm!=0),na.rm=TRUE), summary(dirnw ~ nonzero, response="w"))
tst(sum((undm!=0),na.rm=TRUE)/2, summary(undnw ~ nonzero, response="w"))
tst(sum((bipm!=0),na.rm=TRUE), summary(bipnw ~ nonzero, response="w"))

# diff

posonly <- function(x) pmax(x, 0)
negonly <- function(x) pmin(x, 0)

for(dd in c("t-h", "h-t")){
  for(sa in c("identity", "abs", "posonly", "negonly")){
    saf <- get(sa)
    ddf <- switch(dd, `t-h`=identity, `h-t`=function(x) -x)

    df <- function(x) saf(ddf(x))
    
    tst(sum(df(outer(q,q,"-"))*dirm,na.rm=TRUE), summary(dirnw ~ diff("q", dir=dd, sign.action=sa), response="w"))
    tst(sum(df(outer(q,q,"-"))^2*dirm,na.rm=TRUE), summary(dirnw ~ diff(~q,pow=2, dir=dd, sign.action=sa), response="w"))
    tst(sum(df(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), summary(dirnw ~ diff(function(x) x %v% "q", dir=dd, sign.action=sa, form="nonzero"), response="w"))
    tst(sum(df(outer(q,q,"-"))^2*(dirm!=0),na.rm=TRUE), summary(dirnw ~ diff("q",pow=2, dir=dd, sign.action=sa, form="nonzero"), response="w"))
  }
}

# greaterthan
for(v in dirvt) tst(sum(dirm > v,na.rm=TRUE), summary(dirnw ~ greaterthan(v), response="w"))
for(v in undvt) tst(sum(undm > v,na.rm=TRUE)/2, summary(undnw ~ greaterthan(v), response="w"))
for(v in bipvt) tst(sum(bipm > v,na.rm=TRUE), summary(bipnw ~ greaterthan(v), response="w"))

# equalto
for(v in dirvt) tst(sum(dirm == v,na.rm=TRUE), summary(dirnw ~ equalto(v), response="w"))
for(v in undvt) tst(sum(undm == v,na.rm=TRUE)/2, summary(undnw ~ equalto(v), response="w"))
for(v in bipvt) tst(sum(bipm == v,na.rm=TRUE), summary(bipnw ~ equalto(v), response="w"))

set.seed(123)
for(tol in unique(c(runif(1,0,2), dirvt, undvt, bipvt))){
  for(v in dirvt) tst(sum(abs(dirm-v)<=tol,na.rm=TRUE), summary(dirnw ~ equalto(v,tol), response="w"))
  for(v in undvt) tst(sum(abs(undm-v)<=tol,na.rm=TRUE)/2, summary(undnw ~ equalto(v,tol), response="w"))
  for(v in bipvt) tst(sum(abs(bipm-v)<=tol,na.rm=TRUE), summary(bipnw ~ equalto(v,tol), response="w"))
}

# ininterval
charospec <- function(o1, o2) paste0(if(o1)'('else'[',if(o2)')'else']')
for(o1 in c(FALSE, TRUE)){
  for(o2 in c(FALSE, TRUE)){
    for(lv in c(-Inf,dirvt, Inf))
      for(uv in c(-Inf,dirvt, Inf)){
        truth <- sum(((o1 & dirm>lv) | (!o1 & dirm>=lv)) & ((o2 & dirm<uv) | (!o2 & dirm<=uv)), na.rm=TRUE)
        tst(truth, summary(dirnw ~ ininterval(lv, uv, c(o1,o2)), response="w"))
        tst(truth, summary(dirnw ~ ininterval(lv, uv, charospec(o1,o2)), response="w"))
      }
    
    for(lv in c(-Inf,undvt, Inf))
      for(uv in c(-Inf,undvt, Inf)){
        truth <- sum(((o1 & undm>lv) | (!o1 & undm>=lv)) & ((o2 & undm<uv) | (!o2 & undm<=uv)), na.rm=TRUE)/2
        tst(truth, summary(undnw ~ ininterval(lv, uv, c(o1,o2)), response="w"))
        tst(truth, summary(undnw ~ ininterval(lv, uv, charospec(o1,o2)), response="w"))
      }
    
    for(lv in c(-Inf,bipvt, Inf))
      for(uv in c(-Inf,bipvt, Inf)){
        truth <- sum(((o1 & bipm>lv) | (!o1 & bipm>=lv)) & ((o2 & bipm<uv) | (!o2 & bipm<=uv)), na.rm=TRUE)
        tst(truth, summary(bipnw ~ ininterval(lv, uv, c(o1,o2)), response="w"))
        tst(truth, summary(bipnw ~ ininterval(lv, uv, charospec(o1,o2)), response="w"))
      }
  }
}

# TODO: nodecovar

# smallerthan
for(v in dirvt) tst(sum(dirm < v,na.rm=TRUE), summary(dirnw ~ smallerthan(v), response="w"))
for(v in undvt) tst(sum(undm < v,na.rm=TRUE)/2, summary(undnw ~ smallerthan(v), response="w"))
for(v in bipvt) tst(sum(bipm < v,na.rm=TRUE), summary(bipnw ~ smallerthan(v), response="w"))
 
# nodecov
tst(sum(q*(dirm+t(dirm)),na.rm=TRUE), summary(dirnw ~ nodecov("q"), response="w"))
tst(c(sum(q*(dirm+t(dirm)),na.rm=TRUE),sum(q^2*(dirm+t(dirm)),na.rm=TRUE)), summary(dirnw ~ nodecov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q*((dirm!=0)+t(dirm!=0)),na.rm=TRUE), summary(dirnw ~ nodecov(function(x) x %v% "q", form="nonzero"), response="w"))
tst(c(sum(q*((dirm!=0)+t(dirm!=0)),na.rm=TRUE),sum(q^2*((dirm!=0)+t(dirm!=0)),na.rm=TRUE)), summary(dirnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

tst(sum(q*undm,na.rm=TRUE), summary(undnw ~ nodecov("q"), response="w"))
tst(c(sum(q*undm,na.rm=TRUE),sum(q^2*undm,na.rm=TRUE)), summary(undnw ~ nodecov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q*(undm!=0),na.rm=TRUE), summary(undnw ~ nodecov(function(x) x %v% "q", form="nonzero"), response="w"))
tst(c(sum(q*(undm!=0),na.rm=TRUE),sum(q^2*(undm!=0),na.rm=TRUE)), summary(undnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

tst(sum(q1*bipm+t(q2*t(bipm)),na.rm=TRUE), summary(bipnw ~ nodecov("q"), response="w"))
tst(c(sum(q1*bipm+t(q2*t(bipm)),na.rm=TRUE),sum(q1^2*bipm+t(q2^2*t(bipm)),na.rm=TRUE)), summary(bipnw ~ nodecov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q1*(bipm!=0)+t(q2*t(bipm!=0)),na.rm=TRUE), summary(bipnw ~ nodecov(~q, form="nonzero"), response="w"))
tst(c(sum(q1*(bipm!=0)+t(q2*t(bipm!=0)),na.rm=TRUE),sum(q1^2*(bipm!=0)+t(q2^2*t(bipm!=0)),na.rm=TRUE)), summary(bipnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

# nodefactor
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(dirm+t(dirm)),na.rm=TRUE)), summary(dirnw ~ nodefactor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*((dirm!=0)+t(dirm!=0)),na.rm=TRUE)), summary(dirnw ~ nodefactor(function(x) x %v% "f", base=base, form="nonzero"), response="w"))

  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*undm,na.rm=TRUE)), summary(undnw ~ nodefactor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(undm!=0),na.rm=TRUE)), summary(undnw ~ nodefactor(~f, base=base, form="nonzero"), response="w"))

  tst(sapply(sort(unique(f))[keep], function(x) sum((f1==x)*bipm+t((f2==x)*t(bipm)),na.rm=TRUE)), summary(bipnw ~ nodefactor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum((f1==x)*(bipm!=0)+t((f2==x)*t(bipm!=0)),na.rm=TRUE)), summary(bipnw ~ nodefactor(function(x) x %v% "f", base=base, form="nonzero"), response="w"))
}

# nodeicov
tst(sum(q*t(dirm),na.rm=TRUE), summary(dirnw ~ nodeicov("q"), response="w"))
tst(c(sum(q*t(dirm),na.rm=TRUE),sum(q^2*t(dirm),na.rm=TRUE)), summary(dirnw ~ nodeicov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q*t(dirm!=0),na.rm=TRUE), summary(dirnw ~ nodeicov(~q, form="nonzero"), response="w"))
tst(c(sum(q*t(dirm!=0),na.rm=TRUE),sum(q^2*t(dirm!=0),na.rm=TRUE)), summary(dirnw ~ nodeicov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

# TODO: nodeicovar

# nodeifactor
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*t(dirm),na.rm=TRUE)), summary(dirnw ~ nodeifactor("f", levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*t(dirm!=0),na.rm=TRUE)), summary(dirnw ~ nodeifactor(~f, base=base, form="nonzero"), response="w"))
}

# TODO: nodeisqrtcovar

# nodematch
tst(sum(abs(outer(f,f,"=="))*dirm,na.rm=TRUE), summary(dirnw ~ nodematch("f"), response="w"))
tst(sum(abs(outer(f,f,"=="))*(dirm!=0),na.rm=TRUE), summary(dirnw ~ nodematch(~f, form="nonzero"), response="w"))

for(keep in list(1, 1:2, 1:3)){
  tst(sapply(sort(unique(f))[keep], function(x) sum(abs(outer(f==x,f==x,"&"))*dirm,na.rm=TRUE)), summary(dirnw ~ nodematch("f",diff=TRUE, levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum(abs(outer(f==x,f==x,"&"))*(dirm!=0),na.rm=TRUE)), summary(dirnw ~ nodematch(~f, diff=TRUE, keep=keep, form="nonzero"), response="w"))
}

# TODO: nodemix

# nodeocov
tst(sum(q*dirm,na.rm=TRUE), summary(dirnw ~ nodeocov("q"), response="w"))
tst(c(sum(q*dirm,na.rm=TRUE),sum(q^2*dirm,na.rm=TRUE)), summary(dirnw ~ nodeocov(~poly(q,2,raw=TRUE)), response="w"))
tst(sum(q*(dirm!=0),na.rm=TRUE), summary(dirnw ~ nodeocov(~q, form="nonzero"), response="w"))
tst(c(sum(q*(dirm!=0),na.rm=TRUE),sum(q^2*(dirm!=0),na.rm=TRUE)), summary(dirnw ~ nodeocov(~poly(q,2,raw=TRUE), form="nonzero"), response="w"))

# TODO: nodeocovar

# nodeofactor
for(base in list(0, 1, 2, 1:2, 3)){
  keep <- if(all(base==0)) 1:3 else (1:3)[-base]
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*dirm,na.rm=TRUE)), summary(dirnw ~ nodeofactor(~f, levels=keep), response="w"))
  tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(dirm!=0),na.rm=TRUE)), summary(dirnw ~ nodeofactor("f", base=base, form="nonzero"), response="w"))
}

# TODO: nodeosqrtcovar

# TODO: nodesqrtcovar

# receiver
for(base in list(0, 1, 2, 1:2, 3)){
  i <- seq_len(network.size(dirnw))
  keep <- if(all(base==0)) i else i[-base]
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*t(dirm),na.rm=TRUE)), summary(dirnw ~ receiver(nodes=keep), response="w"))
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*t(dirm!=0),na.rm=TRUE)), summary(dirnw ~ receiver(nodes=keep, form="nonzero"), response="w"))
}

# sender
for(base in list(0, 1, 2, 1:2, 3)){
  i <- seq_len(network.size(dirnw))
  keep <- if(all(base==0)) i else i[-base]
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*dirm,na.rm=TRUE)), summary(dirnw ~ sender(nodes=keep), response="w"))
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*(dirm!=0),na.rm=TRUE)), summary(dirnw ~ sender(nodes=keep, form="nonzero"), response="w"))
}

# sociality
for(base in list(0, 1, 2, 1:2, 3)){
  i <- seq_len(network.size(dirnw))
  keep <- if(all(base==0)) i else i[-base]
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*undm,na.rm=TRUE)), summary(undnw ~ sociality(nodes=keep), response="w"))
  tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*(undm!=0),na.rm=TRUE)), summary(undnw ~ sociality(nodes=keep, form="nonzero"), response="w"))
}


# sum
tst(sum(dirm,na.rm=TRUE), summary(dirnw ~ sum, response="w"))
tst(sum(undm,na.rm=TRUE)/2, summary(undnw ~ sum, response="w"))
tst(sum(bipm,na.rm=TRUE), summary(bipnw ~ sum, response="w"))
