#  File tests/testthat/test-term-valued.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################


tst <- function(truth, fmla){
  test <- summary(fmla, response="w")
  expect_equal(test, truth, ignore_attr=TRUE)
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

bippnw <- bipnw
bippnw %e% "w" <- abs(bippnw %e% "w")
bippm <- as.matrix(bippnw, a="w")

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

dirpnw <- dirnw
dirpnw %e% "w" <- abs(dirpnw %e% "w")
dirpm <- as.matrix(dirpnw, a="w")
diag(dirpm) <- NA


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

undpnw <- undnw
undpnw %e% "w" <- abs(undpnw %e% "w")
undpm <- as.matrix(undpnw, a="w")
diag(undpm) <- NA

test_that("absdiff", {
  tst(sum(abs(outer(q,q,"-"))*dirm,na.rm=TRUE), dirnw ~ absdiff("q"))
  tst(sum(abs(outer(q,q,"-")^2)*dirm,na.rm=TRUE), dirnw ~ absdiff(~q,pow=2))
  tst(sum(abs(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), dirnw ~ absdiff(function(x) x %v% "q", form="nonzero"))
  tst(sum(abs(outer(q,q,"-")^2)*(dirm!=0),na.rm=TRUE), dirnw ~ absdiff("q",pow=2, form="nonzero"))

  tst(sum(abs(outer(q,q,"-"))*dirm,na.rm=TRUE), dirnw ~ B(~absdiff("q"), form="sum"))
  tst(sum(abs(outer(q,q,"-")^2)*dirm,na.rm=TRUE), dirnw ~ B(~absdiff("q",pow=2), form="sum"))
  tst(sum(abs(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), dirnw ~ B(~absdiff("q"), form="nonzero"))
  tst(sum(abs(outer(q,q,"-")^2)*(dirm!=0),na.rm=TRUE), dirnw ~ B(~absdiff("q",pow=2), form="nonzero"))
  tst(sum(abs(outer(q,q,"-"))*(dirm>.5 & dirm<1),na.rm=TRUE), dirnw ~ B(~absdiff("q"), form=~ininterval(.5,1)))
  tst(sum(abs(outer(q,q,"-")^2)*(dirm>.5 & dirm<1),na.rm=TRUE), dirnw ~ B(~absdiff("q",pow=2), form=~ininterval(.5,1)))

  tst(sum(abs(outer(q,q,"-"))*undm,na.rm=TRUE)/2, undnw ~ absdiff("q"))
  tst(sum(abs(outer(q,q,"-")^2)*undm,na.rm=TRUE)/2, undnw ~ absdiff(~q,pow=2))
  tst(sum(abs(outer(q,q,"-"))*(undm!=0),na.rm=TRUE)/2, undnw ~ absdiff(function(x) x %v% "q", form="nonzero"))
  tst(sum(abs(outer(q,q,"-")^2)*(undm!=0),na.rm=TRUE)/2, undnw ~ absdiff("q",pow=2, form="nonzero"))

  tst(sum(abs(outer(q,q,"-"))*undm,na.rm=TRUE)/2, undnw ~ B(~absdiff("q"), form="sum"))
  tst(sum(abs(outer(q,q,"-")^2)*undm,na.rm=TRUE)/2, undnw ~ B(~absdiff("q",pow=2), form="sum"))
  tst(sum(abs(outer(q,q,"-"))*(undm!=0),na.rm=TRUE)/2, undnw ~ B(~absdiff("q"), form="nonzero"))
  tst(sum(abs(outer(q,q,"-")^2)*(undm!=0),na.rm=TRUE)/2, undnw ~ B(~absdiff("q",pow=2), form="nonzero"))
  tst(sum(abs(outer(q,q,"-"))*(undm>.5 & undm<1),na.rm=TRUE)/2, undnw ~ B(~absdiff("q"), form=~ininterval(.5,1)))
  tst(sum(abs(outer(q,q,"-")^2)*(undm>.5 & undm<1),na.rm=TRUE)/2, undnw ~ B(~absdiff("q",pow=2), form=~ininterval(.5,1)))

  tst(sum(abs(outer(q1,q2,"-"))*bipm,na.rm=TRUE), bipnw ~ absdiff("q"))
  tst(sum(abs(outer(q1,q2,"-")^2)*bipm,na.rm=TRUE), bipnw ~ absdiff(~q,pow=2))
  tst(sum(abs(outer(q1,q2,"-"))*(bipm!=0),na.rm=TRUE), bipnw ~ absdiff(function(x) x %v% "q", form="nonzero"))
  tst(sum(abs(outer(q1,q2,"-")^2)*(bipm!=0),na.rm=TRUE), bipnw ~ absdiff("q",pow=2, form="nonzero"))

  tst(sum(abs(outer(q1,q2,"-"))*bipm,na.rm=TRUE), bipnw ~ B(~absdiff("q"), form="sum"))
  tst(sum(abs(outer(q1,q2,"-")^2)*bipm,na.rm=TRUE), bipnw ~ B(~absdiff("q",pow=2), form="sum"))
  tst(sum(abs(outer(q1,q2,"-"))*(bipm!=0),na.rm=TRUE), bipnw ~ B(~absdiff("q"), form="nonzero"))
  tst(sum(abs(outer(q1,q2,"-")^2)*(bipm!=0),na.rm=TRUE), bipnw ~ B(~absdiff("q",pow=2), form="nonzero"))
  tst(sum(abs(outer(q1,q2,"-"))*(bipm>.5 & bipm<1),na.rm=TRUE), bipnw ~ B(~absdiff("q"), form=~ininterval(.5,1)))
  tst(sum(abs(outer(q1,q2,"-")^2)*(bipm>.5 & bipm<1),na.rm=TRUE), bipnw ~ B(~absdiff("q",pow=2), form=~ininterval(.5,1)))
})

test_that("absdiffcat", {
  diffs <- sort(unique(c(abs(outer(q,q,"-")))))
  diffs <- diffs[diffs!=0]
  for(base in c(0, seq_along(diffs))){
    keep <- if(all(base==0)) seq_along(diffs) else seq_along(diffs)[-base]
    tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*dirm,na.rm=TRUE)), dirnw ~ absdiffcat("q",levels=keep))
    tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*(dirm!=0),na.rm=TRUE)), dirnw ~ absdiffcat(~q,base=base, form="nonzero"))

    tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*undm,na.rm=TRUE))/2, undnw ~ absdiffcat("q",levels=keep))
    tst(sapply(diffs[keep], function(x) sum((abs(outer(q,q,"-"))==x)*(undm!=0),na.rm=TRUE))/2, undnw ~ absdiffcat(function(x) x %v% "q",levels=keep, form="nonzero"))

    tst(sapply(diffs[keep], function(x) sum((abs(outer(q1,q2,"-"))==x)*bipm,na.rm=TRUE)), bipnw ~ absdiffcat("q", levels=keep))
    tst(sapply(diffs[keep], function(x) sum((abs(outer(q1,q2,"-"))==x)*(bipm!=0),na.rm=TRUE)), bipnw ~ absdiffcat(~q, base=base, form="nonzero"))
  }
})

test_that("atleast", {
  for(v in dirvt) tst(sum(dirm >= v,na.rm=TRUE), dirnw ~ atleast(v))
  tst(sapply(dirvt, function(v) sum(dirm >= v,na.rm=TRUE)), dirnw ~ atleast(dirvt))
  for(v in undvt) tst(sum(undm >= v,na.rm=TRUE)/2, undnw ~ atleast(v))
  tst(sapply(undvt, function(v) sum(undm >= v,na.rm=TRUE)/2), undnw ~ atleast(undvt))
  for(v in bipvt) tst(sum(bipm >= v,na.rm=TRUE), bipnw ~ atleast(v))
  tst(sapply(bipvt, function(v) sum(bipm >= v,na.rm=TRUE)), bipnw ~ atleast(bipvt))
})

test_that("atmost", {
  for(v in dirvt) tst(sum(dirm <= v,na.rm=TRUE), dirnw ~ atmost(v))
  tst(sapply(dirvt, function(v) sum(dirm <= v,na.rm=TRUE)), dirnw ~ atmost(dirvt))
  for(v in undvt) tst(sum(undm <= v,na.rm=TRUE)/2, undnw ~ atmost(v))
  tst(sapply(undvt, function(v) sum(undm <= v,na.rm=TRUE)/2), undnw ~ atmost(undvt))
  for(v in bipvt) tst(sum(bipm <= v,na.rm=TRUE), bipnw ~ atmost(v))
  tst(sapply(bipvt, function(v) sum(bipm <= v,na.rm=TRUE)), bipnw ~ atmost(bipvt))
})

test_that("b1cov", {
  tst(sum(q1*bipm,na.rm=TRUE), bipnw ~ b1cov("q"))
  tst(c(sum(q1*bipm,na.rm=TRUE),sum(q1^2*bipm,na.rm=TRUE)), bipnw ~ b1cov(~poly(q,2,raw=TRUE)))
  tst(sum(q1*(bipm!=0),na.rm=TRUE), bipnw ~ b1cov(~q, form="nonzero"))
  tst(c(sum(q1*(bipm!=0),na.rm=TRUE),sum(q1^2*(bipm!=0),na.rm=TRUE)), bipnw ~ b1cov(~poly(q,2,raw=TRUE), form="nonzero"))
})

test_that("b1factor", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(sapply(sort(unique(f1))[keep], function(x) sum((f1==x)*bipm,na.rm=TRUE)), bipnw ~ b1factor("f", levels=keep))
    tst(sapply(sort(unique(f1))[keep], function(x) sum((f1==x)*(bipm!=0),na.rm=TRUE)), bipnw ~ b1factor(~f, base=base, form="nonzero"))
  }
})

test_that("b1sociality", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(apply(bipm, 1, sum)[keep], bipnw ~ b1sociality(nodes=keep))
    tst(apply(bipm!=0, 1, sum)[keep], bipnw ~ b1sociality(nodes=keep, form="nonzero"))
  }
})

test_that("b2cov", {
  tst(sum(q2*t(bipm),na.rm=TRUE), bipnw ~ b2cov("q"))
  tst(c(sum(q2*t(bipm),na.rm=TRUE),sum(q2^2*t(bipm),na.rm=TRUE)), bipnw ~ b2cov(~poly(q,2,raw=TRUE)))
  tst(sum(q2*t(bipm!=0),na.rm=TRUE), bipnw ~ b2cov(function(x) x %v% "q", form="nonzero"))
  tst(c(sum(q2*t(bipm!=0),na.rm=TRUE),sum(q2^2*t(bipm!=0),na.rm=TRUE)), bipnw ~ b2cov(~poly(q,2,raw=TRUE), form="nonzero"))
})

test_that("b2factor", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(sapply(sort(unique(f2))[keep], function(x) sum((f2==x)*t(bipm),na.rm=TRUE)), bipnw ~ b2factor("f", levels=keep))
    tst(sapply(sort(unique(f2))[keep], function(x) sum((f2==x)*t(bipm!=0),na.rm=TRUE)), bipnw ~ b2factor(~f, base=base, form="nonzero"))
  }
})

test_that("b2sociality", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(apply(bipm, 2, sum)[keep], bipnw ~ b2sociality(nodes=keep))
    tst(apply(bipm!=0, 2, sum)[keep], bipnw ~ b2sociality(nodes=keep, form="nonzero"))
  }
})

test_that("edgecov", {
  tst(sum(dire*dirm,na.rm=TRUE), dirnw ~ edgecov("e"))
  tst(sum(dire*(dirm!=0),na.rm=TRUE), dirnw ~ edgecov("e", form="nonzero"))

  tst(sum(unde*undm,na.rm=TRUE)/2, undnw ~ edgecov("e"))
  tst(sum(unde*(undm!=0),na.rm=TRUE)/2, undnw ~ edgecov("e", form="nonzero"))

  tst(sum(bipe*bipm,na.rm=TRUE), bipnw ~ edgecov("e"))
  tst(sum(bipe*(bipm!=0),na.rm=TRUE), bipnw ~ edgecov("e", form="nonzero"))
})

test_that("edges", {
  tst(sum((dirm!=0),na.rm=TRUE), dirnw ~ edges)
  tst(sum((undm!=0),na.rm=TRUE)/2, undnw ~ edges)
  tst(sum((bipm!=0),na.rm=TRUE), bipnw ~ edges)
})

test_that("nonzero", {
  tst(sum((dirm!=0),na.rm=TRUE), dirnw ~ nonzero)
  tst(sum((undm!=0),na.rm=TRUE)/2, undnw ~ nonzero)
  tst(sum((bipm!=0),na.rm=TRUE), bipnw ~ nonzero)
})

test_that("diff", {
  posonly <- function(x) pmax(x, 0)
  negonly <- function(x) pmin(x, 0)

  for(dd in c("t-h", "h-t")){
    for(sa in c("identity", "abs", "posonly", "negonly")){
      saf <- get(sa)
      ddf <- switch(dd, `t-h`=identity, `h-t`=function(x) -x)

      df <- function(x) saf(ddf(x))

      tst(sum(df(outer(q,q,"-"))*dirm,na.rm=TRUE), dirnw ~ diff("q", dir=dd, sign.action=sa))
      tst(sum(df(outer(q,q,"-"))^2*dirm,na.rm=TRUE), dirnw ~ diff(~q,pow=2, dir=dd, sign.action=sa))
      tst(sum(df(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), dirnw ~ diff(function(x) x %v% "q", dir=dd, sign.action=sa, form="nonzero"))
      tst(sum(df(outer(q,q,"-"))^2*(dirm!=0),na.rm=TRUE), dirnw ~ diff("q",pow=2, dir=dd, sign.action=sa, form="nonzero"))

      tst(sum(df(outer(q,q,"-"))*dirm,na.rm=TRUE), dirnw ~ B(~diff("q", dir=dd, sign.action=sa), form="sum"))
      tst(sum(df(outer(q,q,"-"))^2*dirm,na.rm=TRUE), dirnw ~ B(~diff("q",pow=2, dir=dd, sign.action=sa), form="sum"))
      tst(sum(df(outer(q,q,"-"))*(dirm!=0),na.rm=TRUE), dirnw ~ B(~diff("q", dir=dd, sign.action=sa), form="nonzero"))
      tst(sum(df(outer(q,q,"-"))^2*(dirm!=0),na.rm=TRUE), dirnw ~ B(~diff("q",pow=2, dir=dd, sign.action=sa), form="nonzero"))

    }
  }
})

test_that("greaterthan", {
  for(v in dirvt) tst(sum(dirm > v,na.rm=TRUE), dirnw ~ greaterthan(v))
  tst(sapply(dirvt, function(v) sum(dirm > v,na.rm=TRUE)), dirnw ~ greaterthan(dirvt))
  for(v in undvt) tst(sum(undm > v,na.rm=TRUE)/2, undnw ~ greaterthan(v))
  tst(sapply(undvt, function(v) sum(undm > v,na.rm=TRUE)/2), undnw ~ greaterthan(undvt))
  for(v in bipvt) tst(sum(bipm > v,na.rm=TRUE), bipnw ~ greaterthan(v))
  tst(sapply(bipvt, function(v) sum(bipm > v,na.rm=TRUE)), bipnw ~ greaterthan(bipvt))
})

test_that("equalto", {
  for(v in dirvt) tst(sum(dirm == v,na.rm=TRUE), dirnw ~ equalto(v))
  for(v in undvt) tst(sum(undm == v,na.rm=TRUE)/2, undnw ~ equalto(v))
  for(v in bipvt) tst(sum(bipm == v,na.rm=TRUE), bipnw ~ equalto(v))

  set.seed(123)
  for(tol in unique(c(runif(1,0,2), dirvt, undvt, bipvt))){
    for(v in dirvt) tst(sum(abs(dirm-v)<=tol,na.rm=TRUE), dirnw ~ equalto(v,tol))
    for(v in undvt) tst(sum(abs(undm-v)<=tol,na.rm=TRUE)/2, undnw ~ equalto(v,tol))
    for(v in bipvt) tst(sum(abs(bipm-v)<=tol,na.rm=TRUE), bipnw ~ equalto(v,tol))
  }
})

test_that("ininterval", {
  charospec <- function(o1, o2) paste0(if(o1)'('else'[',if(o2)')'else']')
  for(o1 in c(FALSE, TRUE)){
    for(o2 in c(FALSE, TRUE)){
      for(lv in c(-Inf,dirvt, Inf))
        for(uv in c(-Inf,dirvt, Inf)){
          truth <- sum(((o1 & dirm>lv) | (!o1 & dirm>=lv)) & ((o2 & dirm<uv) | (!o2 & dirm<=uv)), na.rm=TRUE)
          tst(truth, dirnw ~ ininterval(lv, uv, c(o1,o2)))
          tst(truth, dirnw ~ ininterval(lv, uv, charospec(o1,o2)))
        }

      for(lv in c(-Inf,undvt, Inf))
        for(uv in c(-Inf,undvt, Inf)){
          truth <- sum(((o1 & undm>lv) | (!o1 & undm>=lv)) & ((o2 & undm<uv) | (!o2 & undm<=uv)), na.rm=TRUE)/2
          tst(truth, undnw ~ ininterval(lv, uv, c(o1,o2)))
          tst(truth, undnw ~ ininterval(lv, uv, charospec(o1,o2)))
        }

      for(lv in c(-Inf,bipvt, Inf))
        for(uv in c(-Inf,bipvt, Inf)){
          truth <- sum(((o1 & bipm>lv) | (!o1 & bipm>=lv)) & ((o2 & bipm<uv) | (!o2 & bipm<=uv)), na.rm=TRUE)
          tst(truth, bipnw ~ ininterval(lv, uv, c(o1,o2)))
          tst(truth, bipnw ~ ininterval(lv, uv, charospec(o1,o2)))
        }
    }
  }
})

test_that("nodecovar", {
  tst(sum(apply(undm, 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/2/(length(na.omit(r))-1))), undnw ~ nodecovar(FALSE))
  tst(sum(apply(sqrt(undpm), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/2/(length(na.omit(r))-1))), undpnw ~ nodecovar(FALSE, "sqrt"))
  tst(sum(apply(undm-mean(na.omit(c(undm))), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/2/(length(na.omit(r))-1))), undnw ~ nodecovar(TRUE))
  tst(sum(apply(sqrt(undpm)-mean(na.omit(c(sqrt(undpm)))), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/2/(length(na.omit(r))-1))), undpnw ~ nodecovar(TRUE, "sqrt"))
})

test_that("smallerthan", {
  for(v in dirvt) tst(sum(dirm < v,na.rm=TRUE), dirnw ~ smallerthan(v))
  tst(sapply(dirvt, function(v) sum(dirm < v,na.rm=TRUE)), dirnw ~ smallerthan(dirvt))
  for(v in undvt) tst(sum(undm < v,na.rm=TRUE)/2, undnw ~ smallerthan(v))
  tst(sapply(undvt, function(v) sum(undm < v,na.rm=TRUE)/2), undnw ~ smallerthan(undvt))
  for(v in bipvt) tst(sum(bipm < v,na.rm=TRUE), bipnw ~ smallerthan(v))
  tst(sapply(bipvt, function(v) sum(bipm < v,na.rm=TRUE)), bipnw ~ smallerthan(bipvt))
})

test_that("nodecov", {
  tst(sum(q*(dirm+t(dirm)),na.rm=TRUE), dirnw ~ nodecov("q"))
  tst(c(sum(q*(dirm+t(dirm)),na.rm=TRUE),sum(q^2*(dirm+t(dirm)),na.rm=TRUE)), dirnw ~ nodecov(~poly(q,2,raw=TRUE)))
  tst(sum(q*((dirm!=0)+t(dirm!=0)),na.rm=TRUE), dirnw ~ nodecov(function(x) x %v% "q", form="nonzero"))
  tst(c(sum(q*((dirm!=0)+t(dirm!=0)),na.rm=TRUE),sum(q^2*((dirm!=0)+t(dirm!=0)),na.rm=TRUE)), dirnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"))

  tst(sum(q*undm,na.rm=TRUE), undnw ~ nodecov("q"))
  tst(c(sum(q*undm,na.rm=TRUE),sum(q^2*undm,na.rm=TRUE)), undnw ~ nodecov(~poly(q,2,raw=TRUE)))
  tst(sum(q*(undm!=0),na.rm=TRUE), undnw ~ nodecov(function(x) x %v% "q", form="nonzero"))
  tst(c(sum(q*(undm!=0),na.rm=TRUE),sum(q^2*(undm!=0),na.rm=TRUE)), undnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"))

  tst(sum(q1*bipm+t(q2*t(bipm)),na.rm=TRUE), bipnw ~ nodecov("q"))
  tst(c(sum(q1*bipm+t(q2*t(bipm)),na.rm=TRUE),sum(q1^2*bipm+t(q2^2*t(bipm)),na.rm=TRUE)), bipnw ~ nodecov(~poly(q,2,raw=TRUE)))
  tst(sum(q1*(bipm!=0)+t(q2*t(bipm!=0)),na.rm=TRUE), bipnw ~ nodecov(~q, form="nonzero"))
  tst(c(sum(q1*(bipm!=0)+t(q2*t(bipm!=0)),na.rm=TRUE),sum(q1^2*(bipm!=0)+t(q2^2*t(bipm!=0)),na.rm=TRUE)), bipnw ~ nodecov(~poly(q,2,raw=TRUE), form="nonzero"))
})

test_that("nodefactor", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(dirm+t(dirm)),na.rm=TRUE)), dirnw ~ nodefactor("f", levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*((dirm!=0)+t(dirm!=0)),na.rm=TRUE)), dirnw ~ nodefactor(function(x) x %v% "f", base=base, form="nonzero"))

    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*undm,na.rm=TRUE)), undnw ~ nodefactor("f", levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(undm!=0),na.rm=TRUE)), undnw ~ nodefactor(~f, base=base, form="nonzero"))

    tst(sapply(sort(unique(f))[keep], function(x) sum((f1==x)*bipm+t((f2==x)*t(bipm)),na.rm=TRUE)), bipnw ~ nodefactor("f", levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum((f1==x)*(bipm!=0)+t((f2==x)*t(bipm!=0)),na.rm=TRUE)), bipnw ~ nodefactor(function(x) x %v% "f", base=base, form="nonzero"))
  }
})

test_that("nodeicov", {
  tst(sum(q*t(dirm),na.rm=TRUE), dirnw ~ nodeicov("q"))
  tst(c(sum(q*t(dirm),na.rm=TRUE),sum(q^2*t(dirm),na.rm=TRUE)), dirnw ~ nodeicov(~poly(q,2,raw=TRUE)))
  tst(sum(q*t(dirm!=0),na.rm=TRUE), dirnw ~ nodeicov(~q, form="nonzero"))
  tst(c(sum(q*t(dirm!=0),na.rm=TRUE),sum(q^2*t(dirm!=0),na.rm=TRUE)), dirnw ~ nodeicov(~poly(q,2,raw=TRUE), form="nonzero"))
})

test_that("nodeicovar", {
  tst(sum(apply(dirm, 2, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirnw ~ nodeicovar(FALSE))
  tst(sum(apply(sqrt(dirpm), 2, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirpnw ~ nodeicovar(FALSE, "sqrt"))
  tst(sum(apply(dirm-mean(na.omit(c(dirm))), 2, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirnw ~ nodeicovar(TRUE))
  tst(sum(apply(sqrt(dirpm)-mean(na.omit(c(sqrt(dirpm)))), 2, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirpnw ~ nodeicovar(TRUE, "sqrt"))
})

test_that("nodeifactor", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*t(dirm),na.rm=TRUE)), dirnw ~ nodeifactor("f", levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*t(dirm!=0),na.rm=TRUE)), dirnw ~ nodeifactor(~f, base=base, form="nonzero"))
  }
})

test_that("nodematch", {
  tst(sum(abs(outer(f,f,"=="))*dirm,na.rm=TRUE), dirnw ~ nodematch("f"))
  tst(sum(abs(outer(f,f,"=="))*(dirm!=0),na.rm=TRUE), dirnw ~ nodematch(~f, form="nonzero"))

  for(keep in list(1, 1:2, 1:3)){
    tst(sapply(sort(unique(f))[keep], function(x) sum(abs(outer(f==x,f==x,"&"))*dirm,na.rm=TRUE)), dirnw ~ nodematch("f",diff=TRUE, levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum(abs(outer(f==x,f==x,"&"))*(dirm!=0),na.rm=TRUE)), dirnw ~ nodematch(~f, diff=TRUE, keep=keep, form="nonzero"))
  }
})

# TODO: nodemix

test_that("nodeocov", {
  tst(sum(q*dirm,na.rm=TRUE), dirnw ~ nodeocov("q"))
  tst(c(sum(q*dirm,na.rm=TRUE),sum(q^2*dirm,na.rm=TRUE)), dirnw ~ nodeocov(~poly(q,2,raw=TRUE)))
  tst(sum(q*(dirm!=0),na.rm=TRUE), dirnw ~ nodeocov(~q, form="nonzero"))
  tst(c(sum(q*(dirm!=0),na.rm=TRUE),sum(q^2*(dirm!=0),na.rm=TRUE)), dirnw ~ nodeocov(~poly(q,2,raw=TRUE), form="nonzero"))
})

test_that("nodeocovar", {
  tst(sum(apply(dirm, 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirnw ~ nodeocovar(FALSE))
  tst(sum(apply(sqrt(dirpm), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirpnw ~ nodeocovar(FALSE, "sqrt"))
  tst(sum(apply(dirm-mean(na.omit(c(dirm))), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirnw ~ nodeocovar(TRUE))
  tst(sum(apply(sqrt(dirpm)-mean(na.omit(c(sqrt(dirpm)))), 1, function(r) (sum(na.omit(r)%o%na.omit(r))-sum(na.omit(r)^2))/(length(na.omit(r))-1))), dirpnw ~ nodeocovar(TRUE, "sqrt"))
})

test_that("nodeofactor", {
  for(base in list(0, 1, 2, 1:2, 3)){
    keep <- if(all(base==0)) 1:3 else (1:3)[-base]
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*dirm,na.rm=TRUE)), dirnw ~ nodeofactor(~f, levels=keep))
    tst(sapply(sort(unique(f))[keep], function(x) sum((f==x)*(dirm!=0),na.rm=TRUE)), dirnw ~ nodeofactor("f", base=base, form="nonzero"))
  }
})

# TODO: nodeosqrtcovar

# TODO: nodesqrtcovar

test_that("receiver", {
  for(base in list(0, 1, 2, 1:2, 3)){
    i <- seq_len(network.size(dirnw))
    keep <- if(all(base==0)) i else i[-base]
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*t(dirm),na.rm=TRUE)), dirnw ~ receiver(nodes=keep))
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*t(dirm!=0),na.rm=TRUE)), dirnw ~ receiver(nodes=keep, form="nonzero"))
  }
})

test_that("sender", {
  for(base in list(0, 1, 2, 1:2, 3)){
    i <- seq_len(network.size(dirnw))
    keep <- if(all(base==0)) i else i[-base]
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*dirm,na.rm=TRUE)), dirnw ~ sender(nodes=keep))
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*(dirm!=0),na.rm=TRUE)), dirnw ~ sender(nodes=keep, form="nonzero"))
  }
})

test_that("sociality", {
  for(base in list(0, 1, 2, 1:2, 3)){
    i <- seq_len(network.size(dirnw))
    keep <- if(all(base==0)) i else i[-base]
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*undm,na.rm=TRUE)), undnw ~ sociality(nodes=keep))
    tst(sapply(sort(unique(i))[keep], function(x) sum((i==x)*(undm!=0),na.rm=TRUE)), undnw ~ sociality(nodes=keep, form="nonzero"))
  }
})

test_that("sum", {
  tst(sum(dirm,na.rm=TRUE), dirnw ~ sum)
  tst(sum(undm,na.rm=TRUE)/2, undnw ~ sum)
  tst(sum(bipm,na.rm=TRUE), bipnw ~ sum)
})
