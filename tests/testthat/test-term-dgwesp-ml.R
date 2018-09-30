library(testthat)
context("test-term-dgwesp-ml.R")

# "Correct" transitivity calculators

dediag <- function(m, x=0) {diag(m) <- x; m}

UTP <- function(m1, m2){
  dediag(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2))
}

OTP <- function(m1, m2, in.order=TRUE){
  if(in.order) dediag(m1%*%m2)
  else dediag(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2))
}

ITP <- function(m1, m2, in.order=TRUE){
  if(in.order) dediag(t(m1%*%m2))
  else dediag(t(m1%*%m2+m2%*%m1-(m1*m2)%*%(m1*m2)))
}

OSP <- function(m1, m2){
  dediag(m1%*%t(m2)+m2%*%t(m1)-(m1*m2)%*%t(m1*m2))
}

ISP <- function(m1, m2){
  dediag(t(m1)%*%m2+t(m2)%*%m1-t(m1*m2)%*%(m1*m2))
}

esp <- function(x, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  L.base[L.base==0] <- NA # I.e., those with base=0 don't count at all.
  esp <- dediag(L.base*TP, NA)[upper.tri(TP)]
  tabulate(match(esp, x),length(x))
}

dsp <- function(x, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  dsp <- dediag(TP, NA)[upper.tri(TP)]
  tabulate(match(dsp, x),length(x))
}

nsp <- function(x, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- UTP(Ls.path1, Ls.path2, ...)
  L.base[L.base==1] <- NA # I.e., those with base=1 don't count at all.
  nsp <- dediag((1-L.base)*TP, NA)[upper.tri(TP)]
  tabulate(match(nsp, x),length(x))
}


desp <- function(x, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  L.base[L.base==0] <- NA # I.e., those with base=0 don't count at all.
  esp <- dediag(L.base*TP, NA)
  tabulate(match(esp, x),length(x))
}

ddsp <- function(x, type, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  dsp <- dediag(TP, NA)
  tabulate(match(dsp, x),length(x))
}

dnsp <- function(x, type, L.base, Ls.path1, Ls.path2=Ls.path1, ...){
  TP <- type(Ls.path1, Ls.path2, ...)
  L.base[L.base==1] <- NA # I.e., those with base=1 don't count at all.
  nsp <- dediag((1-L.base)*TP, NA)
  tabulate(match(nsp, x),length(x))
}

library(ergm)
library(purrr)
n <- 5
nw1 <- nw2 <- nw3 <- network.initialize(n,dir=TRUE)
lnw <- Layer(nw1,nw2,nw3)

#### Some code useful for debugging.

## # Construct a transitive triad base-first and then remove the base:
## testseq1 <- list(matrix(c(1,2),ncol=2,byrow=TRUE), # Base
##                  matrix(c(1+n,3+n),ncol=2,byrow=TRUE), # Segment 1
##                  matrix(c(3+2*n,2+2*n),ncol=2,byrow=TRUE), # Segment 2
##                  matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base

## # Construct a transitive triad base-last and then remove the base:
## testseq2 <- list(matrix(c(1+n,3+n),ncol=2,byrow=TRUE), # Segment 1
##                  matrix(c(3+2*n,2+2*n),ncol=2,byrow=TRUE), # Segment 2
##                  matrix(c(1,2),ncol=2,byrow=TRUE), # Base
##                  matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base


## testseq3 <- list(
##   matrix(c(1,2),ncol=2,byrow=TRUE), # Base in L1
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE), # Segment 1 in L3
##   matrix(c(3+1*n,2+1*n),ncol=2,byrow=TRUE), # Segment 2 in L2
##   matrix(c(1,2),ncol=2,byrow=TRUE)) # -Base in L1

## # Construct a transitive triad base-last and then remove the base:
## testseq4 <- list(
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE), # Segment 1 in L3
##   matrix(c(3+1*n,2+1*n),ncol=2,byrow=TRUE), # Segment 2 in L2
##   matrix(c(1,2),ncol=2,byrow=TRUE), # Base in L1
##   matrix(c(1+2*n,3+2*n),ncol=2,byrow=TRUE) # -Segment 1 in L3
##   )

## ergm.godfather(lnw~edges+
##                  ## desp(1,type="OTP",L.base=~`1`,Ls.path=c(~`3`,~`2`),L.in_order=TRUE)+
##                  ## desp(1,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
##                  desp(1,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE),
##                changes = testseq4,stats.start=TRUE)

## summary(lnw~desp(1:18,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE))
## ergm.godfather(lnw~edges+desp(1:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE),
##                changes = testseq1,stats.start=TRUE)
## ergm.godfather(lnw~edges+desp(1:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE),
##                changes = testseq2,stats.start=TRUE)

ctrl <- control.simulate.formula(MCMC.burnin=1, MCMC.interval=1)

for(cache.sp in c(FALSE,TRUE)){
  options(ergm.term=list(cache.sp=cache.sp))
  sptxt <- if(cache.sp) "with shared partner caching" else "without shared partner caching"

test_that(paste("Multilayer dgw*sp statistics for homogeneously directed networks",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # desp distinct layers
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct layers
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct layers
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # desp base and path distinct
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path distinct
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path distinct
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # desp base and path same
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path same
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path same
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # desp distinct base and one layer
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct base and one layer
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct base and one layer
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # desp distinct layers
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddsp distinct layers
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnsp distinct layers
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # desp base and path distinct
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # ddsp base and path distinct
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnsp base and path distinct
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # desp base and path same
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # ddsp base and path same
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnsp base and path same
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # desp distinct base and one layer
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # ddsp distinct base and one layer
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnsp distinct base and one layer
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3)

                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})

### Heterogeneous directedness.
nw1 <- nw3 <- network.initialize(n,dir=TRUE)
nw2 <- network.initialize(n,dir=FALSE)
lnw <- Layer(nw1,nw2,nw3)

test_that(paste("Multilayer dgw*sp statistics for heterogeneously directed networks 1",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # desp distinct layers
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct layers
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct layers
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # desp base and path distinct
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    desp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path distinct
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path distinct
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`1`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    dnsp(0:n,type="ISP",L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # desp base and path same
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path same
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`2`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path same
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`2`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # desp distinct base and one layer
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    desp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    desp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    desp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct base and one layer
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    ddsp(0:n,type="OTP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="ITP",Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    ddsp(0:n,type="OSP",Ls.path=c(~`2`,~`3`))+
                    ddsp(0:n,type="ISP",Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct base and one layer
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=TRUE)+
                    dnsp(0:n,type="OTP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="ITP",L.base=~`2`,Ls.path=c(~`2`,~`3`),L.in_order=FALSE)+
                    dnsp(0:n,type="OSP",L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    dnsp(0:n,type="ISP",L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m2 <- m2 + t(m2)
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # desp distinct layers
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # ddsp distinct layers
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnsp distinct layers
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # desp base and path distinct
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      desp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # ddsp base and path distinct
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnsp base and path distinct
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m1,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      dnsp(0:n,ISP,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # desp base and path same
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # ddsp base and path same
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m2),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m2),
                      # dnsp base and path same
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m2, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # desp distinct base and one layer
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      desp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      desp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      desp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # ddsp distinct base and one layer
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      ddsp(0:n,OTP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,ITP,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      ddsp(0:n,OSP,Ls.path1=m2,Ls.path2=m3),
                      ddsp(0:n,ISP,Ls.path1=m2,Ls.path2=m3),
                      # dnsp distinct base and one layer
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=TRUE),
                      dnsp(0:n,OTP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,ITP,L.base=m2,Ls.path1=m2,Ls.path2=m3, in.order=FALSE),
                      dnsp(0:n,OSP,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      dnsp(0:n,ISP,L.base=m2,Ls.path1=m2,Ls.path2=m3)

                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})


### Undirected.
nw1 <- nw2 <- nw3 <- network.initialize(n,dir=FALSE)
lnw <- Layer(nw1,nw2,nw3)

test_that(paste("Multilayer dgw*sp statistics for undirected networks",sptxt), {
  sim <- suppressWarnings(simulate(lnw~
                    # desp distinct layers
                    desp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct layers
                    ddsp(0:n,Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct layers
                    dnsp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`3`))+
                    # desp base and path distinct
                    desp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path distinct
                    ddsp(0:n,Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path distinct
                    dnsp(0:n,L.base=~`1`,Ls.path=c(~`2`,~`2`))+
                    # desp base and path same
                    desp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # ddsp base and path same
                    ddsp(0:n,Ls.path=c(~`2`,~`2`))+
                    # dnsp base and path same
                    dnsp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`2`))+
                    # desp distinct base and one layer
                    desp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))+
                    # ddsp distinct base and one layer
                    ddsp(0:n,Ls.path=c(~`2`,~`3`))+
                    # dnsp distinct base and one layer
                    dnsp(0:n,L.base=~`2`,Ls.path=c(~`2`,~`3`))
                 ,
                  control=ctrl,
                  nsim=200))

  stats <- sapply(sim,
                  function(nw){
                    n <- network.size(nw)/3
                    m <- as.matrix(nw)
                    m1 <- m[seq_len(n),seq_len(n)]
                    m2 <- m[seq_len(n)+n,seq_len(n)+n]
                    m3 <- m[seq_len(n)+n*2,seq_len(n)+n*2]
                    
                    c(
                      # esp distinct layers
                      esp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # dsp distinct layers
                      dsp(0:n,Ls.path1=m2,Ls.path2=m3),
                      # nsp distinct layers
                      nsp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m3),
                      # esp base and path distinct
                      esp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # dsp base and path distinct
                      dsp(0:n,Ls.path1=m2,Ls.path2=m2),
                      # nsp base and path distinct
                      nsp(0:n,L.base=m1,Ls.path1=m2,Ls.path2=m2),
                      # esp base and path same
                      esp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # dsp base and path same
                      dsp(0:n,Ls.path1=m2,Ls.path2=m2),
                      # nsp base and path same
                      nsp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m2),
                      # esp distinct base and one layer
                      esp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3),
                      # dsp distinct base and one layer
                      dsp(0:n,Ls.path1=m2,Ls.path2=m3),
                      # nsp distinct base and one layer
                      nsp(0:n,L.base=m2,Ls.path1=m2,Ls.path2=m3)
                    )
                  }) %>% t()

  expect_equivalent(attr(sim,"stats"), stats)
})
}
