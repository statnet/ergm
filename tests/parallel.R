#  File tests/parallel.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(ergm)
data(florentine)

for(type in c("SOCK","MPI","MPI2")){
  cat("\n\n======= Testing",type,"=======\n\n")
  
  if(type=="MPI2"){
    cl <- makeCluster(2,"MPI")
    type <- "MPI"
  }
  
  gest <- ergm(flomarriage ~ edges + absdiff("wealth"),
               eval.loglik=FALSE,
               control=control.ergm(MCMC.burnin=1000, MCMC.interval=10, MCMLE.maxit=2, MCMC.samplesize=1000, force.main=TRUE,
                 parallel=2, parallel.type=type))

  print(summary(gest))
  mcmc.diagnostics(gest)

  # FIXME: Set seeds and replace with actual values?
  sim.STAT.SEQ <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), statsonly=TRUE, sequential=TRUE)
  stopifnot(nrow(sim.STAT.SEQ)==5)
  print(sim.STAT.SEQ)

  sim.STAT.seq <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), statsonly=TRUE, sequential=FALSE)
  stopifnot(nrow(sim.STAT.seq)==5)
  print(sim.STAT.seq)

  sim.stat.SEQ <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), statsonly=FALSE, sequential=TRUE)
  stopifnot(length(sim.stat.SEQ)==5)
  print(sim.stat.SEQ)

  sim.stat.seq <- simulate(gest, nsim=5, control=control.simulate.ergm(parallel=2, parallel.type=type), statsonly=FALSE, sequential=FALSE)
  stopifnot(length(sim.stat.seq)==5)
  print(sim.stat.seq)
  
  if(exists("cl")){
    stopCluster(cl)
    rm(cl)
  }
}
  
}, "parallel")
