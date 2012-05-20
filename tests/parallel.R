#  File ergm/tests/parallel.R.off
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
data(florentine)

do.test<-function(type){
  cat("Testing",type,".\n")

  gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
               control=control.ergm(parallel=2, parallel.type=type), verbose=TRUE)

  summary(gest)
  # TODO get MCMC for parallel diagnostics to work.
  cat("Finished",type,".\n")
}

# Note that neither Rmpi nor rpvm are in the suggests list
# for the ergm package.  (Also, rpvm is not currently maintained.)
if(require(Rmpi)) do.test("MPI") else cat("Skipping MPI.\n")
do.test("SOCK")

