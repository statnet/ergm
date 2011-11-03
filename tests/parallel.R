library(ergm)
data(florentine)

do.test<-function(type){
  cat("Testing",type,".\n")

  gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle,
               control=control.ergm(parallel=3, parallel.type=type), verbose=TRUE)

  summary(gest)
  # TODO get MCMC for parallel diagnostics to work.
  cat("Finished",type,".\n")
}

if(require(Rmpi)) do.test("MPI") else cat("Skipping MPI.\n")
try({if(require(rpvm)) do.test("PVM") else cat("Skipping PVM.\n")})
do.test("SOCK")

