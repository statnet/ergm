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

# Commented out because neither Rmpi nor rpvm is in the suggests list
# for the ergm package.  (Also, rpvm is not currently maintained.)
#if(require(Rmpi)) do.test("MPI") else cat("Skipping MPI.\n")
#try({if(require(rpvm)) do.test("PVM") else cat("Skipping PVM.\n")})
do.test("SOCK")

