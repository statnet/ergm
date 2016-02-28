# more comprehensive parallel testing

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
#do.test("SOCK")

# Updated tests

# SOCK
data(faux.magnolia.high)
nw <- faux.magnolia.high
t0 <- proc.time()
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=T), 
                     control=control.ergm(parallel=3, parallel.type="SOCK", 
                                          MCMLE.maxit=100))
proc.time() - t0

# default type
data(faux.magnolia.high)
nw <- faux.magnolia.high
t0 <- proc.time()
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=T), 
                     control=control.ergm(parallel=3, 
                                          MCMLE.maxit=100))
proc.time() - t0

# pre-made SOCK cluster
clus <- makeCluster(3, type='PSOCK')
clus
data(faux.magnolia.high)
nw <- faux.magnolia.high
t0 <- proc.time()
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=T), 
                     control=control.ergm(parallel=clus, 
                                          MCMLE.maxit=100))
proc.time() - t0
stopCluster(clus)

# only works on linux with Rmpi installed, and OpenMPI
# MPI
data(faux.magnolia.high)
nw <- faux.magnolia.high
t0 <- proc.time()
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=T), 
                     control=control.ergm(parallel=3, parallel.type="MPI", 
                                          MCMLE.maxit=100))
proc.time() - t0

# pre-made MPI cluster
clus <- makeCluster(3, type='MPI')
clus
data(faux.magnolia.high)
nw <- faux.magnolia.high
t0 <- proc.time()
fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=T), 
                     control=control.ergm(parallel=clus, 
                                          MCMLE.maxit=100))
proc.time() - t0