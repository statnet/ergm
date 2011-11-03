# Acquires a cluster of specified type.
ergm.getCluster <- function(control, verbose=FALSE){
  capture.output(require(snow, quietly=TRUE, warn.conflicts = FALSE))

  type <- if(is.null(control$parallel.type)) getClusterOption("type") else control$parallel.type

  if(verbose) cat("Using ",type,".\n", sep="")
    #   Start Cluster

  cl <- switch(type,
               PVM={              
                 PVM.running <- try(.PVM.config(), silent=TRUE)
                 if(inherits(PVM.running,"try-error")){
                   hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
                   .PVM.start.pvmd(hostfile)
                   cat("PVM not running. Attempting to start.\n")
                 }
                 makeCluster(control$parallel,type="PVM")
               },
               MPI={
                 # See if a preexisting cluster exists.
                 if(is.null(getMPIcluster()))
                   makeCluster(control$parallel,type="MPI")
                 else
                   getMPIcluster()
               },
               SOCK={
                 makeCluster(control$parallel,type="SOCK")
               }
               )
  
  # Set things up. 
  clusterSetupRNG(cl)
  if("ergm" %in% control$packagenames){
    clusterEvalQ(cl,library(ergm))
  }
  cl
}


# Shuts down clusters. MPI cannot be stopped reliably and restarted
# without closing R.
ergm.stopCluster <- function(object, ...)
  UseMethod("ergm.stopCluster")

#ergm.stopCluster.MPIcluster <- function(object, ...){
#  # Don't stop MPI clusters!
#}

ergm.stopCluster.default <- function(object, ...){
  stopCluster(object)
}



ergm.sample.tomcmc<-function(sample, params){
  library(coda)
  if(params$parallel){

    samplesize<-round(params$samplesize / params$parallel)

    sample<-sapply(seq_len(params$parallel),function(i) {
      # Let mcmc() figure out the "end" from dimensions.
      mcmc(sample[(samplesize*(i-1)+1):(samplesize*i),], start = params$burnin, thin = params$interval)
    }, simplify=FALSE)

    do.call("mcmc.list",sample)

  }else{
    # Let mcmc() figure out the "end" from dimensions.
    mcmc(sample, start = params$burnin, thin = params$interval)
  }
}
