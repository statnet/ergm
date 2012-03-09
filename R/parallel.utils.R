#  File ergm/R/parallel.utils.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
# Set up a flag for whether we are in charge of MPI cluster.
ergm.MPIcluster.started <- FALSE

# Acquires a cluster of specified type.
ergm.getCluster <- function(control, verbose=FALSE){
  capture.output(require(snow, quietly=TRUE, warn.conflicts = FALSE))
# The rpvm package is apparently not being maintained.
#  capture.output(require(rpvm, quietly=TRUE, warn.conflicts = FALSE))

  type <- if(is.null(control$parallel.type)) getClusterOption("type") else control$parallel.type

  if(verbose) cat("Using ",type,".\n", sep="")
    #   Start Cluster

  cl <- switch(type,
# The rpvm package is apparently not being maintained.
#               PVM={              
#                 PVM.running <- try(.PVM.config(), silent=TRUE)
#                 if(inherits(PVM.running,"try-error")){
#                   hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
#                   .PVM.start.pvmd(hostfile)
#                   cat("PVM not running. Attempting to start.\n")
#                 }
#                 makeCluster(control$parallel,type="PVM")
#               },


# FIXME:  Excerpt from an email by Kurt Hornik:
#* checking R code for possible problems ... NOTE
#Found the following possibly unsafe calls:
#File 'R/parallel.utils.R':
# found unlockBinding("ergm.MPIcluster.started",
#   environment(ergm.getCluster))
#
#   ...
#   
#Can you pls fix?  Note that there are three issues: for the first, the
#simplest should be using a dynamic variable, e.g.
#
# ergm.MPIcluster.started <-
# local({
#     started <- FALSE
#     function(new) {
#         if(!missing(new))
#	      started <<- new
#	  else
#	      started
#     }
#     })
#
#and then use as
#
#    ergm.MPIcluster.started(TRUE)
#
#to set and
#
#    ergm.MPIcluster.started()
#
#to get
################ As a result of above problem, "MPI" bit is commented out:
#   
#   
#               MPI={
#                 # See if a preexisting cluster exists.
#                 if(is.null(getMPIcluster())){
#                   # Remember that we are responsible for it.
#                   unlockBinding("ergm.MPIcluster.started", environment(ergm.getCluster))                  
#                   assign("ergm.MPIcluster.started", TRUE, environment(ergm.getCluster))
#                   lockBinding("ergm.MPIcluster.started", environment(ergm.getCluster))
#                   makeCluster(control$parallel,type="MPI")
#                 }else
#                   getMPIcluster()
#               },
               SOCK={
                 makeCluster(control$parallel,type="SOCK")
               }
               )
  
  # Set things up. 
  clusterSetupRNG(cl)
  if("ergm" %in% control$MCMC.packagenames){

    # Try loading from the same location as the master.
    attached <- unlist(clusterCall(cl, require,
                                 package="ergm",
                                 character.only=TRUE))
## FIXME:  setting binding for myLibLoc in the .onLoad function
## is broken; did not have time to fix this.
#                                 ,lib.loc=myLibLoc))
    # If something failed, warn and try loading from anywhere.
    if(!all(attached)){
      if(verbose) cat("Failed to attach ergm on the slave nodes from the same location as the master node. Will try to load from anywhere in the library path.\n")
      attached <- clusterCall(cl, require,
                            package="ergm",
                            character.only=TRUE)      
      if(!all(attached)) stop("Failed to attach ergm on one or more slave nodes. Make sure it's installed on or accessible from all of them and is in the library path.")
    }
    
    if(control$parallel.version.check){
      slave.versions <- clusterCall(cl,packageVersion,"ergm")
      master.version <- packageVersion("ergm")

      if(!all(sapply(slave.versions,identical,master.version)))
        stop("The version of ergm attached on one or more slave nodes is different from from that on the master node (this node). Make sure the same version is installed on all nodes. If you are absolutely certain that this message is in error, override with the parallel.version.check=FALSE control parameter.")
    }
  }
  cl
}


# Shuts down clusters.
ergm.stopCluster <- function(object, ...)
  UseMethod("ergm.stopCluster")

# Only stop the MPI cluster if we were the ones who had started it.
ergm.stopCluster.MPIcluster <- function(object, ...){
  if(ergm.MPIcluster.started) stopCluster(object)
}

ergm.stopCluster.default <- function(object, ...){
  stopCluster(object)
}



ergm.sample.tomcmc<-function(sample, params){
  library(coda)
  samplesize <- nrow(sample)
  if(params$parallel){

    samplesize<-round(samplesize / params$parallel)

    sample<-sapply(seq_len(params$parallel),function(i) {
      # Let mcmc() figure out the "end" from dimensions.
      mcmc(sample[(samplesize*(i-1)+1):(samplesize*i),], start = params$MCMC.burnin, thin = params$MCMC.interval)
    }, simplify=FALSE)

    do.call("mcmc.list",sample)

  }else{
    # Let mcmc() figure out the "end" from dimensions.
    mcmc(sample, start = params$MCMC.burnin, thin = params$MCMC.interval)
  }
}
