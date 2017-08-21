#  File R/parallel.utils.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
# Set up a flag for whether we are in charge of a cluster.
ergm.cluster.started <- local({
  started <- FALSE
  function(new){
    if(!missing(new))
      started <<- new
    else
      started
  }
})

ergm.MCMC.packagenames <- local({
  packagenames <- c("ergm") # It has to include itself.
  function(new){
    if(!missing(new))
      packagenames <<- unique(c(packagenames,new))
    else
      packagenames
  }
})

myLibLoc <- function()
  sub('/ergm/Meta/package.rds','',attr(packageDescription("ergm"),"file"))

# Acquires a cluster of specified type.
ergm.getCluster <- function(control, verbose=FALSE){
  
  if(inherits(control$parallel,"cluster")){
    ergm.cluster.started(FALSE)
    if(verbose) message("Cluster passed by user.")
    cl <- control$parallel
  }else{
    
    #type <- if(is.null(control$parallel.type)) getClusterOption("type") else control$parallel.type
    type <- if(is.null(control$parallel.type)) "PSOCK" else control$parallel.type
    
    if(verbose) message("Using ",type,".")
    
    #   Start Cluster
    cl <- switch(type,
                 # The rpvm package is apparently not being maintained.
                 PVM={              
                   #                capture.output(require(rpvm, quietly=TRUE, warn.conflicts = FALSE))
                   #                PVM.running <- try(.PVM.config(), silent=TRUE)
                   #                if(inherits(PVM.running,"try-error")){
                   #                  hostfile <- paste(Sys.getenv("HOME"),"/.xpvm_hosts",sep="")
                   #                  .PVM.start.pvmd(hostfile)
                   #                  message("PVM not running. Attempting to start.")
                   #                }
                   ergm.cluster.started(TRUE)
                   makeCluster(control$parallel,type="PVM")
                   
                 },
                 MPI={
                   
                   # Remember that we are responsible for it.
                   ergm.cluster.started(TRUE)
                   makeCluster(control$parallel,type="MPI")
                 },
                 SOCK={
                   ergm.cluster.started(TRUE)
                   makeCluster(control$parallel,type="PSOCK")
                   
                 },
                 PSOCK={
                   ergm.cluster.started(TRUE)
                   makeCluster(control$parallel,type="PSOCK")
                   
                 }
    )
  }
  # Set RNG up. 
  clusterSetRNGStream(cl)
  
  # On the off chance that user wants to load extra packages which we don't know about already.
  ergm.MCMC.packagenames(control$MCMC.packagenames)

  for(pkg in ergm.MCMC.packagenames()){
    # Try loading from the same location as the master.
    attached <- unlist(clusterCall(cl, require,
                                   package=pkg,
                                   character.only=TRUE,
                                   lib.loc=myLibLoc()))
    # If something failed, warn and try loading from anywhere.
    if(!all(attached)){
      if(verbose) message("Failed to attach ergm on the slave nodes from the same location as the master node. Will try to load from anywhere in the library path.")
      attached <- clusterCall(cl, require,
                              package=pkg,
                              character.only=TRUE)      
      if(!all(attached)) stop("Failed to attach ergm on one or more slave nodes. Make sure it's installed on or accessible from all of them and is in the library path.")
    }
    
    if(control$parallel.version.check){
      slave.versions <- clusterCall(cl,packageVersion,pkg)
      master.version <- packageVersion(pkg)
      
      if(!all(sapply(slave.versions,identical,master.version)))
        stop("The version of ",pkg, " attached on one or more slave nodes is different from from that on the master node (this node). Make sure that the same version is installed on all nodes. If you are absolutely certain that this message is in error, override with the parallel.version.check=FALSE control parameter.")
    }
  }
  cl
}


# Shuts down clusters.
ergm.stopCluster <- function(object, ...){
  UseMethod("ergm.stopCluster")
}

# Only stop the MPI cluster if we were the ones who had started it.
ergm.stopCluster.MPIcluster <- function(object, ...){
  if(ergm.cluster.started()){
    ergm.cluster.started(FALSE)
    stopCluster(object)
  }
}

ergm.stopCluster.default <- function(object, ...){
  if(ergm.cluster.started()){
    ergm.cluster.started(FALSE)
    stopCluster(object)
  }
}



ergm.sample.tomcmc<-function(sample, params){
  if (inherits(params$parallel,"cluster")) 
    nclus <- nrow(summary(params$parallel))
  else 
    nclus <- params$parallel
  
  samplesize <- nrow(sample)
  if(nclus > 1){
    
    samplesize<-round(samplesize / nclus)
    
    sample<-sapply(seq_len(nclus),function(i) {
      # Let mcmc() figure out the "end" from dimensions.
      mcmc(sample[(samplesize*(i-1)+1):(samplesize*i), , drop=FALSE], start = params$MCMC.burnin, thin = params$MCMC.interval)
    }, simplify=FALSE)
    
    do.call(mcmc.list,sample)
    
  }else{
    # Let mcmc() figure out the "end" from dimensions.
    mcmc(sample, start = params$MCMC.burnin, thin = params$MCMC.interval)
  }
}
