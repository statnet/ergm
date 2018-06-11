#  File R/parallel.utils.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
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

#' @importFrom utils packageDescription
myLibLoc <- function()
  sub('/ergm/Meta/package.rds','',attr(packageDescription("ergm"),"file"))

#' Parallel Processing in the \code{\link[=ergm-package]{ergm}} Package
#' 
#' For estimation that require MCMC, \code{\link[=ergm-package]{ergm}}
#' can take advantage of multiple CPUs or CPU cores on the system on
#' which it runs, as well as computing clusters. It uses package
#' \code{parallel} and \code{snow} to facilitate this, and supports
#' all cluster types that they does. The number of nodes used and the
#' parallel API are controlled using the \code{parallel} and
#' \code{parallel.type} arguments passed to the control functions,
#' such as \code{\link{control.ergm}}.
#' 
#' 
#' Further details on the various cluster types are included below.
#' 
#' 
#' @name ergm-parallel
#' @aliases ergm-parallel parallel ergm.parallel parallel.ergm parallel-ergm
#' 
#' @docType methods
#' @section PSOCK clusters: The \code{parallel} package is used with
#'   PSOCK clusters by default, to utilize multiple cores on a
#'   system. The number of cores on a system can be determined with
#'   the \code{detectCores} function.
#' 
#'   This method works with the base installation of R on all
#'   platforms, and does not require additional software.
#' 
#'   For more advanced applications, such as clusters that span
#'   multiple machines on a network, the clusters can be initialized
#'   manually, and passed into \code{ergm} using the \code{parallel}
#'   control argument. See the second example below.
#'
#' @section MPI clusters: To use MPI to accelerate ERGM sampling, pass
#'   the control parameter \code{parallel.type="MPI"}.
#'   \code{\link[=ergm-package]{ergm}} requires the \code{snow} and
#'   \code{Rmpi} packages to communicate with an MPI cluster.
#'  
#'   Using MPI clusters requires the system to have an existing MPI
#'   installation.  See the MPI documentation for your particular
#'   platform for instructions.
#'
#'   To use `ergm` across multiple machines in a high performance
#'   computing environment, see the section "User initiated clusters"
#'   below.
#'
#'
#' @section User initiated clusters: A cluster can be passed into
#'   \code{ergm} with the \code{parallel} control parameter.
#'   \code{ergm} will detect the number of nodes in the cluster, and
#'   use all of them for MCMC sampling. This method is flexible: it
#'   will accept any cluster type that is compatible with \code{snow}
#'   or \code{parallel} packages. Usage examples for a
#'   multiple-machine high performance MPI cluster can be found at the
#'   statnet wiki:
#'   \url{https://statnet.csde.washington.edu/trac/wiki/ergmParallel}
#'
#'
#' @examples
#' 
#' \donttest{
#' # Uses 2 SOCK clusters for MCMLE estimation
#' data(faux.mesa.high)
#' nw <- faux.mesa.high
#' fauxmodel.01 <- ergm(nw ~ edges + isolates + gwesp(0.2, fixed=TRUE), 
#'                      control=control.ergm(parallel=2, parallel.type="PSOCK"))
#' summary(fauxmodel.01)
#' 
#' }
#'
NULL

#' @rdname ergm-parallel
#' @description The \code{ergm.getCluster} function is usually called
#'   internally by the ergm process (in
#'   \code{\link{ergm_MCMC_sample}}) and will attempt to start the
#'   appropriate type of cluster indicated by the
#'   \code{\link{control.ergm}} settings. It will also check that the
#'   same version of `ergm` is installed on each node.
#' @param control a \code{\link{control.ergm}} (or similar) list of
#'   parameter values from which the parallel settings should be read.
#' @param verbose logical, should detailed status info be printed to
#'   console?
#' 
#' @export ergm.getCluster
ergm.getCluster <- function(control, verbose=FALSE){
  
  if(inherits(control$parallel,"cluster")){
    ergm.cluster.started(FALSE)
    if(verbose) message("Cluster passed by user.")
    cl <- control$parallel
  }else{
    
    #type <- if(is.null(control$parallel.type)) getClusterOption("type") else control$parallel.type
    type <- NVL(control$parallel.type, "PSOCK")
    
    if(verbose) message("Using ",type,".")
    
    #   Start Cluster
    #' @importFrom parallel makeCluster
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
  #' @importFrom parallel clusterSetRNGStream
  clusterSetRNGStream(cl)
  
  # On the off chance that user wants to load extra packages which we don't know about already.
  ergm.MCMC.packagenames(control$MCMC.packagenames)

  for(pkg in ergm.MCMC.packagenames()){
    # Try loading from the same location as the master.
  #' @importFrom parallel clusterCall
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
      #' @importFrom utils packageVersion
      master.version <- packageVersion(pkg)
      
      if(!all(sapply(slave.versions,identical,master.version)))
        stop("The version of ",pkg, " attached on one or more slave nodes is different from from that on the master node (this node). Make sure that the same version is installed on all nodes. If you are absolutely certain that this message is in error, override with the parallel.version.check=FALSE control parameter.")
    }
  }
  cl
}


#' @rdname ergm-parallel
#' @description The \code{ergm.stopCluster} shuts down a
#'   cluster, but only if `ergm.getCluster` was responsible for
#'   starting it.
#'
#' @param object an object, probably of class `cluster`.
#' @param \dots not currently used
#' @export ergm.stopCluster
ergm.stopCluster <- function(object, ...){
  UseMethod("ergm.stopCluster")
}

#' @rdname ergm-parallel
#' @importFrom parallel stopCluster
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
