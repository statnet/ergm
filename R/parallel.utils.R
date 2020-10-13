#  File R/parallel.utils.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
# Save the cluster we are in charge of.
ergm.cluster.started <- local({
  started <- NULL
  function(new){
    if(!missing(new))
      started <<- new
    else
      started
  }
})


#' Keep track of packages that slave nodes need to load.
#' @noRd
#'
#' @param pending if not missing and not [`character`], return the
#'   current list of pending packages; if [`character`], add to list
#'   of pending packages if not already loaded; return the result.
#' @param loaded if not missing and not [`character`], return the
#'   current list of loaded packages; if [`character`], add to list of
#'   loaded packages and remove from the list of pending packages;
#'   return the result.
#' @param reset if `TRUE`, move all loaded packages to pending; return
#'   the result.
ergm.MCMC.packagenames <- local({
  pending.packages <- c("ergm") # It has to include itself.
  loaded.packages <- c()
  function(pending, loaded, reset=FALSE){
    if(!missing(pending)){
      if(is.character(pending))
        pending.packages <<- setdiff(unique(c(pending.packages,pending)),loaded.packages)
      pending.packages
    }else if(!missing(loaded)){
      if(is.character(loaded)){
        loaded.packages <<- unique(c(loaded.packages,loaded))
        pending.packages <<- setdiff(pending.packages,loaded.packages)
      }
      loaded.packages
    }else if(reset){
      pending.packages <<- unique(c(pending.packages,loaded.packages))
      loaded.packages <<- c()
      pending.packages
    }
  }
})

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
#'   or \code{parallel} packages.
## #'   Usage examples for a
## #'   multiple-machine high performance MPI cluster can be found at the
## #'   statnet wiki:
## #'   \url{https://statnet.csde.washington.edu/trac/wiki/ergmParallel}
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
#'   parameter values from which the parallel settings should be read;
#'   can also be [`NULL`], in which case an existing cluster is used
#'   if started, or no cluster otherwise.
#' @param verbose logical, should detailed status info be printed to
#'   console?
#' @param stop_on_exit An [`environment`] or `NULL`. If an
#'   `environment`, defaulting to that of the calling function, the
#'   cluster will be stopped when the calling the frame in question
#'   exits.
#' 
#' @export ergm.getCluster
ergm.getCluster <- function(control=NULL, verbose=FALSE, stop_on_exit=parent.frame()){
  # If we don't want a cluster, just return NULL.
  if (is.numeric(control$parallel) && control$parallel==0) return(NULL)
  
  if(inherits(control$parallel,"cluster")){
    # Control argument *is* a cluster. Overrides everything.
    if(verbose) message("Cluster passed by user.")
    cl <- control$parallel
  }else if(!is.null(ergm.cluster.started())){
    if(verbose) message("Reusing the running cluster.")
    cl <- ergm.cluster.started()
  }else if(is.null(control)){
    return(NULL)
  }else{
    if(verbose) message("Starting a new cluster.")
    if (!is.numeric(control$parallel))
      warning("Unrecognized value passed to parallel= control parameter.")

    #type <- if(is.null(control$parallel.type)) getClusterOption("type") else control$parallel.type
    type <- NVL(control$parallel.type, "PSOCK")
    
    if(verbose) message("Using ",type,".")
    
    #   Start Cluster
    #' @importFrom parallel makeCluster
    cl <- makeCluster(control$parallel, type=type)
    ergm.cluster.started(cl)

    # Based on https://yihui.name/en/2017/12/on-exit-parent/ .
    if(!is.null(stop_on_exit)) do.call(on.exit, list(substitute(ergm.stopCluster(verbose=verbose)), add=TRUE),envir=stop_on_exit)

    # Set RNG up.
    #' @importFrom parallel clusterSetRNGStream
    clusterSetRNGStream(cl)

    # Set all packages to reload.
    ergm.MCMC.packagenames(reset=TRUE)
  }

  ergm.MCMC.packagenames(pending=control$MCMC.packagenames)
  
  # On the off chance that user wants to load extra packages which we don't know about already.
  for(pkg in ergm.MCMC.packagenames(pending=TRUE)){

    # Try loading from the same location as the master.
    #' @importFrom parallel clusterCall
    attached <- unlist(clusterCall(cl, require,
                                   package=pkg,
                                   character.only=TRUE,
                                   lib.loc=.libPaths()))
    # If something failed, warn and try loading from anywhere.
    if(!all(attached)){
      stop("Failed to attach package ", sQuote(pkg), " on one or more slave nodes. Make sure it's installed on or accessible from all of them and is in the library path.")
    }

    if(control$parallel.version.check){
      slave.versions <- clusterCall(cl,packageVersion,pkg)
      #' @importFrom utils packageVersion
      master.version <- packageVersion(pkg)

      wrong <- !sapply(slave.versions,identical,master.version)
      if(any(wrong))
        stop("The version of ", sQuote(pkg), " (", paste(slave.versions[wrong], collapse=", "), ") attached on one or more slave nodes is different from that on the (this) master node (", master.version,"). Make sure that the same version is installed on all nodes. If you are absolutely certain that this message is in error, override with the parallel.version.check=FALSE control parameter.")
    }
    ergm.MCMC.packagenames(loaded=pkg)
  }

  cl
}


#' @rdname ergm-parallel
#' @description The \code{ergm.stopCluster} shuts down a
#'   cluster, but only if `ergm.getCluster` was responsible for
#'   starting it.
#'
#' @param \dots not currently used
#' @export ergm.stopCluster
ergm.stopCluster <- function(..., verbose=FALSE){
  if(...length()) warning("Arguments to ergm.stopCluster() have been deprecated.")
  if(!is.null(ergm.cluster.started())){
    #' @importFrom parallel stopCluster
    if(verbose) message("Stopping the running cluster.")
    stopCluster(ergm.cluster.started())
    ergm.MCMC.packagenames(reset=TRUE)
    ergm.cluster.started(NULL)
  }else if(verbose>1) message("No cluster to stop.")
}

#' @rdname ergm-parallel
#' @description The \code{ergm.restartCluster} restarts and returns a cluster,
#'   but only if `ergm.getCluster` was responsible for starting it.
#'
#' @export ergm.restartCluster
ergm.restartCluster <- function(control=NULL, verbose=FALSE){
  if(!is.null(ergm.cluster.started())){
    if(verbose) message("Restarting the running cluster:")
    ergm.stopCluster(verbose=verbose)
  }else if(verbose>1) message("No cluster to restart.")
  ergm.getCluster(control, verbose=verbose, stop_on_exit=NULL) # stop_on_exit is already set by the initial cluster construction.
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

#' @rdname ergm-parallel
#'
#' @description `nthreads` is a simple generic to obtain the number of
#'   parallel processes represented by its argument, keeping in mind
#'   that having no cluster (e.g., `NULL`) represents one thread.
#'
#' @param clinfo a [`cluster`][parallel::makeCluster] or another object.
#' @export
nthreads <- function(clinfo=NULL, ...){
  UseMethod("nthreads")
}

#' @rdname ergm-parallel
#' @export
nthreads.cluster <- function(clinfo=NULL, ...){
  length(clinfo)
}

#' @rdname ergm-parallel
#' @export
nthreads.NULL <- function(clinfo=NULL, ...){
  NVL3(ergm.cluster.started(), nthreads(.), 1)
}

#' @rdname ergm-parallel
#' @method nthreads control.list
#' @export
nthreads.control.list <- function(clinfo=NULL, ...){
  if(is.numeric(clinfo$parallel)) return(max(1, clinfo$parallel))
  clinfo <- clinfo$parallel
  nthreads(clinfo, ...)
}
