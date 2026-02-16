#  File R/parallel.utils.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
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

#' Parallel Processing in the `ergm` Package
#'
#' Using clusters multiple CPUs or CPU cores to speed up ERGM
#' estimation and simulation.
#'
#' @details
#'
#' For estimation that require MCMC, [ergm][ergm-package] can take
#' advantage of multiple CPUs or CPU cores on the system on which it
#' runs, as well as computing clusters through one of two mechanisms:
#'
#' \describe{\item{Running MCMC chains in parallel}{ Packages
#' `parallel` and `snow` are used to to facilitate this, all cluster
#' types that they support are supported.
#' 
#' The number of nodes used and the parallel API are controlled using
#' the `parallel` and `parallel.type` arguments passed to the control
#' functions, such as [control.ergm()].
#'   
#' The [ergm.getCluster()] function is usually called internally by
#' the ergm process (in [ergm_MCMC_sample()]) and will attempt to
#' start the appropriate type of cluster indicated by the
#' [control.ergm()] settings. The [ergm.stopCluster()] is helpful if
#' the user has directly created a cluster.
#'     
#' Further details on the various cluster types are included below.}
#' 
#' \item{Multithreaded evaluation of model terms}{ Rather than running
#' multiple MCMC chains, it is possible to attempt to accelerate
#' sampling by evaluating qualified terms' change statistics in
#' multiple threads run in parallel. This is done using the
#' [OpenMP](https://www.openmp.org/) API.
#'
#' However, this introduces a nontrivial amont of computational
#' overhead. See below for a list of the major factors affecting
#' whether it is worthwhile.}}
#'
#' Generally, the two approaches should not be used at the same time
#' without caution. In particular, by default, cluster slave nodes
#' will not \dQuote{inherit} the multithreading setting; but
#' `parallel.inherit.MT=` control parameter can override that. Their
#' relative advantages and disadvantages are as follows:
#'
#' * Multithreading terms cannot take advantage of clusters but only
#'   of CPUs and cores.
#'
#' * Parallel MCMC chains produce several independent chains;
#'   multithreading still only produces one.
#' 
#' * Multithreading terms actually accellerates sampling, including
#'   the burn-in phase; parallel MCMC's multiple burn-in runs are
#'   effectively \dQuote{wasted}.
#' 
#' 
#' 
#' @section Different types of clusters:
#'
#' \describe{\item{PSOCK clusters}{ The `parallel` package is used with PSOCK clusters
#' by default, to utilize multiple cores on a system. The number of
#' cores on a system can be determined with the [detectCores()]
#' function.
#'   
#' This method works with the base installation of R on all platforms,
#' and does not require additional software.
#'   
#' For more advanced applications, such as clusters that span multiple
#' machines on a network, the clusters can be initialized manually,
#' and passed into [ergm()] and others using the `parallel` control
#' argument. See the second example below.}
#'
#' \item{MPI clusters}{ To use MPI to accelerate ERGM sampling,
#' pass the control parameter `parallel.type="MPI"`.
#' [ergm][ergm-package] requires the \CRANpkg{snow} and \CRANpkg{Rmpi} packages to
#' communicate with an MPI cluster.
#'   
#' Using MPI clusters requires the system to have an existing MPI
#' installation.  See the MPI documentation for your particular
#' platform for instructions.
#' 
#' To use [ergm()] across multiple machines in a high performance
#' computing environment, see the section "User initiated clusters"
#' below.}
#' 
#' \item{User initiated clusters}{ A cluster can be passed into [ergm()]
#' with the `parallel` control parameter. [ergm()] will detect the
#' number of nodes in the cluster, and use all of them for MCMC
#' sampling. This method is flexible: it will accept any cluster type
#' that is compatible with `snow` or `parallel` packages.
## #' Usage
## #' examples for a multiple-machine high performance MPI cluster can be
## #' found at the [Statnet
## #' wiki](https://statnet.csde.washington.edu/trac/wiki/ergmParallel).
#' }}
#'
#' @section When is multithreading terms worthwhile?:
#' 
#' * The more terms with statistics the model has, the more
#'   benefit from parallel execution.
#' 
#' * The more expensive the terms in the model are, the more benefit
#'   from parallel execution. For example, models with terms like
#'   [`gwdsp`][gwdsp-ergmTerm] will generally get more benefit than models where all
#'   terms are dyad-independent.
#'
#' * Sampling more dense networks will generally get more benefit than
#'   sparse networks. Network size has little, if any, effect.
#'
#' * More CPUs/cores usually give greater speed-up, but only up to a
#'   point, because the amount of overhead grows with the number of
#'   threads; it is often better to \dQuote{batch} the terms into a smaller
#'   number of threads than possible.
#'
#' * Any other workload on the system will have a more severe effect
#'   on multithreaded execution. In particular, do not run more
#'   threads than CPUs/cores that you want to allocate to the tasks.
#'
#' * Under Windows, even compiling with OpenMP appears to introduce
#'   unacceptable amounts of overhead, so it is disabled for Windows
#'   at compile time. To enable, *delete* `src/Makevars.win` and
#'   recompile from scratch.
#' 
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
#' @name ergm-parallel
#' @aliases parallel ergm.parallel parallel.ergm parallel-ergm
#'
NULL

#' @rdname ergm-parallel
#' @description The \code{ergm.getCluster} function is usually called
#'   internally by the ergm process (in
#'   [ergm_MCMC_sample()]) and will attempt to start the
#'   appropriate type of cluster indicated by the
#'   [control.ergm()] settings. It will also check that the
#'   same version of `ergm` is installed on each node.
#' @param control a [control.ergm()] (or similar) list of
#'   parameter values from which the parallel settings should be read;
#'   can also be [`NULL`], in which case an existing cluster is used
#'   if started, or no cluster otherwise.
#' @template verbose
#' @param stop_on_exit An [`environment`] or `NULL`. If an
#'   `environment`, defaulting to that of the calling function, the
#'   cluster will be stopped when the calling the frame in question
#'   exits.
#' 
#' @export ergm.getCluster
ergm.getCluster <- function(control=NULL, verbose=FALSE, stop_on_exit=parent.frame()){
  # If we don't want a cluster, just return NULL.
  if (is.numeric(control$parallel) && control$parallel==0) return(NULL)

  if(ERRVL2(get.MT_terms(), FALSE) && control$parallel.inherit.MT==FALSE) warning("Using term multithreading in combination with parallel MCMC is generally not advised. See help('ergm-parallel') for more information.")
  
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

    if(control$parallel.inherit.MT && ERRVL2(get.MT_terms(), 0)!=0){
      clusterCall(cl, set.MT_terms,
                  get.MT_terms())
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

#' [set.MT_terms()] controls multithreading of model terms.
#'
#' @param n an integer specifying the number of threads to use; 0 (the
#'   starting value) disables multithreading, and \eqn{-1} or
#'   \code{NA} sets it to the number of CPUs detected.
#' 
#' @return [set.MT_terms()] returns the previous setting, invisibly.
#'
#' @note The this is a setting global to the `ergm` package and all of
#'   its C functions, including when called from other packages via
#'   the `Linking-To` mechanism.
#'
#' @rdname ergm-parallel
#' @export
set.MT_terms <- function(n){
  old <- get.MT_terms()
  n <- as.integer(n)
  if(is.na(n) || n < -1) n <- -1L
  .Call("set_ergm_omp_terms", n, PACKAGE="ergm")
  invisible(old)
}

#' [get.MT_terms()] queries the current number of term calculation
#' threads.
#'
#' @return [get.MT_terms()] returns the current setting.
#'
#' @rdname ergm-parallel
#' @export
get.MT_terms <- function(){
  .Call("get_ergm_omp_terms", PACKAGE="ergm")
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


#' A rudimentary cache for large objects
#'
#' This cache is intended to store large, infrequently changing data
#' structures such as [`ergm_model`]s and [`ergm_proposal`]s on worker
#' nodes.
#'
#' @param comm a character string giving the desired function; see the
#'   default argument above for permitted values and Details for
#'   meanings; partial matching is supported.
#' @param key a character string, typically a `digest::digest()` of
#'   the object or a random string.
#' @param object the object to be stored.
#'
#'
#' Supported tasks are, respectively, to do nothing (the default),
#' return all entries (mainly useful for testing), clear the cache,
#' insert into cache, retrieve an object by key, check if a key is
#' present, or list keys defined.
#'
#' Deleting an entry can be accomplished by inserting a `NULL` for
#' that key.
#'
#' Cache is limited to a hard-coded size (currently 4). This should
#' accommodate an [`ergm_model`] and an [`ergm_proposal`] for
#' unconstrained and constrained MCMC. When additional objects are
#' stored, the oldest object is purged and garbage-collected.
#'
#' @note If called via, say, `clusterMap(cl, ergm_state_cache, ...)`
#'   the function will not accomplish anything. This is because
#'   `parallel` package will serialise the `ergm_state_cache()`
#'   function object, send it to the remote node, evaluate it there,
#'   and fetch the return value. This will leave the environment of
#'   the worker's `ergm_state_cache()` unchanged. To actually
#'   evaluate it on the worker nodes, it is recommended to wrap it in
#'   an empty function whose environment is set to [globalenv()]. See
#'   Examples below.
#'
#' @examples
#' \dontrun{
#' # Wrap ergm_state_cache() and call it explicitly from ergm:
#' call_ergm_state_cache <- function(...) ergm::ergm_state_cache(...)
#'
#' # Reset the function's environment so that it does not get sent to
#' # worker nodes (who have their own instance of ergm namespace
#' # loaded).
#' environment(call_ergm_state_cache) <- globalenv()
#'
#' # Now, call the the wrapper function, with ... below replaced by
#' # lists of desired arguments.
#' clusterMap(cl, call_ergm_state_cache, ...)
#' }
#'
#' @export
ergm_state_cache <- local({
  size <- 4
  objects <- list()

  function(comm = c("pass", "all", "clear", "insert", "get", "check", "list"),
           key, object){
    comm <- match.arg(comm)

    if(comm == "all") objects
    else if(comm == "clear") objects <<- list()
    else if(comm == "insert"){
      objects[[key]] <<- object # Add object to cache.
      if(length(objects) > size){
        objects <<- objects[-1] # If too many, clear the oldest.
        gc()
      }
      invisible(NULL)
    }
    else if(comm == "get") objects[[key]]
    else if(comm == "check") key %in% names(objects)
    else if(comm == "pass") invisible(NULL)
    else if(comm == "list") names(objects)
    else stop("Unknown command.")
  }
})

.GUID <- local({
  counter <- -.Machine$integer.max
  function() paste(c(unclass(Sys.time()), Sys.getpid(), (counter<<-counter+1)), collapse="-")
})
