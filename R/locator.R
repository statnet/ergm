#  File R/locator.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#' A simple dictionary to cache recent InitFunction lookups.
#'
#' @param name function name.
#' @param env the environment name for the function; if `NULL`, look
#'   up in cache, otherwise insert or overwrite.
#'
#' @return A character string giving the name of the environment
#'   containing the function, or `NULL` if not in cache.
#' @noRd
locator_cache <- local({
  cache <- list()
  watchlist <- character(0) # Packages being watched for unloading.
  pkglist <- character(0) # Current list of packages.
  # Reset the cache and update the list of watched packages.
  reset <- function(...){
    pkglist <<- .packages()
    new <- setdiff(pkglist, watchlist)
    for(pkg in new){
      setHook(packageEvent(pkg, "detach"), reset)
      setHook(packageEvent(pkg, "onUnload"), reset)
    }
    watchlist <<- c(watchlist, new)
    cache <<- list()
  }
  # Check if new namespaces have been added.
  checknew <- function(){
    if(!setequal(.packages(), pkglist)) reset()
  }
  function(name, env=NULL){
    checknew()
    if(is.null(env)){
      cache[[name]]
    }else{
      cache[[name]] <<- env
    }
  }
})

locate.InitFunction <- function(name, prefix, errname=NULL, env = globalenv()){
  if(is.call(name)) name <- name[[1]]
  name <- as.character(name)
  fname <- paste(prefix,name,sep=".")
  
  # Try the given environment...
  if(!is.null(obj<-get0(fname, mode='function', envir=env))){
    env <- environment(obj)
    envname <- environmentName(env)
    # Check that environment name is not blank or globalenv(), and
    # that the detected environment actually contains the object.
    if(! NVL(envname,"") %in% c("", "R_GlobalEnv") && exists(fname, mode='function', envir=env, inherits=FALSE)) return(call(":::",as.name(envname),as.name(fname)))
    else return(as.name(fname))
  }

  # Try the cache...
  envname <- locator_cache(fname)
  if(!is.null(envname)) return(call(":::",as.name(envname),as.name(fname)))

  # Use getAnywhere()...
  #' @importFrom utils getAnywhere
  m <- getAnywhere(fname)
  if(length(m$objs)){
    ## Prioritise visible over not:
    if(any(m$visible)){
      m <- lapply(m[-1], "[", m$visible)
    }
    if(length(m$objs)>1) warning("Name ",fname," matched by multiple objects; using the first one on the list.")
    envname <- environmentName(environment(m$objs[[1]]))
    locator_cache(fname, envname)
    return(call(":::",as.name(envname),as.name(fname)))
  }

  # If not found, error or NULL...
  if(!is.null(errname)) stop(errname,' ', sQuote(name), " initialization function ", sQuote(fname), " not found.")
  else NULL
}

