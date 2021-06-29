#  File R/check.ErgmTerm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

ergm_Init_inform_once <- once(ergm_Init_inform)
ergm_Init_warn_once <- once(ergm_Init_warn)

#' Ensures an Ergm Term and its Arguments Meet Appropriate Conditions
#'
#' Helper functions for implementing \code{\link[=ergm]{ergm()}}
#' terms, to check whether the term can be used with the specified
#' network.  For information on ergm terms, see
#' \link{ergm-terms}. \code{ergm.checkargs},
#' \code{ergm.checkbipartite}, and \code{ergm.checkderected} are
#' helper functions for an old API and are deprecated. Use
#' \code{check.ErgmTerm}.
#'
#' The \code{check.ErgmTerm} function ensures for the
#' \code{\link{InitErgmTerm}.X} function that the term X: \itemize{
#' \item is applicable given the 'directed' and 'bipartite' attributes
#' of the given network \item is not applied to a directed bipartite
#' network \item has an appropiate number of arguments \item has
#' correct argument types if arguments where provided \item has
#' default values assigned if defaults are available } by halting
#' execution if any of the first 3 criteria are not met.
#'
#' @param nw the network that term X is being checked against
#' @param arglist the list of arguments for term X
#' @param directed logical, whether term X requires a directed
#'   network; default=NULL
#' @param bipartite whether term X requires a bipartite network (T or
#'   F); default=NULL
#' @param nonnegative whether term X requires a network with only
#'   nonnegative weights; default=FALSE
#' @param varnames the vector of names of the possible arguments for
#'   term X; default=NULL
#' @param vartypes the vector of types of the possible arguments for
#'   term X, separated by commas; an empty string (`""`) or `NA` disables the check for that argument; default=NULL
#' @param defaultvalues the list of default values for the possible
#'   arguments of term X; default=list()
#' @param required the logical vector of whether each possible
#'   argument is required; default=NULL
#' @param dep.inform,dep.warn a list of length equal to the number of
#'   arguments the term can take; if the corresponding element of the
#'   list is not `FALSE`, a [message()] or a [warning()] respectively
#'   will be issued if the user tries to pass it; if the element is a
#'   character string, it will be used as a suggestion for
#'   replacement.
#' @return A list of the values for each possible argument of term X;
#'   user provided values are used when given, default values
#'   otherwise. The list also has an `attr(,"missing")` attribute
#'   containing a named logical vector indicating whether a particular
#'   argument had been set to its default.
#'
#' @import network
#' @export check.ErgmTerm
check.ErgmTerm <- function(nw, arglist, directed=NULL, bipartite=NULL, nonnegative=FALSE,
                           varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL, dep.inform=rep(FALSE, length(required)), dep.warn=rep(FALSE, length(required))){
  # Ensure that all inputs are of the correct type.
  ergm_Init_try(arglist <- as.list(arglist))
  varnames <- as.character(varnames)
  vartypes <- as.character(vartypes)
  defaultvalues <- as.list(defaultvalues)
  required <- as.logical(required)
  dep.inform <- as.list(dep.inform)
  dep.warn <- as.list(dep.warn)

  stopifnot(all_identical(c(length(varnames), length(vartypes), length(defaultvalues), length(required), length(dep.inform), length(dep.warn))))
  message <- NULL
  if (!is.null(directed) && directed != (dnw<-is.directed(nw))) {
    message <- paste("networks with directed==", dnw, sep="")
  }
  
  bnw<- nw %n% "bipartite"
  # check for bipartite 1st partition size zero (not yet supported by ergm)
  if(is.numeric(bnw)){
    if (bnw==0){
      message <- "networks with a bipartite first partition of size 0 (bipartite=0)"
    }
  }
  
  if(is.null(bnw)) bnw <- 0
  if (!is.null(bipartite) && bipartite != (bnw > 0)) {
    #bipartite != (bnw <- eval(expression(nw %n% "bipartite"),parent.frame()) > 0)) {
    message <- paste("networks with bipartite", 
                     ifelse(bnw>0, " > 0", "==FALSE"), sep="")
  }
  if (is.directed(nw) && bnw > 0) {
    message <- "directed bipartite networks"
  }
  if (is.null(message) && nonnegative && any(nw %e% (nw%ergmlhs%"response") < 0)){
    message <- "networks with negative dyad weights"
  }
  if (!is.null(message)) {
    ergm_Init_abort("Term may not be used with ",message,".")
  }

  # Construct a dummy function that copies all its arguments into a
  # list and sets an attribute indicating whether they are missing.
  f <- function(){
    ..n <- names(formals())
    ..l <- structure(vector("list", length(..n)), names=..n)
    ..m <- structure(logical(length(..n)), names=..n)
    for(..arg in ..n){
      ..m[..arg] <- do.call(missing,list(as.name(..arg)))
      ..l[..arg] <- list(get(..arg, inherits=FALSE))
    }
    structure(..l, missing=..m)
  }

  # Set the argument names and their defaults (if not required).
  formals(f) <- replace(structure(defaultvalues, names = varnames), required, list(quote(expr=)))
  # Now, try calling it with the arglist.
  ergm_Init_try(out <- do.call(f, arglist, envir=baseenv(), quote=TRUE))
  # out is now a list containing elements of arglist in the correct order and defaults filled in.

  for(m in seq_along(out)){
    name <- names(out)[m]
    miss <- attr(out, "missing")[m]
    val <- out[[m]]

    # Check type
    if(!is.na(vartypes[m]) && nchar(vartypes[m])
       && all(sapply(strsplit(vartypes[m],",",fixed=TRUE)[[1]], function(vartype) !is.null(val) && !is(val, vartype))))
      ergm_Init_abort(sQuote(name), " argument is not of any of the expected (", vartypes[m], ") types.")

    # Check deprecation (but only if passed explicitly)
    if(!miss){
      if(!isFALSE(dep.inform[[m]])) {
        if(is.character(dep.inform[[m]]))
          ergm_Init_inform_once("Argument ", sQuote(varnames[m]), " has been superseded by ", sQuote(dep.inform[[m]]), ", and it is recommended to use the latter.  Note that its interpretation may be different.")
        else
          ergm_Init_inform_once("Argument ", sQuote(varnames[m]), " has been deprecated and may be removed in a future version.")
      }
      if(!isFALSE(dep.warn[[m]])) {
        if(is.character(dep.warn[[m]]))
          ergm_Init_warn_once("Argument ", sQuote(varnames[m]), " has been deprecated and may be removed in a future version.  Use ", sQuote(dep.warn[[m]]), " instead.  Note that its interpretation may be different.")
        else
          ergm_Init_warn_once("Argument ", sQuote(varnames[m]), " has been deprecated and may be removed in a future version.")
      }
    }
  }

  out
}
