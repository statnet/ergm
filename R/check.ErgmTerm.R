#  File R/check.ErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
#====================================================================================
# This file contains the following 6 files that help check the validity of ergm terms
#       <check.ErgmTerm>               <get.InitErgm.fname>
#       <assignvariables>              <zerowarnings>
#       <check.ErgmTerm.summarystats>  <extremewarnings>
#====================================================================================




######################################################################################
# The <check.ErgmTerm> function ensures for the <InitErgmTerm.X> function that the
# term X:
#   1) is applicable given the 'directed' and 'bipartite' attributes of the given
#      network
#   1.5) is not applied to a directed bipartite network
#   2) has an appropiate number of arguments
#   3) has correct argument types if arguments where provided
#   4) has default values assigned if defaults are available
# by halting execution if any of the first 3 criteria are not met
#
# --PARAMETERS--
#  nw           : the network that term X is being checked against  
#  arglist      : the list of arguments for term X
#  directed     : whether term X requires a directed network (T or F); default=NULL
#  bipartite    : whether term X requires a bipartite network (T or F); default=NULL
#  nonnegative  : whether term X requires a network with only nonnegative weights; default=FALSE
#  varnames     : the vector of names of the possible arguments for term X;
#                 default=NULL 
#  vartypes     : the vector of types of the possible arguments for term X;
#                 default=NULL 
#  defaultvalues: the list of default values for the possible arguments of term X;
#                 default=list()
#  required     : the logical vector of whether each possible argument is required;
#                 default=NULL
#
# --RETURNED--
#   out: a list of the values for each possible argument of term X; user provided 
#        values are used when given, default values otherwise.
#
######################################################################################
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
#'   term X; default=NULL
#' @param defaultvalues the list of default values for the possible
#'   arguments of term X; default=list()
#' @param required the logical vector of whether each possible
#'   argument is required; default=NULL
#' @template response
#' @return A list of the values for each possible argument of term X;
#'   user provided values are used when given, default values
#'   otherwise.
#'
#' @import network
#' @export check.ErgmTerm
check.ErgmTerm <- function(nw, arglist, directed=NULL, bipartite=NULL, nonnegative=FALSE,
                           varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL, response=NULL) {
  stopifnot(all_identical(c(length(varnames), length(vartypes), length(defaultvalues), length(required))))
  message <- NULL
  if (!is.null(directed) && directed != (dnw<-is.directed(nw))) {
    #directed != (dnw<-eval(expression(nw$gal$dir),parent.frame()))) {
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
  if (is.null(message) && nonnegative && any(nw %e% response < 0)){
    message <- "networks with negative dyad weights"
  }
  if (!is.null(message)) {
    ergm_Init_abort("Term may not be used with ",message,".")
  }

  sr=sum(required)
  lv=length(varnames)
  la=length(arglist)
  if(la < sr || la > lv) {
    if (sr < lv)
      expected = paste("from",sr,"to",lv,"arguments,")
    else if(sr==1)
      expected = "1 argument,"
    else
      expected = paste(sr,"arguments,")
    ergm_Init_abort("Model term expected ", expected, " got ", la, '.')
  }
# The correctness of what the user typed is checked, but it is assumed
# that each InitErgmTerm function faithfully passes in what the user typed;
# thus, the correctness of input from the InitErgmTerm function isn't checked.
  out = defaultvalues
  names(out)=varnames
  m=NULL
  if (la>0) {
    for(i in 1:la) { # check each arglist entry
      if (!is.null(names(arglist)) && (name <- names(arglist)[i]) != "") {
        m = pmatch(name, varnames)# try to match user-typed name if applicable
        if(is.na(m)) { # User typed an unrecognizable name
          ergm_Init_abort("Model term does not recognize ", sQuote(name), " argument.")
        }
        # valid name match with mth variable if we got to here
        if (all(sapply(strsplit(vartypes[m],",",fixed=TRUE)[[1]], function(vartype) !is.null(arglist[[i]]) && !is(arglist[[i]], vartype)))) {
          # Wrong type
          ergm_Init_abort(sQuote(name), " argument is not of the expected ", sQuote(vartypes[m]), " type.")
        }
        # correct type if we got to here
        out[m] <- list(arglist[[i]])
      } else { # no user-typed name for this argument
        if (!is.null(m)) {
          ergm_Init_abort("Unnamed argument follows named argument.")
        }
        if (all(sapply(strsplit(vartypes[i],",",fixed=TRUE)[[1]], function(vartype) !is.null(arglist[[i]]) && !is(arglist[[i]], vartype)))) {
          # Wrong type
          ergm_Init_abort("Argument number ", i, " is not of the expected ", sQuote(vartypes[i]), " type.")
        }
        # correct type if we got to here
        out[i] <- list(arglist[[i]])
      }
    }
  }
  #  c(.conflicts.OK=TRUE,out)
  out
}
