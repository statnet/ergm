#  File R/get.node.attr.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
###############################################################################
# The <get.node.attr> function returns the vector of covariates for the given
# network and specified attribute if the attribute exists - execution will
# halt if the attribute is not correctly given as a single string or is not 
# found in the vertex attribute list; optionally <get.node.attr> will also 
# check that return vector is numeric, halting execution if not
#
# --PARAMETERS--
#   nw          : a network object
#   attrname    : the name of a nodal attribute, as a character string
#   functionname: the name of the calling function; this is only used for
#                 the warning messages that accompany a halt
#   numeric     : whether to halt execution if the return vector is not
#                 numeric; default=FALSE
#   
# --RETURNED--
#   out:  the vector of 'attrname' covariates
#
###############################################################################



#' Retrieve and check assumptions about vertex attributes (nodal covariates) in
#' a network
#' 
#' The \code{get.node.attr} function returns the vector of nodal covariates for
#' the given network and specified attribute if the attribute exists -
#' execution will halt if the attribute is not correctly given as a single
#' string or is not found in the vertex attribute list; optionally
#' \code{get.node.attr} will also check that return vector is numeric, halting
#' execution if not. The purpose is to validate assumptions before passing
#' attribute data into an ergm term.
#' 
#' 
#' @param nw a \code{\link{network}} object
#' @param attrname the name of a nodal attribute, as a character string
#' @param functionname the name of the calling function a character string;
#' this is only used for the warning messages that accompany a halt
#' @param numeric logical, whether to halt execution if the return vector is
#' not numeric; default=FALSE
#' @return returns the vector of 'attrname' covariates for the vertices in the
#' network
#' @seealso \code{\link[network]{get.vertex.attribute}} for a version without
#' the checking functionality
#' @examples
#' 
#' data(faux.mesa.high)
#' get.node.attr(faux.mesa.high,'Grade')
#' 
#' @export get.node.attr
get.node.attr <- function(nw, attrname, functionname=NULL, numeric=FALSE) {  

# This is a kludge, which has been patched to bring it in line with the
# corrected class definitions.  -CTB

  if (is.null(functionname)) {
    # Assume it's being called from InitErgm.* or InitErgmTerm.*
    # Otherwise,  get.InitErgm.fname() will return NULL
    functionname <- get.InitErgm.fname()
    functionname <- ifelse (is.null(functionname), "unknown function",
                        sub('.*[.]', '', functionname)) # truncate up to last '.'
  }
  if (!is.character(attrname) || length(attrname)>1){
    stop(paste("The argument", attrname, "passed to", functionname,
               "must be a single character string naming a nodal attribute."),
         call.=FALSE)
  }
  #We'll assume that every vertex must have a value, so checking the 
  #first is reasonable.  #skye:  I think this is not a correct assumption
	
  if(NVL(get.network.attribute(nw,"bipartite"),FALSE)){}  
  # not sure what this code above was supposed to do, but it was only checking if the attribute existed if the network	was bipartite
  # maybe, if it was bipartite, it should check if appropriate values exist for the mode in question? and should use 'is.bipartite()' for the check
  
  if (!any(attrname==unique(unlist(lapply(nw$val,names))))){
    stop(paste("Attribute '", attrname, "' named in", functionname,
               "model term is not a vertex attribute  of the network."),
         call.=FALSE)
  }
  #"[["(nw$val,attrname)
  out <- unlist(get.vertex.attribute(nw,attrname))
  if(numeric && !is.numeric(out)) {
    stop("The ", attrname, " attribute for the ", functionname, 
         " term is not numeric as required.", call.=FALSE)
  }
  out
}

.despace <- function(s) gsub("[[:space:]]", "", s)

#' @name node-attr
#' @title Specifying nodal attributes and their levels
#'
#' This document describes both the ways in which to specify nodal
#' attribute or functions and which levels for categorical factors to
#' include, and the helper functions for use in `InitErgmTerm`
#' implementations.
#'
#' @details
#' 
#' Term nodal attribute arguments, typically called `attrname`, `by`,
#' `on`, etc. are interpreted as follows: \describe{
#' 
#' \item{a single character string}{Extract the vertex attribute.}
#' 
#' \item{a character vector of length > 1}{Extract the vertex
#' attributes and pastes them together, separated by dots.}
#' 
#' \item{a function}{The function is called on the LHS network,
#' expected to return a vector of appropriate length. (Shorter vectors
#' will be recycled as needed.)}
#' 
#' \item{a formula}{The expression on the RHS of the formula is
#' evaluated in an environment of the vertex attributes of the
#' network, expected to return a vector of appropriate
#' length. (Shorter vectors will be recycled as needed.) Within this
#' expression, the network itself accessible as either `.` or
#' `.nw`. For example,
#' `nodecov(~abs(Grade-mean(Grade))/network.size(.))` would return the
#' absolute difference of each actor's "Grade" attribute from its
#' network-wide mean, divided by the network size.}
#' 
#' }
#'
#' For categorical attributes, to select which levels are of interest
#' and their ordering, use the argument `levels`.  It is interpreted
#' as follows: \describe{
#'
#' \item{a vector}{Use these levels in the order specified; discard
#' others.}
#' 
#' \item{a function}{The function is called on the list of unique
#' values of the attribute, the values of the attribute themselves,
#' and the network itself, depending on its arity. Its return value
#' may be a vector of levels to use.}
#'
#' \item{a formula}{The expression on the RHS of the formula is
#' evaluated in an environment in which the network itself is
#' accessible as `.nw`, the list of unique values of the attribute as
#' `.` or as `.levels`, and the attribute vector itself as
#' `.attr`. For example, `` would return the absolute difference of
#' each actor's "Grade" attribute from its network-wide mean, divided
#' by the network size. Its return value may be a vector of levels to
#' use.}
#' 
#' }
#' 
#' Note that `levels` often has a default that is sensible for the
#' term in question.
#' 
#' @aliases attrname on by
NULL

#' @name node-attr-api
#' @title Helper functions for specifying nodal attribute levels
#'
#' These functions are meant to be used in `InitErgmTerm` and other
#' implementations to provide the user with a way to extract nodal
#' attributes and select their levels in standardized and flexible
#' ways described under [`node-attr`].
#'
#' @param object An argument specifying the nodal attribute to select
#'   or which levels to include.
#' @param nw Network on the LHS of the formula.
#' @param attr A vector of length equal to the number of nodes,
#'   specifying the attribute vector.
#' @param levels Starting set of levels to use; defaults to the sorted
#'   list of unique attributes.
#' @param type Permitted types of the attribute; no checking by
#'   default.
#' @param bip Bipartedness mode: affects either length of attribute
#'   vector returned or the length permited: `"n"` for full network, `"b1"`
#'   for first mode of a bipartite network, and `"b2"` for the second.
#' @param ... Additional argument to the functions of network or to
#'   the formula's environment.
#'
#'
NULL

#' @rdname node-attr-api
#'
#' @description `ergm_get_vattr` extracts and processes the specified
#'   nodal attribute vector. It is strongly recommended that
#'   [check.ErgmTerm()]'s corresponding
#'   `vartype="function,formula,character"`.
#' 
#' @return A vector of length equal to the number of nodes giving the
#'   selected attribute function. It may also have an attribute
#'   `"name"`, which controls the suggested name of the attribute
#'   combination.
#' @export
ergm_get_vattr <- function(object, nw, bip=c("n","b1","b2"), ...){
  bip <- match.arg(bip)
  UseMethod("ergm_get_vattr")
}

.rightsize_vattr <- function(a, nw, bip){
  rep_len_warn <- function(x, length.out){
    if(length.out%%length(x)) warning("Length of vertex attribute vector is not a multiple of network size or bipartite group size.")
    rep_len(x, length.out)
  }
  if(!is.bipartite(nw) || bip=="n") rep_len_warn(a, network.size(nw))
  else if(length(a)==network.size(nw)) # Input vector is n-long, need to trim.
    switch(bip,
           b1 = a[seq_len(nw%n%"bipartite")],
           b2 = a[-seq_len(nw%n%"bipartite")])
  else # Othewise, recycle until the right length.
    rep_len_warn(a, switch(bip,
                           b1=nw%n%"bipartite",
                           b2=network.size(nw)-nw%n%"bipartite"))
}


#' @importFrom purrr "%>%" "map" "pmap_chr"
#' @export
ergm_get_vattr.character <- function(object, nw, bip=c("n","b1","b2"), ...){
  if(!all(object%in%list.vertex.attributes(nw))) stop(sQuote(object), " is not a valid nodal attribute.")

  o <- (if(length(object)==1) nw%v%object
        else object %>% map(~nw%v%.) %>% pmap_chr(paste, sep=".")) %>%
    .rightsize_vattr(nw, bip)
  
  attr(o, "name") <- paste(object, collapse=".")
  o
}


#' @export
ergm_get_vattr.function <- function(object, nw, bip=c("n","b1","b2"), ...){
  object(nw, ...) %>%
    .rightsize_vattr(nw, bip)
}


#' @importFrom purrr "%>%" "map" "set_names"
#' @importFrom tibble "lst"
#' @export
ergm_get_vattr.formula <- function(object, nw, bip=c("n","b1","b2"), ...){
  a <- list.vertex.attributes(nw)
  vlist <- c(a %>% map(~nw%v%.) %>% set_names(a),
             lst(`.`=nw, .nw=nw, ...))

  e <- object[[length(object)]]
  o <- eval(e, envir=vlist, enclos=environment(object)) %>%
    .rightsize_vattr(nw, bip)
  
  if(length(object)>2) attr(o, "name") <- eval_lhs.formula(e)
  o
}

#' @rdname node-attr-api
#'
#' @description `ergm_attr_levels` filters the levels of the
#'   attribute.  It is strongly recommended that [check.ErgmTerm()]'s
#'   corresponding
#'   `vartype="function,formula,character,numeric,logical,NULL"`
#' 
#' @return A vector of levels to use and their order.
#' @export
ergm_attr_levels <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  UseMethod("ergm_attr_levels")
}

#' @export
ergm_attr_levels.default <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  object
}

#' @export
ergm_attr_levels.NULL <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  levels
}

#' @export
ergm_attr_levels.function <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  object <- if('...' %in% names(formals(object))) object(levels, attr, nw, ...)
            else switch(length(formals(object)),
                        object(levels),
                        object(levels, attr),
                        object(levels, attr, nw))
  ergm_attr_levels(object, attr, nw, levels, ...)
}

#' @export
ergm_attr_levels.formula <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  vlist <- lst(`.`=levels, .levels=levels, .attr=attr, .nw=nw, ...)
  e <- object[[length(object)]]
  object <- eval(e, envir=vlist, enclos=environment(object))  
  ergm_attr_levels(object, attr, nw, levels, ...)
}
