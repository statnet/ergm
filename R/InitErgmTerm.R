#  File R/InitErgmTerm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch


#===========================================================================
# This file contains the following 74 new, easier-to-write ERGM term
# initialization functions (each prepended with "InitErgmTerm"):
#   A:   <absdiff>          <absdiffcat>      <altkstar>
#        <asymmetric>       <adegcor>
#   B:   <b1concurrent>     <b1degree>        <b1factor>
#        <b1star>           <b1starmix>       <b1twostar>
#        <b2concurrent>     <b2degree>        <b2factor>         
#        <b2star>           <b2starmix>       <b2twostar>
#        <balance>
#   C:   <concurrent>       <cycle>           <ctriple>=<ctriad> 
#   D:   <degree>           <degreepopularity><density>         <dsp>
#        <dyadcov>          <degcrossprod>    <degcor>
#   E:   <edgecov>          <edges>           <esp>
#   G:   <gwb1degree>       <gwb2degree>      <gwdegree>
#        <gwdsp>            <gwesp>           <gwidegree>
#        <gwnsp>            <gwodegree>
#   H:   <hamming>
#   I:   <idegree>          <intransitive>    <idegreepopularity> 
#        <isolates>         <istar>
#   K:   <kstar>
#   L:   <localtriangle>
#   M:   <m2star>           <meandeg>         <mutual>
#   N:   <nearsimmelian>    <nodefactor>      <nodecov>=<nodemain> 
#        <nodeicov>         <nodeifactor>     <nodematch>=<match>
#        <nodemix>          <nodeocov>        <nodeofactor>       
#        <nsp>
#   O:   <odegree>          <opentriad>       <ostar>
#        <odegreepopularity>  
#   P:   <pdegcor>
#   R:   <receiver>         <rdegcor>
#   S:   <sender>           <simmelian>       <simmelianties>
#        <smalldiff>        <sociality>
#   T:   <threepath>        <transitive>      <triangles>=<triangle>
#        <triadcensus>      <tripercent>      <ttriple>=<ttriad>
#        <transitiveties>   <twopath
#==========================================================================

################################################################################
# The <InitErgmTerm.X> functions initialize each ergm term, X, by
#   1) checking the validity of X and its arguments via <check.ErgmTerm> and
#   2) setting appropiate values for each of the components in the returned list
# X is initialized for inclusion into a model that is specified by formula F and
# built via <ergm_model>
# 
# --PARAMETERS--
#   nw        : the network given in formula F
#   arglist   : the arguments given with term X in formula F
#
# --IGNORED PARAMETERS--
#   ... : ignored, but necessary to accomodate other arguments
#         passed by <ergm_model>
#
# --RETURNED--
#   a list of term-specific elements required by the C changestats
#   functions and other R rountines; the first two components of this
#   list are required*, the remaining components are optional:
#     *name      : the name of term X; this is used to locate the C function
#                  calculating the change statistics for X, which will be
#                  'name' prepended with "d_"; for example if X=absdiff,
#                  'name'="absdiff", and the C function is "d_absdiff"
#     *coef.names: the vector of names for the coefficients (parameters)
#                  as they will be reported in the output
#     inputs     : the vector of (double-precision numeric) inputs that the 
#                  changestat function called d_<name> will require
#                  (see WHAT THE C CHANGESTAT FUNCTION RECEIVES below);
#                  this MUST be a vector!; thus, if the inputs are  matrices,
#                  they must be "flattened" to vectors; if they are categorical
#                  character-valued variables, they must be converted to numbers;
#                  optionally, 'inputs' may have an attribute named
#                  "ParamsBeforeCov",which is the number that used to be the
#                  old Element 1 and is needed for backwards compatability
#                  (see the old <InitErgm> for details); default=NULL
#     soname     : the name of the package containing the C function called
#                  d_'name'; default="ergm"
#     dependence : whether the addition of term X to the model makes the model
#                  into a dyadic dependence model (T or F); if all terms have
#                  'dependence' set FALSE, the model is assumed to be a
#                  dyadic independence model; default=TRUE
#    emptynwstats: the vector of values (if nonzero) for the statistics evaluated
#                  on the empty network; if all are zero for this term, this
#                  argument may be omitted.  Example:  If the degree0 term is
#                  among the statistics, this argument is unnecessary because
#                  degree0 = number of nodes for the empty network
#    minpar      : the vector of minimal valid values for each of the model's parameters
#    maxpar      : the vector of maximal valid values for each of the model's parameters
#    params      : a list whose names correspond to parameter values for curved exponential family model
#                  terms only; the items in the list are there for historical reasons and are ignored;
#    map         : a function taking two arguments, theta and length('params'), which
#                  gives the map from the canonical parameters, theta, to the curved
#                  parameters, eta; 'map' is only necessary for curved exponential
#                  family model terms
#   gradient     : a function taking two arguments, theta and length('params'), which
#                  gives the gradient of the eta map above as a p by q matrix, where
#                  p=length(theta), q=length(params); 'gradient' is only necessary
#                  for curved exponential family model terms
#   offset       : a logical value; if TRUE, forces the term to be an offset
#   offsettheta  : a logical vector length equal to the number of parameters; if TRUE,
#                  the corresponding parameter is forced to be an offset
#   offsetmap    : a logical vector length equal to the number of statistics; if TRUE,
#                  the corresponding statistic is forced to be an offset
# WHAT THE C CHANGESTAT FUNCTION RECEIVES:
#                The changestat function, written in C and called d_'name',
#                will have access to 'inputs'; this array will be called INPUT_PARAMS
#                in the C code and its entries may accessed as INPUT_PARAMS[0],
#                INPUT_PARAMS[1], and so on; the size of INPUT_PARAMS=N_INPUT_PARAMS,
#                a value which is automatically set for you and which is available
#                inside the C function; thus INPUT_PARAMS[N_INPUT_PARAMS-1] is the last
#                element in the vector; note in particular that it is NOT necessary 
#                to add the number of inputs to 'inputs' since this is done automatically
#
################################################################################

#' @importFrom utils hasName

#=======================InitErgmTerm utility functions============================#

#' @title Curved settings for geometric weights for the `gw*` terms
#'
#' @description This is a list containing `map` and `gradient` for the weights described by Hunter (2007).
#'
#' @references
#' David R. Hunter (2007) Curved Exponential Family Models for Social Networks. *Social Networks*, 29: 216-230. \doi{10.1016/j.socnet.2006.08.005}
#' @keywords internal
#' @name ergm_GWDECAY
#' @export ergm_GWDECAY
ergm_GWDECAY <- list(
  map = function(x,n,...) {
    i <- 1:n
    x[1] * ifelse(i==1, 1, (exp(x[2])*(1-(1-exp(-x[2]))^i)))
  },
  gradient = function(x,n,...) {
    i <- 1:n
    e2 <- exp(x[2])
    a <- 1-exp(-x[2])
    rbind((1-a^i)*e2, ifelse(i==1, 0, x[1] * ( (1-a^i)*e2 - i*a^(i-1) ) ) )
  },
  minpar = c(-Inf, 0)
)

GWDECAY <- ergm_GWDECAY

.spcache.aux <- function(type){
  type <- toupper(type)
  trim_env(as.formula(as.call(list(as.name('~'), as.call(list(as.name('.spcache.net'),type=if(type=='ITP')'OTP' else type))))))
}

nodecov_names <- function(nodecov, prefix=NULL){
  cn <- if(is.matrix(nodecov)){
          cn <- colnames(nodecov)
          if(is.null(cn) || all(cn==seq_along(cn))) paste(attr(nodecov, "name"), seq_len(ncol(nodecov)), sep=".")
          else cn
        }else attr(nodecov, "name")
  NVL3(prefix, paste0(prefix,".",cn), cn)
}

#' A common pattern for obtaining an edge covariate
#'
#' @param name a string containing the name of the calling term.
#' @param nw the LHS network.
#' @param a list returned by [check.ErgmTerm()].
#'
#' @return A list with two elements: `xm` for the obtained predictor
#'   matrix and `cn` for the standard coefficient name.
#' @keywords internal
#' @export
ergm_edgecov_args <- function(name, nw, a){
  # Process the arguments
  if(is.network(a$x)){
    if(!is.null(a$attrname) && !a$attrname %in% list.edge.attributes(a$x))
      ergm_Init_abort("Specified network ", sQuote(deparse1(attr(a,"exprs")$x)), " does not have an edge attribute ", sQuote(a$attrname), ".")
    xm <- as.matrix(a$x, matrix.type="adjacency", a$attrname)
  }else if(is.character(a$x)){
    xm <- get.network.attribute(nw, a$x)
    if(is.null(xm) || !is.matrix(xm)) ergm_Init_abort("There is no network attribute named ", sQuote(a$x), " or it is not a matrix.")
    if(is.network(xm)){
      if(!is.null(a$attrname) && !a$attrname %in% list.edge.attributes(xm)) ergm_Init_abort("Network at attribute named ", sQuote(a$x), " does not have an edge attribute ", sQuote(a$attrname), ".")
      xm <- as.matrix(xm, matrix.type="adjacency", attrname=a$attrname)
      name <- paste(name, a$x, sep=".")
    }
  }else xm <- as.matrix(a$x)

  cn <- if(!is.null(a$attrname)) paste(name, a$attrname, sep = ".")
        else paste(name, if(is.character(a$x)) a$x else deparse1(attr(a,"exprs")$x), sep = ".")

  list(xm=xm, cn=cn)
}

# LEVELS_BASE1 is a placeholder for whatever the value of levels= or
# nodes= should be when base==1. For now, it's NULL to prevent the two
# arguments from interfering. Eventually, when base= is removed, it
# will need to be set to -1 either here or by search-and-replace.
LEVELS_BASE1 <- NULL

decay_vs_fixed <- function(a, name, no_curved_attrarg=TRUE){
  if(!is.null(a$alpha)){
    ergm_Init_abort("For consistency with ", sQuote("gw*degree"), " terms, in all ", sQuote("gw*sp"), " and ", sQuote("dgw*sp"), " terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".")
  }

  if(a$fixed){
    if(is.null(a$decay)) ergm_Init_abort("Using ", sQuote('fixed=TRUE')," requires a decay parameter ", sQuote('decay'), ".")
  }else{
    if(!is.null(a$decay)) ergm_Init_warn("Decay parameter ", sQuote('decay')," passed with ", sQuote('fixed=FALSE'), ". ", sQuote('decay')," will be ignored. To specify an initial value for ", sQuote('decay'),", use the ", sQuote('control.ergm()'), " parameter ", sQuote('init='), ".")

    if(no_curved_attrarg && !is.null(NVL(a$attrname,a$attr))) ergm_Init_abort("Using ", sQuote('fixed=FALSE'), " with an attribute is not implemented at this time. Use ", sQuote('fixed=TRUE'), ".")
  }
}

.degrange_impl <- function(deg, dir, bip, nw, arglist, ..., version){
  termname <- paste0(deg, "degrange")
  coefpre <- paste0(deg, "deg")

  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=dir, bipartite=bip,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=dir, bipartite=bip,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels
  }
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term ", sQuote(termname), " must have arguments either of the same length, or one of them must have length 1.")

  to <- ifelse(to==Inf, pmax(from, network.size(nw)) + 1, to)
  if(any(from>=to)) ergm_Init_abort("Term ", sQuote(termname), " must have ", sQuote("from"), "<", sQuote("to"), ".")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- if(NVL(bip, FALSE)) ergm_get_vattr(byarg, nw, bip = if(homophily) "n" else deg) else ergm_get_vattr(byarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degrange and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3,du[1,]==0]
      for(i in seq_along(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(from==0)) {
      emptynwstats <- rep(0, length(from))
      emptynwstats[from==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(from)==0){return(NULL)}
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste(coefpre, from,"+",sep=""),
                         paste(coefpre, from,"to",to,sep=""))
    name <- paste0(deg, "degrange")
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in c_*degrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste(coefpre, from,"+", ".homophily.",attrname,sep=""),
                         paste(coefpre, from,"to",to, ".homophily.",attrname,sep=""))
    name <- paste0(deg, "degrange_w_homophily")
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in c_*degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste(coefpre, du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste(coefpre, du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- paste0(deg, "degrange_by_attr")
    inputs <- c(as.vector(du), nodecov)
  }

  list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints=paste0(deg, "degreedist"), emptynwstats=emptynwstats)
}


.degree_impl <- function(deg, dir, bip, nw, arglist, ..., version){
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=dir, bipartite=bip,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=dir, bipartite=bip,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- a$levels
  }
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- if(NVL(bip, FALSE)) ergm_get_vattr(byarg, nw, bip = if(homophily) "n" else deg) else ergm_get_vattr(byarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in seq_along(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste0(deg, "degree", d)
    name <- paste0(deg, "degree")
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    # See comment in d_degree_w_homophily function
    coef.names <- paste0(deg, "deg", d, ".homophily.", attrname)
    name <- paste0(deg, "degree_w_homophily")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste0(deg, "deg", du[1,], ".", attrname,u[du[2,]])
    name <- paste0(deg, "degree_by_attr")
    inputs <- c(as.vector(du), nodecov)
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=network.size(nw), conflicts.constraints=paste0(deg, "degreedist"))
}

#=======================InitErgmTerm functions:  A============================#


################################################################################

#' @templateVar name absdiff
#' @title Absolute difference in nodal attribute
#' @description This term adds one network statistic to the model equaling
#'   the sum of `abs(attr[i]-attr[j])^pow` for all edges `(i,j)` in
#'   the network.
#'
#' @usage
#' # binary: absdiff(attr,
#' #                 pow=1)
#'
#' @template ergmTerm-attr
#' @param pow power to which to take the absolute difference
#'
#' @template ergmTerm-args-3.9.4
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept quantitative nodal attribute
InitErgmTerm.absdiff <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","pow"),
                        vartypes = c("character","numeric"),
                        defaultvalues = list(NULL,1),
                        required = c(TRUE,FALSE))
    ### Process the arguments
    nodecov <- get.node.attr(nw, a$attrname)
    covname <- a$attrname
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","pow"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric"),
                        defaultvalues = list(NULL,1),
                        required = c(TRUE,FALSE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
    covname <- attr(nodecov, "name")
  }
  ### Construct the list to return
  list(name="absdiff",                                     #name: required
       coef.names = paste(paste("absdiff",if(a$pow!=1) a$pow else "",sep=""), covname, sep="."), #coef.names: required
       inputs = c(a$pow,nodecov),  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}


################################################################################

#' @templateVar name absdiffcat
#' @title Categorical absolute difference in nodal attribute
#' @description This term adds one statistic for every possible nonzero distinct
#'	 value of `abs(attr[i]-attr[j])` in the network. The value of each such
#'	 statistic is the number of edges in the network with the corresponding
#'	 absolute difference. 
#'
#' @usage
#' # binary: absdiffcat(attr,
#' #                 base=NULL,
#' #                 levels=NULL)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain specifies which nonzero difference to include in or exclude from the model.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-args-3.9.4
#' @template ergmTerm-base-dep
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.absdiffcat <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","base"),
                        vartypes = c("character","numeric"),
                        defaultvalues = list(NULL,NULL),
                        required = c(TRUE,FALSE),
                        dep.inform = list(FALSE, "levels"))
    attrarg <- a$attrname
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","base","levels"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric",ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL,NULL,NULL),
                        required = c(TRUE,FALSE,FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr
  }
  ### Process the arguments
  nodecov <- ergm_get_vattr(attrarg, nw, accept = "numeric")
  attrname <- attr(nodecov, "name")

  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  
  u <- ergm_attr_levels(a$levels, nodecov, nw, levels = u)
  
  if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && any(NVL(a$base,0)!=0)) u <- u[-a$base]
  if (length(u)==0)
    ergm_Init_abort ("Argument to absdiffcat() has too few distinct differences")

  ### Construct the list to return
  inputs <- c(u, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(u) # See comment at top of file
  list(name="absdiffcat",                                  #name: required
       coef.names = paste("absdiff", attrname, u, sep="."), #coef.names: required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}


################################################################################

#' @templateVar name altkstar
#' @title Alternating \eqn{k}-star
#' @description Add one network statistic to the model equal to a weighted alternating
#'   sequence of \eqn{k}-star statistics with weight parameter `lambda`.
#' @details This is the version given in Snijders et al. (2006). The `gwdegree` and
#'   `altkstar` produce mathematically equivalent models, as long as they are used
#'   together with the `edges` (or `kstar(1)`) term, yet the interpretation of the
#'   `gwdegree` parameters is slightly more straightforward than the interpretation
#'   of the `altkstar` parameters. For this reason, we recommend the use of the
#'   `gwdegree` instead of `altkstar`. See Section 3 and especially equation (13)
#'   of Hunter (2007) for details.
#'
#' @usage
#' # binary: altkstar(lambda,
#' #                 fixed=FALSE)
#'
#' @param lambda weight parameter to model
#' @param fixed indicates whether the `decay` parameter is
#'   fixed at the given value, or is to be fit as a curved exponential family model
#'   (see Hunter and Handcock, 2006).  The default is `FALSE`, which means the scale
#'   parameter is not fixed and thus the model is a CEF model.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-undirected
#'
#' @concept curved
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.altkstar <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=NULL,
                      varnames = c("lambda","fixed"),
                      vartypes = c("numeric","logical"),
                      defaultvalues = list(NULL,FALSE),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  if(!a$fixed){ # This is a curved exponential family model
    if(!is.null(a$lambda)) warning("In term 'altkstar': decay parameter 'lambda' passed with 'fixed=FALSE'. 'lambda' will be ignored. To specify an initial value for 'lambda', use the 'init' control parameter.", call.=FALSE)
    d <- 1:(network.size(nw)-1)
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(x[2]*((1-1/x[2])^i + i) - 1)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(x[2]*((1-1/x[2])^i + i) - 1,
            x[1]*(i - 1 + (x[2]*x[2]-x[2]+i)*((1-1/x[2])^(i-1))/(x[2]*x[2])))
    }
    ### Construct the list to return
    outlist <- list(name="degree",                 #name: required
       coef.names = paste("altkstar#", d, sep=""), #coef.names: required
       inputs = d, map=map, gradient=gradient,
       params=list(altkstar=NULL, altkstar.lambda=a$lambda),
       minpar = c(-Inf, 0)
       )
  } else {
    if(is.null(a$lambda)) stop("Term 'altkstar' with 'fixed=TRUE' requires a decay parameter 'lambda'.", call.=FALSE)
    coef.names = paste("altkstar", a$lambda, sep=".")
    outlist <- list (name="altkstar",                      #name: required
                     coef.names = coef.names,
                     inputs=a$lambda
                     )
  }
  outlist
}



################################################################################

#' @templateVar name asymmetric
#' @title Asymmetric dyads
#' @description This term adds one network statistic to the model equal to the
#'   number of pairs of actors for which exactly one of
#'   \eqn{(i{\rightarrow}j)}{(i,j)} or \eqn{(j{\rightarrow}i)}{(j,i)} exists.
#'   
#' @usage
#' # binary: asymmetric(attr=NULL, diff=FALSE, keep=NULL, levels=NULL)
#'
#' @param attr quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If specified, only symmetric pairs that match on the vertex attribute are counted.
#'
#' @param diff Used in the same way as for the `nodematch` term. (See `nodematch` (`ergmTerm?nodematch`) for details.)
#' @param keep deprecated
#' @param level Used in the same way as for the `nodematch` term. (See `nodematch` (`ergmTerm?nodematch`) for details.)
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @template ergmTerm-keep-dep
#'
#' @concept directed
#' @concept dyad-independent
#' @concept triad-related

InitErgmTerm.asymmetric <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                        varnames = c("attrname", "diff", "keep"),
                        vartypes = c("character", "logical", "numeric"),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels"))
    attrarg <- a$attrname
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                        varnames = c("attr", "diff", "keep", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels", FALSE))
    attrarg <- a$attr
  }
  ### Process the arguments
  if (!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(a$levels, nodecov, nw, levels = sort(unique(nodecov)))
    if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && !is.null(a$keep)) u <- u[a$keep]
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }
  ### Construct the list to return
  out <- list(name="asymmetric",                      #name: required
              coef.names = "asymmetric",              #coef.names: required
              minval = 0,
              maxval = network.dyadcount(nw,FALSE)/2
              ) 
  if (!is.null(attrarg)) {
    if (a$diff) {
      out$coef.names <- paste("asymmetric", attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("asymmetric", attrname, sep=".")
      out$inputs <- nodecov
    }
  }

  out
}

################################################################################

#' @templateVar name attrcov
#' @title Edge covariate by attribute pairing
#' 
#' @description This term adds one statistic to the model, equal to the sum of the covariate values
#'   for each edge appearing in the network, where the covariate value for a given edge is determined by its mixing type on
#'   `attr`. Undirected networks are regarded as having undirected mixing, and it is assumed that `mat` is symmetric
#'   in that case.
#'   
#'   This term can be useful for simulating large networks with many mixing types, where `nodemix` would be slow due to
#'   the large number of statistics, and `edgecov` cannot be used because an adjacency matrix would be too big.
#'
#' @usage
#' # binary: attrcov(attr, mat)
#'
#' @template ergmTerm-attr
#'
#' @param mat a matrix of covariates with the same dimensions as a mixing matrix for `attr`
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmTerm.attrcov <- function (nw, arglist, ..., version=packageVersion("ergm")) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "mat"),
                      vartypes = c(ERGM_VATTR_SPEC, "matrix"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))


  if(is.bipartite(nw)) {
    b1nodecov <- ergm_get_vattr(a$attr, nw, bip="b1")
    b2nodecov <- ergm_get_vattr(a$attr, nw, bip="b2")

    attrname <- attr(b1nodecov, "name")
    
    b1levels <- sort(unique(b1nodecov))
    b2levels <- sort(unique(b2nodecov))
    
    nodecov <- c(match(b1nodecov, b1levels), match(b2nodecov, b2levels))
  
    if(NROW(a$mat) != length(b1levels) || NCOL(a$mat) != length(b2levels)) {
      ergm_Init_abort("mat has wrong dimensions for attr")
    }
  } else {
    nodecov <- ergm_get_vattr(a$attr, nw)
    attrname <- attr(nodecov, "name")
    
    levels <- sort(unique(nodecov))
    nodecov <- match(nodecov, levels)
      
    if(NROW(a$mat) != length(levels) || NCOL(a$mat) != length(levels)) {
      ergm_Init_abort("mat has wrong dimensions for attr")
    }
  }
  
  list(name = "attrcov", 
       coef.names = paste("attrcov", attrname, sep = "."),
       dependence = FALSE,
       inputs = NULL, # passed by name below
       nr = NROW(a$mat),
       nc = NCOL(a$mat),
       mat = if(is.double(a$mat)) a$mat else as.double(a$mat), # do not transpose, and only coerce if we have to, as a$mat can be quite large
       nodecov = c(0L, nodecov) - 1L # two shifts to make the C code cleaner
       )
}



#=======================InitErgmTerm functions:  B============================#

################################################################################

#' @templateVar name b1concurrent
#' @title Concurrent node count for the first mode in a bipartite network
#' @description This term adds one
#'   network statistic to the model, equal to the number of nodes in the first
#'   mode of the network with degree 2 or higher. The first mode of a bipartite
#'   network object is sometimes known as the "actor" mode. 
#'   This term can only be
#'   used with undirected bipartite networks.
#'
#' @usage
#' # binary: b1concurrent(by=NULL, levels=NULL)
#'
#' @param by optional argument specifying a vertex attribute (see Specifying
#'   Vertex attributes and Levels (`?nodal_attributes`) for details).
#'   It functions just like the `by` argument of the `b1degree` term.
#'   Without the optional argument, this statistic is equivalent to `b1mindegree(2)` .
#'
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute

InitErgmTerm.b1concurrent<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("by", "levels"),
                        vartypes = c("character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("by", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- a$levels    
  }

  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite") 
  byarg <- a$by
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw, bip = "b1")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_b1concurrent_by_attr function
    name <- "b1concurrent_by_attr"
    coef.names<-paste("b1concurrent",".", attrname, u, sep="")
    inputs <- c(ui, nodecov)
  }else{
    name <- "b1concurrent"
    coef.names<-paste("b1concurrent",sep="")
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval=0, maxval=nb1)
}

################################################################################

#' @templateVar name b1degrange
#' @title Degree range for the first mode in a bipartite network
#' @description This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i}th
#'   such statistic equals the number of nodes of the first mode
#'   ("actors") in the network of degree greater than or equal to
#'   `from[i]` but strictly less than `to[i]` , i.e. with edge count
#'   in semiopen interval `[from,to)` . 
#'   
#'   This term can only be used with bipartite networks; for directed networks
#'   see `idegrange` and `odegrange` . For undirected networks,
#'   see `degrange` , and see `b2degrange`
#'   for degrees of the second mode ("events").
#'
#' @usage
#' # binary: b1degrange(from, to=`+Inf`, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @template ergmTerm-from-to
#'
#' @template ergmTerm-by
#'
#' @template ergmTerm-general
#'
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @concept bipartite
#' @concept undirected
InitErgmTerm.b1degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degrange_impl("b1", NULL, TRUE, nw, arglist, ..., version=version)
}

################################################################################

#' @templateVar name b1cov
#' @title Main effect of a covariate for the first mode in a bipartite network
#' @description This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
#'   value of `attr(i)` for all edges /eqn{(i,j)} in the network. This
#'   term may only be used with bipartite networks. For categorical attributes,
#'   see `b1factor` .
#'   
#' @usage
#' # binary: b1cov(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept undirected
#' @concept bipartite
#' @concept dyad-independent
#' @concept quantitative nodal attribute
#' @concept frequently-used
InitErgmTerm.b1cov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                        varnames = c("attrname","transform","transformname"),
                        vartypes = c("character","function","character"),
                        defaultvalues = list(NULL,function(x)x,""),
                        required = c(TRUE,FALSE,FALSE))
    ### Process the arguments
    attrname<-a$attrname
    f<-a$transform
    f.name<-a$transformname
    coef.names <- paste(paste("b1cov",f.name,sep=""),attrname,sep=".")
    nb1 <- get.network.attribute(nw, "bipartite")
    nodecov <- f(get.node.attr(nw, attrname, "b1cov", numeric=TRUE)[1:nb1])
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                        varnames = c("attr"),
                        vartypes = c(ERGM_VATTR_SPEC),
                        defaultvalues = list(NULL),
                        required = c(TRUE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip = "b1", multiple="matrix")
    coef.names <- nodecov_names(nodecov, "b1cov")
  }
  # C implementation is identical
  list(name="nodeocov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}


################################################################################

#' @templateVar name b1degree
#' @title Degree for the first mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i}th such statistic equals the number of
#'   nodes of degree `d[i]` in the first mode of a bipartite network, i.e.
#'   with exactly `d[i]` edges. The first mode of a bipartite network object
#'   is sometimes known as the "actor" mode. 
#'   
#' @usage
#' # binary: b1degree(d, by=NULL, levels=NULL)
#'
#' @param d a vector of distinct integers. 
#'
#' @param by a vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details). If this is specified
#'   then each node's degree is tabulated only with other nodes having the same
#'   value of the `by` attribute.
#'
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b1degree <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degree_impl("b1", NULL, TRUE, nw, arglist, ..., version=version)
}


################################################################################

#' @templateVar name b1dsp
#' @title Dyadwise shared partners for dyads in the first bipartition
#' @description This term adds one
#'   network statistic to the model for each element in `d` ; the \eqn{i}th
#'   such statistic equals the number of dyads in the first bipartition with exactly
#'   `d[i]` shared partners. (Those shared partners, of course, must be members
#'   of the second bipartition.) This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: b1dsp(d)
#'
#' @param d a vector of distinct integers. 
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
InitErgmTerm.b1dsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d <- a$d
  
  if(length(d) == 0)
    return(NULL)

  emptynwstats <- rep(0, length(d))
  nb1 <- get.network.attribute(nw, "bipartite")

  type <- "OSP"
  typecode <- 4 # OSP type
  
  # in an empty network, the number of b1 dyads with zero shared
  # partners is just the number of b1 dyads, which is nb1*(nb1-1)/2
  emptynwstats[d==0] <- nb1*(nb1-1)/2 
  
  list(name="ddspbwrap", coef.names=paste("b1dsp",d,sep=""), inputs=c(typecode, d), 
       emptynwstats=emptynwstats, minval = 0, maxval = nb1*(nb1-1)/2, dependence = TRUE, auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
}


################################################################################

#' @templateVar name b1factor
#' @title Factor attribute effect for the first mode in a bipartite network
#' @description This term adds multiple network statistics to the model, one for each of (a subset of) the
#'   unique values of the `attr` attribute. Each of these statistics
#'   gives the number of times a node with that attribute in the first mode of
#'   the network appears in an edge. The first mode of a bipartite network object
#'   is sometimes known as the "actor" mode.
#'   
#' @usage
#' # binary: b1factor(attr, base=1, levels=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent
#' @concept frequently-used
#' @concept categorical nodal attribute
InitErgmTerm.b1factor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr", "base", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr                        
    levels <- a$levels    
  }

  ### Process the arguments  
  nodecov <- ergm_get_vattr(attrarg, nw, bip = "b1")
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

  if (attr(a,"missing")["levels"] && any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
  }

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  } 
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)-1

  ### Construct the list to return
  inputs <- nodepos
  list(name="nodeofactor", coef.names=paste("b1factor", attrname, paste(u), sep="."), inputs=inputs, dependence=FALSE, minval=0)
}

################################################################################

#' @templateVar name b1sociality
#' @title Degree
#' @description This term adds one network statistic for each node in the first bipartition, equal to the number of
#'   ties of that node. This term can only be used with bipartite networks. For directed networks, see `sender` and
#'   `receiver`. For unipartite networks, see `sociality`.
#'
#' @usage
#' # binary: b1sociality(nodes=-1)
#'
#' @param nodes By default, `nodes=-1` means that the statistic for the
#'   first node (in the second bipartition) will be omitted, but this argument may be changed to control
#'   which statistics are included. The `nodes` argument is interpreted using the new UI for level specification
#'   (see Specifying Vertex Attributes and Levels (`?nodal_attributes`) for details), where both the attribute and the sorted
#'   unique values are the vector of vertex indices `(nb1 + 1):n` , where
#'   `nb1` is the size of the first bipartition and `n` is the total number of nodes in the network. Thus `nodes=120` will include only the statistic
#'   for the 120th node in the second biparition, while `nodes=I(120)` will include only the statistic for the 120th node in the entire network.
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent

InitErgmTerm.b1sociality<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("nodes"),
                      vartypes = c(ERGM_LEVELS_SPEC),
                      defaultvalues = list(-1),
                      required = c(FALSE))
                      
  nb1 <- get.network.attribute(nw, "bipartite")
  d <- ergm_attr_levels(a$nodes, 1:nb1, nw, 1:nb1)
  
  ld<-length(d)
  if(ld==0){return(NULL)}

  coef.names <- paste("b1sociality",d,sep="")
  inputs <- c(d,0) # Input requires a "guard" value.

  list(name="sociality", coef.names=coef.names, inputs=inputs, minval=0, maxval=network.size(nw)-nb1, conflicts.constraints="b1degrees", dependence=FALSE)
}

################################################################################

#' @templateVar name b1star
#' @title \eqn{k}-stars for the first mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `k` . The \eqn{i} th such statistic counts the number of
#'   distinct `k[i]` -stars whose center node is in the first mode of the
#'   network. The first mode of a bipartite network object is sometimes known as
#'   the "actor" mode. A \eqn{k} -star is defined to be a center node \eqn{N} and
#'   a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the
#'   ties \eqn{\{N, O_i\}}{\{N, O[i]\}} exist for \eqn{i=1, \dots, k}. If `args` is specified then the count is over
#'   the number of \eqn{k}-stars (with center node in the first mode) where all
#'   nodes have the same value of the attribute. This term can only be used for
#'   undirected bipartite networks. 
#'
#' @usage
#' # binary: b1star(k, attr=NULL, levels=NULL)
#'
#' @param k a vector of distinct integers
#' @template ergmTerm-attr
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @note `b1star(1)` is equal to `b2star(1)` and to `edges` .
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b1star <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attrname", "levels"),
                         vartypes = c("numeric", "character", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attr", "levels"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))  
    attrarg <- a$attr
    levels <- a$levels  
  }
  ### Process the arguments
  if (!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "ostar"
    coef.names <- paste("b1star", a$k, ".", attrname, sep="")
    inputs <- c(a$k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- length(a$k)
  } 
  else {
    name <- "ostar"
    coef.names <- paste("b1star",a$k,sep="")
    inputs <- a$k
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval = 0)
}

################################################################################

#' @templateVar name b1starmix
#' @title Mixing matrix for \eqn{k}-stars centered on the first mode of a bipartite network
#' @description This term counts all \eqn{k}-stars in which
#'   the b2 nodes (called events in some contexts) are homophilous in the sense
#'   that they all share the same value of `attr` . However, the b1 node
#'   (in some contexts, the actor) at the center of the \eqn{k}-star does NOT have to
#'   have the same value as the b2 nodes; indeed, the values taken by the b1
#'   nodes may be completely distinct from those of the b2 nodes, which allows
#'   for the use of this term in cases where there are two separate nodal
#'   attributes, one for the b1 nodes and another for the b2 nodes (in this case,
#'   however, these two attributes should be combined to form a single nodal
#'   attribute, `attr`). A different statistic is created for each
#'   value of `attr` seen in a b1 node, even if no \eqn{k}-stars are observed
#'   with this value.
#'
#' @usage
#' # binary: b1starmix(k, attr, base=NULL, diff=TRUE)
#'
#' @param k only a single value of \eqn{k} is allowed
#' @template ergmTerm-attr
#' @param base deprecated
#'
#' @param diff whether a different statistic is created for each value seen in a b2 node. When `diff=TRUE`,
#'    the default, a different statistic is created for each value and thus the behavior of this term is reminiscent of the
#'   `nodemix` term, from which it takes its name; when `diff=FALSE` ,
#'   all homophilous \eqn{k}-stars are counted together, though these \eqn{k}-stars are still
#'   categorized according to the value of the central b1 node.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b1starmix <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attrname", "base", "diff"),
                         vartypes = c("numeric", "character", "numeric", "logical"),
                         defaultvalues = list(NULL, NULL, NULL, TRUE),
                         required = c(TRUE, TRUE, FALSE, FALSE))
    attrarg <- a$attrname
  } else {
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attr", "base", "diff"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, "numeric", "logical"),
                         defaultvalues = list(NULL, NULL, NULL, TRUE),
                         required = c(TRUE, TRUE, FALSE, FALSE))
    attrarg <- a$attr
  }
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- sort(unique(nodecov))
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { ergm_Init_abort("Only a single scalar k may be used with each b1starmix term") }
  b1namescov <- sort(unique(nodecov[1:nb1]))
  b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
  b1nodecov <- match(nodecov[1:nb1],b1namescov)
  b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
  namescov <- u[c(b1namescov, b2namescov)]
  nr <- length(b1namescov)
  nc <- length(b2namescov)
  nodecov <- c(b1nodecov, b2nodecov + nr)
  if (a$diff) {
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base,] }
    name <- "b1starmix"
    coef.names <- paste("b1starmix", a$k, attrname,
                        apply(matrix(namescov[u],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- 1:nr
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base] }
    name <- "b1starmixhomophily"
    coef.names <- paste("b1starmix", a$k, attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval = 0)
}

################################################################################

#' @templateVar name b1twostar
#' @title Two-star census for central nodes centered on the first mode of a bipartite network
#' @description This term takes two nodal attributes. Assuming that there are
#'   \eqn{n_1} values of `b1attr` among the b1 nodes and \eqn{n_2}
#'   values of `b2attr` among the b2 nodes, then the total number of
#'   distinct categories of two stars according to these two attributes is
#'   \eqn{n_1(n_2)(n_2+1)/2}. By default, this model term creates a distinct statistic
#'   counting each of these categories.
#'   
#' @usage
#' # binary: b1twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)
#'
#' @param b1attr b1 nodes (actors in some contexts) (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#' @param b2attr b2 nodes (events in some contexts). If `b2attr` is not passed, it is assumed to be the same as `b1attr` .
#' @param b1levels,b2levels,base,levels2 used to leave some of the categories out (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-base-dep2
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b1twostar <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("b1attrname", "b2attrname", "base", "b1levels", "b2levels"),
                         vartypes = c("character", "character", "numeric", "character,numeric,logical", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE, FALSE, FALSE),
                         dep.inform = list(FALSE, FALSE, "levels2", FALSE, FALSE))
    b1attrarg <- a$b1attrname
    b2attrarg <- a$b2attrname
    b1levels <- if(!is.null(a$b1levels)) I(a$b1levels) else NULL
    b2levels <- if(!is.null(a$b2levels)) I(a$b2levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("b1attr", "b2attr", "base", "b1levels", "b2levels", "levels2"),
                         vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL, NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                         dep.inform = list(FALSE, FALSE, "levels2", FALSE, FALSE, FALSE))
    b1attrarg <- a$b1attr
    b2attrarg <- a$b2attr
    b1levels <- a$b1levels
    b2levels <- a$b2levels
  }
  
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  
  b1nodecov <- ergm_get_vattr(b1attrarg, nw, bip = "b1")
  b1attrname <- attr(b1nodecov, "name")
  b1u <- ergm_attr_levels(b1levels, b1nodecov, nw, sort(unique(b1nodecov)))
  
  if(is.null(b2attrarg)) { b2attrarg <- b1attrarg }
  b2nodecov <- ergm_get_vattr(b2attrarg, nw, bip = "b2")
  b2attrname <- attr(b2nodecov, "name")  
  b2u <- ergm_attr_levels(b2levels, b2nodecov, nw, sort(unique(b2nodecov)))

  nr <- length(b1u)
  nc <- length(b2u)
  
  levels2.grid <- expand.grid(row = b1u, col = b2u, col2 = b2u, stringsAsFactors=FALSE)
  indices2.grid <- expand.grid(row = 1:nr, col = 1:nc, col2 = 1:nc)
  
  levels2.list <- transpose(levels2.grid[indices2.grid$col <= indices2.grid$col2,])
  indices2.grid <- indices2.grid[indices2.grid$col <= indices2.grid$col2,]
  
  levels2.sel <- if((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(a$base != 0)) levels2.list[-a$base]
                 else ergm_attr_levels(a$levels2, list(row = b1nodecov, col = b2nodecov, col2 = b2nodecov), nw, levels2.list)
  
  rows2keep <- match(levels2.sel,levels2.list, NA)
  rows2keep <- rows2keep[!is.na(rows2keep)]
  
  u <- indices2.grid[rows2keep,]
  
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  
  coef.names <- paste("b1twostar", b1attrname, b1u[u[,1]],  b2attrname,
                      apply(matrix(b2u[cbind(u[,2], u[,3])],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b1twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]), minval = 0)
}

################################################################################

#' @templateVar name b2concurrent
#' @title Concurrent node count for the second mode in a bipartite network
#' @description This term adds one
#'   network statistic to the model, equal to the number of nodes in the second
#'   mode of the network with degree 2 or higher. The second mode of a bipartite
#'   network object is sometimes known as the "event" mode. 
#'   Without the optional argument, this statistic is equivalent to `b2mindegree(2)`.
#'   
#' @usage
#' # binary: b2concurrent(by=NULL)
#'
#' @param by This optional argument specifie a vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details);
#'   it functions just like the `by` argument of the `b2degree` term.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept frequently-used
InitErgmTerm.b2concurrent<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("by", "levels"),
                        vartypes = c("character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("by", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- a$levels    
  }

  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  byarg <- a$by

  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw, bip = "b2")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_b2concurrent_by_attr function
    coef.names <- paste("b2concurrent",".", attrname,u, sep="")
    name <- "b2concurrent_by_attr"
    inputs <- c(ui, nodecov)
  }else{
    coef.names <- "b2concurrent"
    name <- "b2concurrent"
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw)-nb1)
}

################################################################################

#' @templateVar name b2cov
#' @title Main effect of a covariate for the second mode in a bipartite  network
#' @description This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
#'   value of `attr(j)` for all edges \eqn{(i,j)} in the network. This
#'   term may only be used with bipartite networks. For categorical attributes, see `b2factor`.
#'
#' @usage
#' # binary: b2cov(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept undirected
#' @concept bipartite
#' @concept dyad-independent
#' @concept quantitative nodal attribute
#' @concept frequently-used
InitErgmTerm.b2cov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname","transform","transformname"),
                        vartypes = c("character","function","character"),
                        defaultvalues = list(NULL,function(x)x,""),
                        required = c(TRUE,FALSE,FALSE))
    ### Process the arguments
    attrname<-a$attrname
    f<-a$transform
    f.name<-a$transformname
    coef.names <- paste(paste("b2cov",f.name,sep=""),attrname,sep=".")
    nb1 <- get.network.attribute(nw, "bipartite")
    nodecov <- f(get.node.attr(nw, attrname, "b2cov", numeric=TRUE)[(nb1+1):network.size(nw)])
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr"),
                        vartypes = c(ERGM_VATTR_SPEC),
                        defaultvalues = list(NULL),
                        required = c(TRUE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip = "b2", multiple="matrix")
    coef.names <- nodecov_names(nodecov, "b2cov")
  }
  list(name="b2cov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}


################################################################################

#' @templateVar name b2degrange
#' @title Degree range for the second mode in a bipartite network
#' @description This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i} th
#'   such statistic equals the number of nodes of the second mode
#'   ("events") in the network of degree greater than or equal to
#'   `from[i]` but strictly less than `to[i]` , i.e. with edge count
#'   in semiopen interval `[from,to)` . 
#'   
#'   This term can only be used with bipartite networks; for directed networks
#'   see `idegrange` and `odegrange` . For undirected networks,
#'   see `degrange` , and see `b1degrange`
#'   for degrees of the first mode ("actors").
#'
#' @usage
#' # binary: b2degrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @template ergmTerm-from-to
#'
#' @template ergmTerm-by
#'
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
InitErgmTerm.b2degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degrange_impl("b2", NULL, TRUE, nw, arglist, ..., version=version)
}



################################################################################

#' @templateVar name b2degree
#' @title Degree for the second mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `d` ; the \eqn{i} th such statistic equals the number of
#'   nodes of degree `d[i]` in the second mode of a bipartite network, i.e.
#'   with exactly `d[i]` edges. The second mode of a bipartite network
#'   object is sometimes known as the "event" mode. 
#'
#' @usage
#' # binary: b2degree(d, by=NULL)
#'
#' @param d a vector of distinct integers
#'
#' @param by this optional term specifies
#'   a vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details). If this is specified
#'   then each node's degree is tabulated only with other nodes having the same
#'   value of the `by` attribute.
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b2degree <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degree_impl("b2", NULL, TRUE, nw, arglist, ..., version=version)
}

################################################################################

#' @templateVar name b2dsp
#' @title Dyadwise shared partners for dyads in the second bipartition
#' @description This term adds one network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the number of dyads in the second bipartition with exactly
#'   `d[i]` shared partners. (Those shared partners, of course, must be members
#'   of the first bipartition.) This term can only be used with bipartite networks.
#'
#' @usage
#' # binary: b2dsp(d)
#'
#' @param d a vector of distinct integers
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
InitErgmTerm.b2dsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d <- a$d
  
  if(length(d) == 0)
    return(NULL)

  emptynwstats <- rep(0, length(d))
  nb2 <- network.size(nw) - get.network.attribute(nw, "bipartite")
  
  type <- "ISP"
  typecode <- 5 # ISP type
  
  # in an empty network, the number of b2 dyads with zero shared
  # partners is just the number of b2 dyads, which is nb2*(nb2-1)/2
  emptynwstats[d==0] <- nb2*(nb2-1)/2 
  
  list(name="ddspbwrap", coef.names=paste("b2dsp",d,sep=""), inputs=c(typecode, d), 
       emptynwstats=emptynwstats, minval = 0, maxval = nb2*(nb2-1)/2, dependence = TRUE, auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
}

################################################################################

#' @templateVar name b2factor
#' @title Factor attribute effect for the second mode in a bipartite network
#' @description This term adds multiple network statistics to the model, one for each of (a subset of) the
#'   unique values of the `attr` attribute. Each of these statistics
#'   gives the number of times a node with that attribute in the second mode of
#'   the network appears in an edge. The second mode of a bipartite network
#'   object is sometimes known as the "event" mode.
#'   
#' @usage
#' # binary: b2factor(attr, base=1, levels=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.b2factor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("attr", "base", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr                        
    levels <- a$levels    
  }

  ### Process the arguments
  nodecov <- ergm_get_vattr(attrarg, nw, bip = "b2")
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    if (attr(a,"missing")["levels"] && any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
  }

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  }
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)-1

  ### Construct the list to return
  inputs <- nodepos
  list(name="b2factor", coef.names=paste("b2factor", attrname, paste(u), sep="."), inputs=inputs, dependence=FALSE, minval=0)
}


################################################################################

#' @templateVar name b2sociality
#' @title Degree
#' @description This term adds one network statistic for each node in the second bipartition, equal to the number of
#'   ties of that node. For directed networks, see `sender` and
#'   `receiver` . For unipartite networks, see `sociality` .
#'
#' @usage
#' # binary: b2sociality(nodes=-1)
#'
#' @param nodes By default, `nodes=-1` means that the statistic for the
#'   first node (in the second bipartition) will be omitted, but this argument may be changed to control
#'   which statistics are included. The `nodes` argument is interpreted using the new UI for level specification
#'   (see Specifying Vertex Attributes and Levels (`?nodal_attributes`) for details), where both the attribute and the sorted
#'   unique values are the vector of vertex indices `(nb1 + 1):n` , where
#'   `nb1` is the size of the first bipartition and `n` is the total number of nodes in the network. Thus `nodes=120` will include only the statistic
#'   for the 120th node in the second biparition, while `nodes=I(120)` will include only the statistic for the 120th node in the entire network.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-bipartite
#'
#' @concept bipartite
#' @concept undirected
#' @concept dyad-independent
InitErgmTerm.b2sociality<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("nodes"),
                      vartypes = c(ERGM_LEVELS_SPEC),
                      defaultvalues = list(-1),
                      required = c(FALSE))
                      
  nb1 <- get.network.attribute(nw, "bipartite")
  d <- ergm_attr_levels(a$nodes, (1 + nb1):network.size(nw), nw, (1 + nb1):network.size(nw))
  
  ld<-length(d)
  if(ld==0){return(NULL)}

  coef.names <- paste("b2sociality",d,sep="")
  inputs <- c(d,0) # Input requires a "guard" value.
  
  list(name="sociality", coef.names=coef.names, inputs=inputs, minval=0, maxval=nb1, conflicts.constraints="b2degrees", dependence=FALSE)
}


################################################################################

#' @templateVar name b2star
#' @title \eqn{k}-stars for the second mode in a bipartite network
#' @description This term adds one network statistic to the model for
#'   each element in `k` . The \eqn{i} th such statistic counts the number of
#'   distinct `k[i]` -stars whose center node is in the second mode of the
#'   network. The second mode of a bipartite network object is sometimes known as
#'   the "event" mode. A \eqn{k} -star is defined to be a center node \eqn{N} and
#'   a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the
#'   ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k} . This term can only be used for
#'   undirected bipartite networks. 
#'
#' @usage
#' # binary: b2star(k, attr=NULL, levels=NULL)
#'
#' @param k a vector of distinct integers
#' @param attr quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) then the count is over
#'   the number of \eqn{k} -stars (with center node in the second mode) where all
#'   nodes have the same value of the attribute. 
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @note `b2star(1)` is equal to `b1star(1)` and to `edges` .
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b2star <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attrname", "levels"),
                         vartypes = c("numeric", "character", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attr", "levels"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))  
    attrarg <- a$attr
    levels <- a$levels  
  }
  ### Process the arguments
  if (!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "istar"
    coef.names <- paste("b2star", a$k, ".", attrname, sep="")
    inputs <- c(a$k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- length(a$k)
  } 
  else {
    name <- "istar"
    coef.names <- paste("b2star",a$k,sep="")
    inputs <- a$k
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval=0)
}

################################################################################

#' @templateVar name b2starmix
#' @title Mixing matrix for \eqn{k}-stars centered on the second mode of a bipartite network
#' @description This term is exactly the same as `b1starmix` except that the roles of
#'   b1 and b2 are reversed.
#'
#' @usage
#' # binary: b2starmix(k, attr, base=NULL, diff=TRUE)
#'
#' @param k only a single value of \eqn{k} is allowed
#' @template ergmTerm-attr
#' @param base deprecated
#' @param diff whether a different statistic is created for each value seen in a b1 node. When `diff=TRUE`,
#'    the default, a different statistic is created for each value and thus the behavior of this term is reminiscent of the
#'   `nodemix` term, from which it takes its name; when `diff=FALSE` ,
#'   all homophilous \eqn{k}-stars are counted together, though these \eqn{k}-stars are still
#'   categorized according to the value of the central b1 node.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b2starmix <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attrname", "base", "diff"),
                         vartypes = c("numeric", "character", "numeric", "logical"),
                         defaultvalues = list(NULL, NULL, NULL, TRUE),
                         required = c(TRUE, TRUE, FALSE, FALSE))
    attrarg <- a$attrname
  } else {
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("k", "attr", "base", "diff"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, "numeric", "logical"),
                         defaultvalues = list(NULL, NULL, NULL, TRUE),
                         required = c(TRUE, TRUE, FALSE, FALSE))
    attrarg <- a$attr
  }
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- sort(unique(nodecov))
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { ergm_Init_abort("Only a single scalar k may be used with each b2starmix term") }
  b1namescov <- sort(unique(nodecov[1:nb1]))
  b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
  b1nodecov <- match(nodecov[1:nb1],b1namescov)
  b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
  namescov <- u[c(b1namescov, b2namescov)]
  nr <- length(b1namescov)
  nc <- length(b2namescov)
  nodecov <- c(b1nodecov, b2nodecov + nr)
  if (a$diff) {
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base,] }
    name <- "b2starmix"
    coef.names <- paste("b2starmix", a$k, attrname,
                        apply(matrix(namescov[u[,2:1]],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- nr+(1:nc)
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base] }
    name <- "b2starmixhomophily"
    coef.names <- paste("b2starmix", a$k, attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval=0)
}

################################################################################

#' @templateVar name b2twostar
#' @title Two-star census for central nodes centered on the second mode of a bipartite network
#' @description This term is exactly the same as `b1twostar` except that the
#'   roles of b1 and b2 are reversed.
#'
#' @usage
#' # binary: b2twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)
#'
#' @param b1attr b1 nodes (actors in some contexts) (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#' @param b2attr b2 nodes (events in some contexts). If `b1attr` is not passed, it is assumed to be the same as `b2attr` .
#' @param b1levels,b2levels,base,levels2 used to leave some of the categories out (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-base-dep2
#'
#' @concept bipartite
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.b2twostar <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("b1attrname", "b2attrname", "base", "b1levels", "b2levels"),
                         vartypes = c("character", "character", "numeric", "character,numeric,logical", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels2", FALSE, FALSE))
    b1attrarg <- a$b1attrname
    b2attrarg <- a$b2attrname
    b1levels <- if(!is.null(a$b1levels)) I(a$b1levels) else NULL
    b2levels <- if(!is.null(a$b2levels)) I(a$b2levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("b1attr", "b2attr", "base", "b1levels", "b2levels", "levels2"),
                         vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL, NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels2", FALSE, FALSE, FALSE))
    b1attrarg <- a$b1attr
    b2attrarg <- a$b2attr
    b1levels <- a$b1levels
    b2levels <- a$b2levels
  }
  
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  
  b1nodecov <- ergm_get_vattr(b1attrarg, nw, bip = "b1")
  b1attrname <- attr(b1nodecov, "name")
  b1u <- ergm_attr_levels(b1levels, b1nodecov, nw, sort(unique(b1nodecov)))
  
  if(is.null(b2attrarg)) { b2attrarg <- b1attrarg }
  b2nodecov <- ergm_get_vattr(b2attrarg, nw, bip = "b2")
  b2attrname <- attr(b2nodecov, "name")  
  b2u <- ergm_attr_levels(b2levels, b2nodecov, nw, sort(unique(b2nodecov)))

  nr <- length(b1u)
  nc <- length(b2u)
  
  levels2.grid <- expand.grid(row = b2u, col = b1u, col2 = b1u, stringsAsFactors=FALSE)
  indices2.grid <- expand.grid(row = 1:nc, col = 1:nr, col2 = 1:nr)
  
  levels2.list <- transpose(levels2.grid[indices2.grid$col <= indices2.grid$col2,])
  indices2.grid <- indices2.grid[indices2.grid$col <= indices2.grid$col2,]
  
  levels2.sel <- if((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) levels2.list[-a$base]
                 else ergm_attr_levels(a$levels2, list(row = b2nodecov, col = b1nodecov, col2 = b1nodecov), nw, levels2.list)
  
  rows2keep <- match(levels2.sel,levels2.list, NA)
  rows2keep <- rows2keep[!is.na(rows2keep)]
  
  u <- indices2.grid[rows2keep,]
  
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  
  coef.names <- paste("b2twostar", b2attrname, b2u[u[,1]],  b1attrname,
                      apply(matrix(b1u[cbind(u[,2], u[,3])],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b2twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]), minval=0)
}

################################################################################

#' @templateVar name balance
#' @title Balanced triads
#' @description This term adds one network statistic to the model equal to the number of
#'   triads in the network that are balanced. The balanced triads are those of
#'   type `102` or `300` in the categorization of Davis and Leinhardt (1972). For details on the 16 possible triad types, see
#'   `?triad.classify` in the `{sna}` package. For an undirected
#'   network, the balanced triads are those with an odd number of ties (i.e., 1
#'   and 3).
#'
#' @usage
#' # binary: balance
#'
#' @template ergmTerm-general
#'
#' @concept triad-related
#' @concept directed
#' @concept undirected
InitErgmTerm.balance<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist)
    
  list(name="balance", coef.names="balance", dependence=TRUE, minval=0)
}


#=======================InitErgmTerm functions:  C============================#


################################################################################

#' @templateVar name concurrent
#' @title Concurrent node count
#' @description This term adds one network statistic to the model, equal to the number of
#'   nodes in the network with degree 2 or higher. 
#'   This term can only be used with undirected
#'   networks.
#'
#' @usage
#' # binary: concurrent(by=NULL, levels=NULL)
#' 
#' @param by this optional argument specifies a vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.)
#'   It functions just like the `by` argument of the `degree` term.
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.concurrent<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("by", "levels"),
                        vartypes = c("character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("by", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(FALSE, FALSE))
    levels <- a$levels    
  }
  byarg <- a$by
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)

  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_concurrent_by_attr function
    coef.names <- paste("concurrent",".", attrname,u, sep="")
    name <- "concurrent_by_attr"
    inputs <- c(ui, nodecov)
  }else{
    coef.names <- "concurrent"
    name <- "concurrent"
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw))
}




################################################################################

#' @templateVar name ctriple
#' @title Cyclic triples
#' @description By default, this term adds one
#'   statistic to the model, equal to the number of cyclic triples in the
#'   network, defined as a set of edges of the form \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k), (k{\rightarrow}i)\}}{\{(i,j), (j,k), (k,i)\}} . 
#'
#' @usage
#' # binary: ctriple(attr=NULL, diff=FALSE, levels=NULL)
#'
#' @param attr,diff quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If `attr` is specified and `diff` is `FALSE` , then the statistic is the number of cyclic triples where all
#'   three nodes have the same value of the attribute. If `attr` is specified and `diff` is `TRUE` , then one statistic is added to the model for each value of `attr`, equal to the number of cyclic triples where all
#'   three nodes have that value of the attribute.
#' @templateVar explain specifies the value of `attr` to consider if `attr` is passed and `diff=TRUE`.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @note for all directed networks, `triangle` is equal to
#'   `ttriple+ctriple` , so at most two of these three terms can be in a
#'   model. 
#'
#' @concept directed
#' @concept triad-related
#' @concept categorical nodal attribute
InitErgmTerm.ctriple<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attrname","diff", "levels"),
                        vartypes = c("character","logical", "character,numeric,logical"),
                        defaultvalues = list(NULL,FALSE,NULL),
                        required = c(FALSE,FALSE,FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                    
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attr","diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL,FALSE,NULL),
                        required = c(FALSE,FALSE,FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
  diff <- a$diff;
  if(!is.null(attrarg)){
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (!diff) {
      coef.names <- paste("ctriple",attrname,sep=".")
      inputs <- c(nodecov)
    } else {
      #  Number of input parameters before covariates equals number of
      #  unique elements in nodecov, namely length(u)
      coef.names <- paste("ctriple", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)
    }
  }else{
#    No attributes (or diff)
#    No covariates, so no need for "ParamsBeforeCov"
    coef.names <- "ctriple"
    inputs <- NULL
  }
  list(name="ctriple", coef.names=coef.names, inputs=inputs, minval = 0)
}

#' @templateVar name ctriple
#' @template ergmTerm-rdname
#' @aliases ctriad-ergmTerm
#' @usage
#' # binary: ctriad
InitErgmTerm.ctriad<-InitErgmTerm.ctriple

################################################################################

#' @templateVar name cycle
#' @title k-Cycle Census
#' @description This term adds one network statistic to the model for each value of `k` ,
#'   corresponding to the number of `k` -cycles (or, alternately, semicycles)
#'   in the graph.
#'   
#'   This term can be used with either directed or undirected networks.
#'
#' @usage
#' # binary: cycle(k, semi=FALSE)
#'
#' @param k a vector of integers giving the cycle lengths to count.
#'   Directed cycle lengths may range from `2` to `N` (the network size); undirected
#'   cycle lengths and semicycle lengths may range from `3` to `N` ; length 2 semicycles
#'   are not currently supported. 
#'   
#' @param semi an optional logical indicating whether semicycles
#'   (rather than directed cycles) should be counted; this is ignored in the
#'   undirected case.
#'   
#' @template ergmTerm-general
#'
#' @param directed 2-cycles are equivalent to mutual dyads.
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.cycle <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                     varnames = c("k","semi"),
                     vartypes = c("numeric","logical"),
                     defaultvalues = list(NULL,FALSE),
                     required = c(TRUE,FALSE))
  ### Process the arguments
  if(any(a$k > network.size(nw))) {
    ergm_Init_warn("cycles of length greater than the network size cannot exist and their statistics will be omitted")
    a$k <- a$k[a$k <= network.size(nw)]
  }

  if(!is.directed(nw) && any(a$k < 3)) {
    ergm_Init_warn("cycles of length less than 3 cannot exist in an undirected network and their statistics will be omitted")
    a$k <- a$k[a$k >= 3]  
  }

  if(any(a$k < 2)) {
    ergm_Init_warn("cycles of length less than 2 cannot exist and their statistics will be omitted")
    a$k <- a$k[a$k >= 2]
  }
  
  if(is.directed(nw) && a$semi && any(a$k == 2)) {
    ergm_Init_warn("semicycles of length 2 are not currently supported and their statistics will be omitted")
    a$k <- a$k[a$k >= 3]  
  }

  if (length(a$k)==0) return(NULL)

  semi<-is.directed(nw)&&a$semi  #Are we computing semicycles?
  ### Construct the list to return
  if(semi)
    basenam<-"semicycle"
  else
    basenam<-"cycle"
  list(name="cycle",                            #name: required
       coef.names = paste(basenam, a$k, sep=""),  #coef.names: required
       inputs = c(a$semi, max(a$k), (2:max(a$k)) %in% a$k),
       minval = 0)
}



#=======================InitErgmTerm functions:  D============================#


################################################################################

#' @templateVar name degcor
#' @title Degree Correlation
#' @description This term adds one network statistic equal to the correlation
#'   of the degrees of all pairs of nodes in the network which are tied.
#'   Only coded for undirected networks.
#'
#' @usage
#' # binary: degcor
#'
#' @template ergmTerm-general
#'
#' @concept undirected
InitErgmTerm.degcor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(nodes=TRUE))
  el=as.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="degcor",                            #name: required
       coef.names = "degcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}

################################################################################

#' @templateVar name degcrossprod
#' @title Degree Cross-Product
#' @description This term adds one network statistic equal to the mean of the cross-products
#'   of the degrees of all pairs of nodes in the network which are tied.
#'   Only coded for undirected networks.
#'
#' @usage
#' # binary: degcrossprod
#'
#' @template ergmTerm-general
#'
#' @concept undirected
InitErgmTerm.degcrossprod<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 
  ### Construct the list to return
  list(name="degcrossprod",                            #name: required
       coef.names = "degcrossprod",                    #coef.names: required
       inputs=2*summary(nw ~ edges),
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}

################################################################################


#' @templateVar name degrange
#' @title Degree range
#' @description This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i} th
#'   such statistic equals the number of nodes in the network of degree
#'   greater than or equal to
#'   `from[i]` but strictly less than `to[i]` , i.e. with edges
#'   in semiopen interval `[from,to)` .
#'
#' @details This term can only be used with undirected networks; for directed networks
#'   see `idegrange` and `odegrange` . This term can be used
#'   with bipartite networks, and will count nodes of both first and second mode in
#'   the specified degree range. To count only nodes of the first mode ("actors"), use `b1degrange`
#'   and to count only those fo the second mode ("events"), use `b2degrange` .
#'
#' @usage
#' # binary: degrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @template ergmTerm-from-to
#'
#' @template ergmTerm-by
#'
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degrange_impl("", FALSE, NULL, nw, arglist, ..., version=version)
}



################################################################################

#' @templateVar name degree
#' @title Degree
#' @description This term adds one
#'   network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the number of nodes in the network of degree
#'   `d[i]` , i.e. with exactly `d[i]` edges. 
#'   This term can only be used with undirected networks; for directed networks
#'   see `idegree` and `odegree` .
#'
#' @usage
#' # binary: degree(d, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @param d vector of distinct integers
#' @template ergmTerm-by
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.degree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degree_impl("", FALSE, NULL, nw, arglist, ..., version=version)
}


################################################################################

#' @templateVar name degree1.5
#' @title Degree to the 3/2 power
#' @description This term adds one network statistic to the model equaling the sum over
#'   the actors of each actor's degree taken to the 3/2 power (or,
#'   equivalently, multiplied by its square root). This term is an
#'   undirected analog to the terms of Snijders et al. (2010), equations
#'   (11) and (12). This term can only be used with undirected networks.
#'
#' @usage
#' # binary: degree1.5
#'
#' @template ergmTerm-general
#'
#' @concept undirected
InitErgmTerm.degree1.5<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="degreepopularity", coef.names="degree1.5",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="degreedist")
}


################################################################################
#' @include ergm-deprecated.R
#' @describeIn ergm-deprecated Use [`degree1.5`][degree1.5-ergmTerm] instead.
InitErgmTerm.degreepopularity<-function (nw, arglist, ...) {
  .Deprecated("degree1.5")
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="degreepopularity", coef.names="degreepopularity",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="degreedist")
}


################################################################################

#' @templateVar name density
#' @title Density
#' @description This term adds one network statistic equal to the density of the network.
#'   For undirected networks, `density` equals `kstar(1)` or
#'   `edges` divided by \eqn{n(n-1)/2} ; for directed networks,
#'   `density` equals `edges` or `istar(1)` or `ostar(1)`
#'   divided by \eqn{n(n-1)} .
#'
#' @usage
#' # binary: density
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmTerm.density<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="density", coef.names="density", dependence=FALSE, minval = 0, maxval = 1, conflicts.constraints="edges")
}

################################################################################

#' @templateVar name diff
#' @title Difference
#' @description For values of `pow` other than
#'   `0` , this term adds one network statistic to the model,
#'   equaling the sum, over directed edges \eqn{(i,j)} , of
#'   `sign.action(attr[i]-attr[j])^pow` if `dir` is
#'   `"t-h"` and of `sign.action(attr[j]-attr[i])^pow` if
#'   `"h-t"` . That is, the
#'   argument `dir` determines which vertex's attribute is
#'   subtracted from which, with tail being the origin of a directed edge
#'   and head being its destination, and bipartite networks' edges being
#'   treated as going from the first part (b1) to the second (b2).
#'   
#'   If `pow==0` , the exponentiation is replaced by the signum
#'   function: `+1` if the difference is positive, `0` if there
#'   is no difference, and `-1` if the difference is negative. Note
#'   that this function is applied after the
#'   `sign.action` . The comparison is exact, so when using
#'   calculated values of `attr` , ensure that values that you
#'   want to be considered equal are, in fact, equal.
#'   
#' @usage
#' # binary: diff(attr, pow=1, dir="t-h", sign.action="identity")
#'
#' @template ergmTerm-attr
#' @param pow exponent for the node difference
#' @param dir determines which vertix's attribute is subtracted from which. Accepts: `"t-h"` (the default), `"tail-head"` , `"b1-b2"`, `"h-t"` , `"head-tail"` , and `"b2-b1"` .
#' @param sign.action one of `"identity"`, `"abs"`, `"posonly"`, `"negonly"`. The following `sign.actions` are possible:
#'   
#'   - `"identity"` (the default) no transformation of the
#'   difference regardless of sign
#'
#'   - `"abs"` absolute value of the difference: equivalent
#'   to the absdiff term
#'
#'   - `"posonly"` positive differences are kept, negative
#'   differences are replaced by 0
#'
#'   - `"negonly"` negative differences are kept, positive
#'   differences are replaced by 0
#'
#' @template ergmTerm-general
#'
#' @note this term may not be meaningful for unipartite undirected
#'   networks unless `sign.action=="abs"` . When used on such a
#'   network, it behaves as if all edges were directed, going from the
#'   lower-indexed vertex to the higher-indexed vertex.
#'
#' @concept dyad-independent
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept bipartite
#' @concept quantitative nodal attribute
InitErgmTerm.diff <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attrname","pow", "dir", "sign.action"),
                        vartypes = c("character","numeric", "character", "character"),
                        defaultvalues = list(NULL,1, "t-h", "identity"),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr","pow", "dir", "sign.action"),
                        vartypes = c(ERGM_VATTR_SPEC,"numeric", "character", "character"),
                        defaultvalues = list(NULL,1, "t-h", "identity"),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
  }  
                        
  ### Process the arguments
  nodecov <- ergm_get_vattr(attrarg, nw, accept="numeric")
  attrname <- attr(nodecov, "name")
  DIRS <- c("t-h", "tail-head", "b1-b2",
            "h-t", "head-tail", "b2-b1")
  dir <- match.arg(tolower(a$dir), DIRS)
  dir.mul <- if(match(dir, DIRS)<=3) +1 else -1
  
  SIGN.ACTIONS <- c("identity", "abs", "posonly", "negonly")
  sign.action <- match.arg(tolower(a$sign.action), SIGN.ACTIONS)
  sign.code <- match(sign.action, SIGN.ACTIONS)

  if(sign.action!="abs" && !is.directed(nw) && !is.bipartite(nw)) ergm_Init_inform(paste("Note that behavior of term diff() on unipartite, undirected networks may be unexpected. See", sQuote("ergmTerm?diff"),"for more information."))
  
  # 1 and 4 are sign codes that allow negative differences.
  if(sign.code %in% c(1, 4) &&  a$pow!=round(a$pow)) ergm_Init_abort("In term diff(attr, pow, sign=",a$sign,"), pow must be an integer.")
  
  ### Construct the list to return
  list(name="diff",                                     #name: required
       coef.names = paste0("diff", if(a$pow!=1) a$pow else "", if(sign.action!="identity") paste0(".", sign.action), if(sign.action!="abs") paste0(".", dir), ".", attrname), #coef.names: required
       inputs = c(a$pow, dir.mul, sign.code, nodecov),  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

################################################################################

#' @templateVar name dsp
#' @title Dyadwise shared partners
#' @description This term adds one
#'   network statistic to the model for each element in `d` ; the \eqn{i} th
#'   such statistic equals the number of dyads in the network with exactly
#'   `d[i]` shared partners. This term can be used with directed and
#'   undirected networks.
#'   
#' @usage
#' # binary: dsp(d)
#'
#' @param d a vector of distinct integers
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn dsp
#' @templateVar kind (directed) dyad `(i,j)`
#' @templateVar see ddsp
#' @template ergmTerm-sp-to-dsp
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.dsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d <- a$d
  if (any(d==0)) {
    emptynwstats <- rep(0, length(d))
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw, "bipartite")
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw,FALSE)
    }
  }else{
    emptynwstats <- NULL
  }
  ld<-length(d)
  if(ld==0){return(NULL)}
  if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
  if (!is.null(emptynwstats)){
    list(name=dname, coef.names=paste("dsp",d,sep=""),
         inputs=c(d), emptynwstats=emptynwstats, minval = 0, auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }else{
    list(name=dname, coef.names=paste("dsp",d,sep=""),inputs=c(d), minval = 0, auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}


################################################################################

#' @templateVar name dyadcov
#' @title Dyadic covariate
#' @description This term adds three statistics to the model, each equal to the sum of the
#'   covariate values for all dyads occupying one of the three possible non-empty
#'   dyad states (mutual, upper-triangular asymmetric, and lower-triangular
#'   asymmetric dyads, respectively), with the empty or null state serving as a
#'   reference category. If the network is undirected, `x` is either a
#'   matrix of edgewise covariates, or a network; if the latter, optional
#'   argument `attrname` provides the name of the edge attribute to use for
#'   edge values. This term adds one statistic to the model, equal to the sum of
#'   the covariate values for each edge appearing in the network. The
#'   `edgecov` and `dyadcov` terms are equivalent for undirected
#'   networks.
#'
#' @usage
#' # binary: dyadcov(x, attrname=NULL)
#'
#' @template ergmTerm-x-attrname
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept quantitative dyadic attribute
InitErgmTerm.dyadcov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x","attrname"),
                      vartypes = c("matrix,network,character","character"),
                      defaultvalues = list(NULL,NULL),
                      required = c(TRUE,FALSE),
                      argexpr = substitute(arglist))

  l <- ergm_edgecov_args("dyadcov", nw, a); xm <- l$xm; cn <- l$cn

  if(is.directed(nw)){
   #Check for symmetry
   # DH:  Since nw is directed, why are we testing for symmetry here?  
   if (any(xm[upper.tri(xm)]!=t(xm)[upper.tri(xm)])){
     xm[lower.tri(xm)]<-t(xm)[lower.tri(xm)]
     ergm_Init_warn("asymmetric covariate in dyadcov; using upper triangle only")
   }
   coef.names <- paste(cn, c("mutual","utri","ltri"),sep=".")
  }else{
#  So it is undirected
    coef.names <- cn
  }
   
#  There is 1 input parameter before the covariate vector, so "ParamsBeforeCov"
#  is set to 1 (although in this case, this is actually arbitrary since
#  d_dyadcov ignores the value of inp->attrib).

   inputs = c(NCOL(xm), as.double(xm))
   attr(inputs, "ParamsBeforeCov") <- 1
   
   list(name = "dyadcov", coef.names=coef.names, inputs=inputs, dependence=FALSE)
}




#=======================InitErgmTerm functions:  E============================#


################################################################################

#' @templateVar name edgecov
#' @title Edge covariate
#' @description This term adds one statistic to the model, equal to the sum
#'   of the covariate values for each edge appearing in the network. The
#'   `edgecov` term applies to both directed and undirected networks. For
#'   undirected networks the covariates are also assumed to be undirected. The
#'   `edgecov` and `dyadcov` terms are equivalent for undirected
#'   networks.
#'
#' @usage
#' # binary: edgecov(x, attrname=NULL)
#'
#' @template ergmTerm-x-attrname
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept frequently-used
#' @concept quantitative dyadic attribute
InitErgmTerm.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network,character", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE),
                      argexpr = substitute(arglist))

  l <- ergm_edgecov_args("edgecov", nw, a); xm <- l$xm; cn <- l$cn

  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov", coef.names = cn, inputs = inputs, dependence=FALSE,
       minval = sum(c(xm)[c(xm)<0]),
       maxval = sum(c(xm)[c(xm)>0])
       )
}

################################################################################


#' @templateVar name edges
#' @title Number of edges in the network
#' @description This term adds one network statistic equal to the number of
#'   edges (i.e. nonzero values) in the network. For undirected networks, `edges`
#'   is equal to `kstar(1)`; for directed networks, edges is equal to both
#'   `ostar(1)` and `istar(1)`.
#'
#' @usage
#' # binary: edges
#' 
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmTerm.edges<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="edges", coef.names="edges", dependence=FALSE,
       minval = 0, maxval = network.dyadcount(nw,FALSE), conflicts.constraints="edges")
}



################################################################################

#' @templateVar name esp
#' @title Edgewise shared partners
#' @description This is just like the `dsp` term, except this term adds one network
#'   statistic to the model for each element in `d` where the \eqn{i} th such
#'   statistic equals the number of edges (rather than dyads) in the
#'   network with exactly `d[i]` shared partners. This term can be used with
#'   directed and undirected networks.
#'   
#' @usage
#' # binary: esp(d)
#'
#' @param d a vector of distinct integers
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn esp
#' @templateVar kind (directed) edge `i -> j`
#' @templateVar see desp
#' @template ergmTerm-sp-to-dsp
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.esp<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}

  list(name=dname, coef.names=paste("esp",d,sep=""), inputs=c(d), minval=0, auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
}



#=======================InitErgmTerm functions:  G============================#

################################################################################

#' @templateVar name gwb1degree
#' @title Geometrically weighted degree distribution for the first mode in a bipartite network
#' @description This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the `decay` parameter, which should be non-negative,
#'   for nodes in the first mode of a bipartite network. The first mode of a bipartite network
#'   object is sometimes known as the "actor" mode.
#'   
#'   This term can only be used with undirected bipartite
#'   networks.
#'
#' @usage
#' # binary: gwb1degree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)
#'
#' @templateVar multiplicand first mode degree frequencies
#' @template ergmTerm-gw-decay-fixed
#' @template ergmTerm-attr
#' @templateVar underlying degree
#' @template ergmTerm-gw-cutoff
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb1degree<-function(nw, arglist, gw.cutoff=30, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
    # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                        varnames = c("decay", "fixed", "attrname","cutoff", "levels"),
                        vartypes = c("numeric", "logical", "character","numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
    # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                        varnames = c("decay", "fixed", "attr","cutoff", "levels"),
                        vartypes = c("numeric", "logical", ERGM_VATTR_SPEC,"numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels
  }
  ### Process the arguments
  decay_vs_fixed(a, 'gwb1degree')
  decay<-a$decay; fixed<-a$fixed
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
  maxesp <- min(cutoff, network.size(nw)-nb1)

  d <- 1:maxesp
  if(!fixed){# This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="b1degree", coef.names=paste("gwb1degree#",d,sep=""), inputs=c(d),
           conflicts.constraints="b1degreedist", params=list(gwb1degree=NULL,gwb1degree.decay=decay)), GWDECAY)
  } else {
    if(!is.null(attrarg)) {
      nodecov <- ergm_get_vattr(attrarg, nw, bip="b1")
      attrname <- attr(nodecov, "name")
      u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      # See comment in d_gwb1degree_by_attr function
      name <- "gwb1degree_by_attr"
      coef.names <- paste("gwb1deg", decay, ".",attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwb1degree"
      coef.names <- paste("gwb1deg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    list(minval=0, maxval=network.size(nw), dependence=TRUE, name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="b1degreedist")
  }
}


################################################################################

#' @templateVar name gwb1dsp
#' @title Geometrically weighted dyadwise shared partner distribution for dyads in the first bipartition
#' @description This term adds one network statistic to the model equal to the geometrically
#'   weighted dyadwise shared partner distribution for dyads in the first bipartition with decay parameter
#'   `decay` parameter, which should be non-negative. This term can only be used with bipartite networks.
#'   
#' @usage
#' # binary: gwb1dsp(decay=0, fixed=FALSE, cutoff=30)
#'
#' @templateVar multiplicand shared partner counts
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying b1dsp
#' @template ergmTerm-gw-cutoff
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb1dsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("decay","fixed","cutoff"),
                      vartypes = c("numeric","logical","numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff),
                      required = c(FALSE, FALSE, FALSE))
  decay<-a$decay
  fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...

  type <- "OSP"
  typecode <- 4 # OSP type
  
  basenam <- "gwb1dsp"
  
  maxdsp <- min(cutoff, network.size(nw) - nw %n% "bipartite")

  if(!fixed){ # This is a curved exponential family model
    d <- 1:maxdsp

    if(length(d) == 0)
      return(NULL)
    
    # first name must match `basenam`
    params<-list(gwb1dsp=NULL,gwb1dsp.decay=decay)
    
    c(list(name="ddspbwrap", coef.names=paste("b1dsp#",d,sep=""), 
         inputs=c(if(!cache.sp) -1,typecode,d), params=params, auxiliaries=if(cache.sp) .spcache.aux(type) else NULL), GWDECAY)
  }else{
    coef.names <- paste("gwb1dsp.fixed",decay,sep=".")
    list(name="dgwdspbwrap", coef.names=coef.names, inputs=c(if(!cache.sp) -1,decay,typecode,maxdsp), auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
  }
}

################################################################################

#' @templateVar name gwb2degree
#' @title Geometrically weighted degree distribution for the second mode in a bipartite network
#' @description This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the which should be non-negative,
#'   for nodes in the
#'   second mode of a bipartite network. The second mode of a bipartite network
#'   object is sometimes known as the "event" mode.
#'   
#' @usage
#' # binary: gwb2degree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)
#'
#' @templateVar multiplicand second mode degree frequencies
#' @template ergmTerm-gw-decay-fixed
#' @template ergmTerm-attr
#' @templateVar underlying degree
#' @template ergmTerm-gw-cutoff
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb2degree<-function(nw, arglist, gw.cutoff=30, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
    # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                        varnames = c("decay", "fixed", "attrname","cutoff", "levels"),
                        vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
    # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                        varnames = c("decay", "fixed", "attr","cutoff", "levels"),
                        vartypes = c("numeric", "logical", ERGM_VATTR_SPEC,"numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels
  }
  
  ### Process the arguments  
  decay_vs_fixed(a, 'gwb2degree')
  decay<-a$decay; fixed<-a$fixed
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:nb1
  maxesp <- min(cutoff,nb1)
  d <- 1:maxesp
  if(!fixed){# This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="b2degree", coef.names=paste("gwb2degree#",d,sep=""), inputs=c(d),
           conflicts.constraints="b2degreedist", params=list(gwb2degree=NULL,gwb2degree.decay=decay)), GWDECAY)
  } else {
    if(!is.null(attrarg)) {
      nodecov <- ergm_get_vattr(attrarg, nw, bip="b2")
      attrname <- attr(nodecov, "name")
      u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
     #  No covariates here, so "ParamsBeforeCov" unnecessary
     # See comment in d_gwb2degree_by_attr function
      name <- "gwb2degree_by_attr"
      coef.names <- paste("gwb2deg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwb2degree"
      coef.names <- paste("gwb2deg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    list(minval=0, maxval=network.size(nw), dependence=TRUE, name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="b2degreedist")
  }
}

################################################################################

#' @templateVar name gwb2dsp
#' @title Geometrically weighted dyadwise shared partner distribution for dyads in the second bipartition
#' @description This term adds one network statistic to the model equal to the geometrically
#'   weighted dyadwise shared partner distribution for dyads in the second bipartition with decay parameter
#'   `decay` parameter, which should be non-negative. This term can only be used with bipartite networks.
#'   
#' @usage
#' # binary: gwb2dsp(decay=0, fixed=FALSE, cutoff=30)
#'
#' @templateVar multiplicand shared partner counts
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying b2dsp
#' @template ergmTerm-gw-cutoff
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @concept bipartite
#' @concept undirected
#' @concept curved
InitErgmTerm.gwb2dsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("decay","fixed","cutoff"),
                      vartypes = c("numeric","logical","numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff),
                      required = c(FALSE, FALSE, FALSE))
  decay<-a$decay
  fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...

  type <- "ISP"
  typecode <- 5 # ISP type
  
  basenam <- "gwb2dsp"
  
  maxdsp <- min(cutoff, nw %n% "bipartite")

  if(!fixed){ # This is a curved exponential family model
    d <- 1:maxdsp

    if(length(d) == 0)
      return(NULL)
    
    # first name must match `basenam`
    params<-list(gwb2dsp=NULL,gwb2dsp.decay=decay)
    
    c(list(name="ddspbwrap", coef.names=paste("b2dsp#",d,sep=""), 
         inputs=c(if(!cache.sp) -1,typecode,d), params=params, auxiliaries=if(cache.sp) .spcache.aux(type) else NULL), GWDECAY)
  }else{
    coef.names <- paste("gwb2dsp.fixed",decay,sep=".")    
    list(name="dgwdspbwrap", coef.names=coef.names, inputs=c(if(!cache.sp) -1,decay,typecode,maxdsp), auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
  }
}

################################################################################

#' @templateVar name gwdegree
#' @title Geometrically weighted degree distribution
#' @description This term adds one network statistic to the model equal to the weighted
#'   degree distribution with decay controlled by the `decay` parameter, which should be non-negative.
#'   
#' @usage
#' # binary: gwdegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)
#'
#' @templateVar multiplicand degree frequencies
#' @template ergmTerm-gw-decay-fixed
#' @template ergmTerm-attr
#' @templateVar underlying degree
#' @template ergmTerm-gw-cutoff
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept curved
#' @concept frequently-used
InitErgmTerm.gwdegree<-function(nw, arglist, gw.cutoff=30, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("decay", "fixed", "attrname","cutoff", "levels"),
                        vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("decay", "fixed", "attr","cutoff", "levels"),
                        vartypes = c("numeric", "logical", ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }

  decay_vs_fixed(a, 'gwdegree')
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if(!fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="degree", coef.names=paste("gwdegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="degreedist", params=list(gwdegree=NULL,gwdegree.decay=decay)), GWDECAY)
  } else {
    if(!is.null(attrarg)) {
      nodecov <- ergm_get_vattr(attrarg, nw)
      attrname <- attr(nodecov, "name")
      u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwdegree_by_attr"
      coef.names <- paste("gwdeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwdegree"
      coef.names <- paste("gwdeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    list(minval=0, maxval=network.size(nw), dependence=TRUE, name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="degreedist")
  }
}


################################################################################

#' @templateVar name gwdsp
#' @title Geometrically weighted dyadwise shared partner distribution
#' @description This term adds one network statistic to the model equal to the geometrically
#'   weighted dyadwise shared partner distribution with decay parameter
#'   `decay` parameter, which should be non-negative. 
#'   This term can be used with directed and undirected networks.
#'   
#' @usage
#' # binary: gwdsp(decay, fixed=FALSE, cutoff=30)
#'
#' @templateVar multiplicand shared partner or directed 2-path count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying DSP
#' @template ergmTerm-gw-cutoff
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn gwdsp
#' @templateVar kind (directed) dyad `(i,j)`
#' @templateVar see dgwdsp
#' @template ergmTerm-sp-to-dsp
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept directed
#' @concept undirected
#' @concept curved
InitErgmTerm.gwdsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","alpha"),
                      vartypes = c("numeric","logical","numeric","numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  decay_vs_fixed(a, 'gwdsp')
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  if(!fixed){ # This is a curved exponential family model
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
    c(list(name=dname, coef.names=paste("gwdsp#",d,sep=""), 
           inputs=c(d), params=list(gwdsp=NULL,gwdsp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    if (!fixed)
      coef.names <- "gwdsp"
    else  # fixed == TRUE
      coef.names <- paste("gwdsp.fixed.",decay,sep="")
  if(is.directed(nw)){dname <- "gwtdsp"}else{dname <- "gwdsp"}
  list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}



################################################################################

#' @templateVar name gwesp
#' @title Geometrically weighted edgewise shared partner distribution
#' @description This term is just like `gwdsp` except it adds a statistic equal to the
#'   geometrically weighted edgewise (not dyadwise) shared partner
#'   distribution with decay parameter
#'   `decay` parameter, which should be non-negative. This term can be used with directed and
#'   undirected networks.
#'   
#' @usage
#' # binary: gwesp(decay, fixed=FALSE, cutoff=30)
#'
#' @templateVar multiplicand shared partner or directed 2-path count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying ESP
#' @template ergmTerm-gw-cutoff
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn gwesp
#' @templateVar kind (directed) edge `i -> j`
#' @templateVar see dgwesp
#' @template ergmTerm-sp-to-dsp
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept curved
InitErgmTerm.gwesp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  decay_vs_fixed(a, 'gwesp')
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!fixed){ # This is a curved exponential family model
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
    c(list(name=dname, coef.names=paste("esp#",d,sep=""), 
         inputs=c(d), params=list(gwesp=NULL,gwesp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    coef.names <- paste("gwesp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtesp"}else{dname <- "gwesp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}



################################################################################

#' @templateVar name gwidegree
#' @title Geometrically weighted in-degree distribution
#' @description This term adds one network statistic to the model
#'   equal to the weighted in-degree distribution with decay parameter
#'   `decay` parameter, which should be non-negative. This
#'   term can only be used with directed networks.
#'   
#' @usage
#' # binary: gwidegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)
#'
#' @templateVar multiplicand indegree frequencies
#' @template ergmTerm-gw-decay-fixed
#' @template ergmTerm-attr
#' @templateVar underlying degree
#' @template ergmTerm-gw-cutoff
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept curved
InitErgmTerm.gwidegree<-function(nw, arglist, gw.cutoff=30, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("decay", "fixed", "attrname","cutoff", "levels"),
                        vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("decay", "fixed", "attr","cutoff", "levels"),
                        vartypes = c("numeric", "logical", ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
  decay_vs_fixed(a, 'gwidegree')
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
  maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if(!fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="idegree", coef.names=paste("gwidegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="idegreedist", params=list(gwidegree=NULL,gwidegree.decay=decay)), GWDECAY)
  } else { 
    if(!is.null(attrarg)) {
      nodecov <- ergm_get_vattr(attrarg, nw)
      attrname <- attr(nodecov, "name")
      u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwidegree_by_attr"
      coef.names <- paste("gwideg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwidegree"
      coef.names <- paste("gwideg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    list(minval=0, maxval=network.size(nw), dependence=TRUE, name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="idegreedist")
  }
}


################################################################################

#' @templateVar name gwnsp
#' @title Geometrically weighted nonedgewise shared partner distribution
#' @description This term is just like
#'   `gwesp` and `gwdsp` except it adds a statistic equal to
#'   the geometrically weighted *nonedgewise* (that is, over dyads
#'   that do not have an edge) shared partner distribution with weight
#'   parameter `decay` parameter, which should be non-negative. This term can be used with
#'   directed and undirected networks.
#'   
#' @usage
#' # binary: gwnsp(decay, fixed=FALSE, cutoff=30)
#'
#' @templateVar multiplicand shared partner or directed 2-path count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying NSP
#' @template ergmTerm-gw-cutoff
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn gwnsp
#' @templateVar kind (directed) non-edge `(i,j)`
#' @templateVar see dgwnsp
#' @template ergmTerm-sp-to-dsp
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept directed
#' @concept undirected
#' @concept curved
InitErgmTerm.gwnsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  decay_vs_fixed(a, 'gwnsp')
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!fixed){ # This is a curved exponential family model
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tnsp"}else{dname <- "nsp"}
    c(list(name=dname, coef.names=paste("nsp#",d,sep=""),
           inputs=c(d), params=list(gwnsp=NULL,gwnsp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    coef.names <- paste("gwnsp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtnsp"}else{dname <- "gwnsp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)    
  }
}


################################################################################

#' @templateVar name gwodegree
#' @title Geometrically weighted out-degree distribution
#' @description This term adds one network statistic to the model
#'   equal to the weighted out-degree distribution with decay parameter
#'   `decay` parameter, which should be non-negative. This
#'   term can only be used with directed networks.
#'   
#' @usage
#' # binary: gwodegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)
#'
#' @templateVar multiplicand outdegree frequencies
#' @template ergmTerm-gw-decay-fixed
#' @template ergmTerm-attr
#' @templateVar underlying degree
#' @template ergmTerm-gw-cutoff
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept curved
InitErgmTerm.gwodegree<-function(nw, arglist, gw.cutoff=30, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("decay", "fixed", "attrname","cutoff", "levels"),
                        vartypes = c("numeric", "logical", "character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("decay", "fixed", "attr","cutoff", "levels"),
                        vartypes = c("numeric", "logical", ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, gw.cutoff, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
  decay_vs_fixed(a, 'gwodegree')
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
  maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if(!fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="odegree", coef.names=paste("gwodegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="odegreedist", params=list(gwodegree=NULL,gwodegree.decay=decay)), GWDECAY)
  } else {
    if(!is.null(attrarg)) {
      nodecov <- ergm_get_vattr(attrarg, nw)
      attrname <- attr(nodecov, "name")
      u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwodegree_by_attr"
      coef.names <- paste("gwodeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwodegree"
      coef.names <- paste("gwodeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    list(minval=0, maxval=network.size(nw), dependence=TRUE, name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="odegreedist")
  }
}



#=======================InitErgmTerm functions:  H============================#


################################################################################

#' @templateVar name hamming
#' @title Hamming distance
#' @description This term adds one statistic to the model equal to the weighted or
#'   unweighted Hamming distance of the network from the network specified by
#'   `x` . Unweighted Hamming distance is defined as the total
#'   number of pairs \eqn{(i,j)} (ordered or unordered, depending on whether the
#'   network is directed or undirected) on which the two networks differ. If the
#'   optional argument `cov` is specified, then the weighted Hamming
#'   distance is computed instead, where each pair \eqn{(i,j)} contributes a
#'   pre-specified weight toward the distance when the two networks differ on
#'   that pair.
#'
#' @usage
#' # binary: hamming(x, cov, attrname=NULL)
#'
#' @param x defaults to be the observed
#'   network, i.e., the network on the left side of the \eqn{\sim} in the formula
#'   that defines the ERGM.
#' @param cov either a matrix of edgewise weights or a network
#' @param attrname option argument that provides the name of the edge attribute
#'   to use for weight values when a network is specified in `cov`
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmTerm.hamming<-function (nw, arglist, ...) {
  a <- check.ErgmTerm (nw, arglist,
        varnames = c("x","cov","attrname","defaultweight"),
        vartypes = c("matrix,network","matrix,network","character","numeric"),
        defaultvalues = list(nw, NULL, NULL, NULL),
        required = c(FALSE, FALSE, FALSE, FALSE))

  ## Process hamming network ##
  if(is.network(a$x)){                                                    # Arg to hamming is a network
    # check for attribute existance before creating matrix
  
    if( is.null(a$attrname) || is.null(get.edge.attribute(a$x,a$attrname))){ 
      xm<-as.edgelist(a$x)  # so call the non attribute version
    } else {
      xm<-as.edgelist(a$x,a$attrname)
    }
    
    
  }else if(is.character(a$x)){                                                # Arg to hamming is the name of an attribute in nw
    xm<-get.network.attribute(nw,a$x)
    xm<-as.edgelist(xm)
  }else if(is.null(a$x)){
    xm<-as.edgelist(nw)                                # Arg to hamming does not exist; uses nw
  }else if(is.matrix(a$x) && ncol(a$x)!=2){
    xm<-as.edgelist(update(nw,a$x,matrix.type="adjacency"))
  }else{
    xm<-as.matrix(a$x)                                                    # Arg to hamming is anything else; attempts to coerce
  }
  if (is.vector(xm)) xm <- matrix(xm, ncol=2)

  ## Process case without dyadcov (i.e. unweighted) ##
  sc03 <- sys.call(0)[[3]]
  coef.names <- "hamming"  # This might be modified later
  if (is.null(a$cov)) {
    minval <- 0
    maxval <- network.dyadcount(nw,FALSE)
    if (length(sc03)>1) 
      coef.names <- paste("hamming", as.character(sc03[[2]]), sep=".")
    covm <- NULL
    if (is.null(a$defaultweight))
      a$defaultweight <- 1.0
    emptynwstats <- NROW(xm) * a$defaultweight

  ## Process case with dyadcov (i.e. weighted) ##
  } else {
    # Extract dyadic covariate
    if(is.network(a$cov)){
      covm<-as.edgelist(a$cov,a$attrname)
      if(length(covm)==2){covm <- matrix(covm,ncol=2)}
      if(length(covm)==3){covm <- matrix(covm,ncol=3)}
      if (NCOL(covm)==2)
        covm <- cbind(covm,1)
    }else if(is.character(a$cov)){
      covm<-get.network.attribute(nw,a$cov)
      covm<-as.edgelist(covm) # DH:  Not really sure what should happen here
    }else{
      covm<-as.matrix(a$cov)
    }
    if (is.null(covm) || !is.matrix(covm) || NCOL(covm)!=3){
      ergm_Init_abort("Improper dyadic covariate passed to hamming()")
    }
    emptynwstats <- sum(apply(xm, 1, function(a,b) sum(b[(a[1]==b[,1] & a[2]==b[,2]),3]), covm))
    if (is.null(a$defaultweight))
      a$defaultweight <- 0
    if(!is.null(a$attrname) && length(sc03)>1){
      coef.names<-paste("hamming", as.character(sc03[2]), "wt",
                as.character(a$attrname), sep = ".")
    }else if (length(sc03)>1) {
      coef.names<-paste("hamming", as.character(sc03[2]), "wt",
                as.character(sys.call(0)[[3]][3]), sep = ".")
    }
    minval <- sum(c(covm)[c(covm)<0])
    maxval <- sum(c(covm)[c(covm)<0])
  }
  ## Return ##
  if (!is.null(xm)) {
    xm <- to_ergm_Cdouble(xm, prototype=nw)
  }
  if (!is.null(covm)) {
    covm <- to_ergm_Cdouble(covm, prototype=nw)
  }else covm <- 0
  inputs <- c(xm, a$defaultweight, covm)
  list(name="hamming", coef.names=coef.names, #name and coef.names: required 
       inputs = inputs, emptynwstats = emptynwstats, dependence = FALSE,
       minval = minval, maxval = maxval)
}

################################################################################
#' @rdname ergm-deprecated
#' @aliases hammingmix
InitErgmTerm.hammingmix<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  .Deprecate_once(msg="hammingmix() has been deprecated due to disuse.")
  if(version <= as.package_version("3.9.4")){
    # There is no reason hammingmix should be directed-only, but for now
    # the undirected version does not seem to work properly, so:
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attrname","x","base","contrast"),
                        vartypes = c("character","matrix,network","numeric","logical"),
                        defaultvalues = list(NULL,nw,NULL,FALSE),
                        required = c(TRUE,FALSE,FALSE,FALSE),
                        dep.inform = list(FALSE, FALSE, "levels2", FALSE))
    attrarg <- a$attrname
  }else{
    # There is no reason hammingmix should be directed-only, but for now
    # the undirected version does not seem to work properly, so:
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attr", "x", "base", "levels", "levels2","contrast"),
                        vartypes = c(ERGM_VATTR_SPEC, "matrix,network", "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC,"logical"),
                        defaultvalues = list(NULL,nw,NULL,NULL,NULL,FALSE),
                        required = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE),
                        dep.inform = list(FALSE, FALSE, "levels2", FALSE, FALSE, FALSE))
    attrarg <- a$attr
  }

  x<-a$x

  if (a$contrast) {
    ergm_Init_abort("The 'contrast' argument of the hammingmix term is deprecated.  Use 'levels2' instead")
  }
  if(is.network(x)){
    xm<-as.edgelist(x)
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.edgelist(xm)
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    ergm_Init_abort("hammingmix() requires an edgelist")
  }

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")

  u <- ergm_attr_levels(a$levels, nodecov, nw, sort(unique(nodecov)))
  namescov <- u
  
  nr <- length(u)
  nc <- length(u)

  levels2.list <- transpose(expand.grid(row = u, col = u, stringsAsFactors=FALSE))
  indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
    
  levels2.sel <- if((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) levels2.list[-a$base]
                 else ergm_attr_levels(a$levels2, list(row = nodecov, col = nodecov), nw, levels2.list)
  
  rows2keep <- match(levels2.sel,levels2.list, NA)
  rows2keep <- rows2keep[!is.na(rows2keep)]
  
  u <- indices2.grid[rows2keep,]

  nodecov.indices <- match(nodecov, namescov, nomatch=length(namescov) + 1)

  coef.names <- paste("hammingmix",attrname,
                      apply(matrix(namescov[as.matrix(u)],ncol=2),1,paste,collapse="."), 
                      sep=".")
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui),
  inputs=c(to_ergm_Cdouble(xm, prototype=nw), u[,1], u[,2], nodecov.indices)
  attr(inputs, "ParamsBeforeCov") <- nrow(u)
  # The emptynwstats code below does not work right for
  # undirected networks, mostly since hammingmix doesn't work 
  # in this case anyway.
  nw %v% "_tmp_nodecov" <- as.vector(nodecov)
  if(version <= as.package_version("3.9.4")){
    emptynwstats <- summary(nw ~ nodemix("_tmp_nodecov", base=a$base))
  }else{
    nodemix.call <- c(list(as.name("nodemix"),"_tmp_nodecov"), list(base=a$base, levels=a$levels, levels2=a$levels2)[!attr(a,"missing")[c("base","levels","levels2")]])
    nodemix.call <- as.call(nodemix.call)
    nodemix.form <- as.formula(call("~", nw, nodemix.call))
    emptynwstats <- summary(nodemix.form)
  }
  list(name="hammingmix", coef.names=coef.names, inputs=inputs, 
       emptynwstats=emptynwstats, dependence=FALSE)
}






#=======================InitErgmTerm functions:  I============================#

################################################################################

#' @templateVar name idegrange
#' @title In-degree range
#' @description This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i} th
#'   such statistic equals the number of nodes in the network of in-degree
#'   greater than or equal to `from[i]` but strictly less than `to[i]` , i.e. with
#'   in-edge count in semiopen interval `[from,to)` .
#'   
#'   This term can only be used with directed networks; for undirected
#'   networks (bipartite and not)
#'   see `degrange` . For degrees of specific modes of bipartite
#'   networks, see `b1degrange` and `b2degrange` . For
#'   in-degrees, see `idegrange` .
#'
#' @usage
#' # binary: idegrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @template ergmTerm-from-to
#' @template ergmTerm-by
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.idegrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degrange_impl("i", TRUE, NULL, nw, arglist, ..., version=version)
}

################################################################################

#' @templateVar name idegree
#' @title In-degree
#' @description This term adds one network statistic to
#'   the model for each element in `d` ; the \eqn{i} th such statistic equals
#'   the number of nodes in the network of in-degree `d[i]` , i.e. the number
#'   of nodes with exactly `d[i]` in-edges. 
#'   This term can only be used with directed networks; for undirected networks
#'   see `degree` .
#'
#' @usage
#' # binary: idegree(d, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-by
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.idegree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degree_impl("i", TRUE, NULL, nw, arglist, ..., version=version)
}



################################################################################

#' @templateVar name idegree1.5
#' @title In-degree to the 3/2 power
#' @description This term adds one network statistic to the model equaling the sum
#'   over the actors of each actor's indegree taken to the 3/2 power
#'   (or, equivalently, multiplied by its square root). This term is
#'   analogous to the term of Snijders et al. (2010), equation (12). This
#'   term can only be used with directed networks.
#'
#' @usage
#' # binary: idegree1.5
#'
#' @template ergmTerm-general
#'
#' @concept directed
InitErgmTerm.idegree1.5<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="idegreepopularity", coef.names="idegree1.5",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="idegreedist")
}


################################################################################
#' @describeIn ergm-deprecated Use [`idegree1.5`][idegree1.5-ergmTerm] instead.
InitErgmTerm.idegreepopularity<-function (nw, arglist, ...) {
  .Deprecated("idegree1.5")
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="idegreepopularity", coef.names="idegreepopularity",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="idegreedist")
}



################################################################################

#' @templateVar name intransitive
#' @title Intransitive triads
#' @description This term adds one statistic to the model, equal to the number of triads in
#'   the network that are intransitive. The intransitive triads are those of type
#'   `111D` , `201` , `111U` , `021C` , or `030C` in the
#'   categorization of Davis and Leinhardt (1972). For details on the 16 possible
#'   triad types, see `triad.classify` in the
#'   \CRANpkg{sna} package. Note the distinction from the `ctriple`
#'   term.
#'
#' @usage
#' # binary: intransitive
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept directed
#' @concept triad-related
InitErgmTerm.intransitive<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="intransitive", coef.names="intransitive", minval = 0)
}

################################################################################

#' @templateVar name isolatededges
#' @title Isolated edges
#' @description This term adds one statistic to the
#'   model equal to the number of isolated edges in the network, i.e., the number
#'   of edges each of whose endpoints has degree 1. This term can only be used
#'   with undirected networks.
#'
#' @usage
#' # binary: isolatededges
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept bipartite
InitErgmTerm.isolatededges <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=NULL,
                     varnames = NULL,
                     vartypes = NULL,
                     defaultvalues = list(),
                     required = NULL)
  ### Construct the list to return
  list(name="isolatededges",                               #name: required
       coef.names = "isolatededges",                       #coef.names: required
       emptynwstats = 0,                                   #When nw is empty, isolatededges=0
       minval = 0,
       maxval = if(is.bipartite(nw)) min(nw%n%"bipartite", network.size(nw) - nw%n%"bipartite") else floor(network.size(nw)/2),
       dependence = TRUE
       )                                                               
}

################################################################################

#' @templateVar name isolates
#' @title Isolates
#' @description This term adds one statistic to the
#'   model equal to the number of isolates in the network. For an undirected
#'   network, an isolate is defined to be any node with degree zero. For a
#'   directed network, an isolate is any node with both in-degree and out-degree
#'   equal to zero.
#'
#' @usage
#' # binary: isolates
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
#' @concept frequently-used
InitErgmTerm.isolates <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                     varnames = NULL,
                     vartypes = NULL,
                     defaultvalues = list(),
                     required = NULL)
  ### Construct the list to return
  list(name="isolates",                               #name: required
       coef.names = "isolates",                       #coef.names: required
       emptynwstats = network.size(nw), # When nw is empty, isolates=n, not 0,
       minval = 0,
       maxval = network.size(nw),
       conflicts.constraints="degreedist"
       )                                                               
}

################################################################################

#' @templateVar name istar
#' @title In-stars
#' @description This term adds one network statistic to the
#'   model for each element in `k` . The \eqn{i} th such statistic counts the
#'   number of distinct `k[i]` -instars in the network, where a
#'   \eqn{k} -instar is defined to be a node \eqn{N} and a set of \eqn{k}
#'   different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the ties
#'   \eqn{(O_j{\rightarrow}N)}{(O_j, N)} exist for \eqn{j=1, \dots, k} . If `attr` is specified
#'   then the count is over the number of \eqn{k} -instars where all nodes have
#'   the same value of the attribute. This term can only be used for directed
#'   networks; for undirected networks see `kstar` . Note that
#'   `istar(1)` is equal to both `ostar(1)` and `edges` .
#'
#' @usage
#' # binary: istar(k, attr=NULL, levels=NULL)
#'
#' @param k a vector of distinct integers
#' @template ergmTerm-attr
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.istar<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("k", "attrname", "levels"),
                        vartypes = c("numeric", "character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL    
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("k", "attr", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels    
  }
  k <- a$k
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
#     Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
  }else{
  }
  lk<-length(k)
  if(lk==0){return(NULL)}
  if(!is.null(attrarg)){
    coef.names <- paste("istar",k,".",attrname,sep="")
    inputs <- c(k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lk
  }else{
    coef.names <- paste("istar",k,sep="")
    inputs <- c(k)
  }
  list(name="istar", coef.names=coef.names, inputs=inputs, minval = 0, conflicts.constraints="idegreedist")
}





#=======================InitErgmTerm functions:  K============================#

################################################################################

#' @templateVar name kstar
#' @title \eqn{k}-stars
#' @description This term adds one
#'   network statistic to the model for each element in `k` . The \eqn{i} th
#'   such statistic counts the number of distinct `k[i]` -stars in the
#'   network, where a \eqn{k} -star is defined to be a node \eqn{N} and a set of
#'   \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the ties
#'   \eqn{\{N, O_i\}}{\{N, O[i]\}} exist for \eqn{i=1, \dots, k} . If this is specified then the count is over
#'   the number of \eqn{k} -stars where all nodes have the same value of the
#'   attribute. This term can only be used for undirected networks; for directed
#'   networks, see `istar` , `ostar` , `twopath` and `m2star` .
#'   Note that `kstar(1)` is equal to `edges` .
#'
#' @usage
#' # binary: kstar(k, attr=NULL, levels=NULL)
#'
#' @param k a vector of distinct integers
#' @template ergmTerm-attr
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.kstar<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("k", "attrname", "levels"),
                        vartypes = c("numeric", "character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("k", "attr", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
  k<-a$k
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
#    Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
  }
  lk<-length(k)
  if(lk==0){return(NULL)}
  if(!is.null(attrarg)){
    coef.names <- paste("kstar",k,".",attrname,sep="")
    inputs <- c(k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lk
  }else{
    coef.names <- paste("kstar",k,sep="")
    inputs <- c(k)
  }
  list(name="kstar", coef.names=coef.names, inputs=inputs, minval = 0, conflicts.constraints="degreedist")
}




#=======================InitErgmTerm functions:  L============================#

################################################################################

#' @templateVar name localtriangle
#' @title Triangles within neighborhoods
#' @description This term adds one statistic to the model equal to the number of triangles
#'   in the network between nodes "close to" each other. For an undirected
#'   network, a local triangle is defined to be any set of three edges between
#'   nodal pairs \eqn{\{(i,j), (j,k), (k,i)\}} that are in the same neighborhood.
#'   For a directed network, a triangle is defined as any set of three edges
#'   \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)} and either
#'   \eqn{(k{\rightarrow}i)} or \eqn{(k{\leftarrow}i)} where again all nodes are
#'   within the same neighborhood.
#'
#' @usage
#' # binary: localtriangle(x)
#'
#' @param x an undirected
#'   network or an symmetric adjacency matrix that specifies whether the two nodes
#'   are in the same neighborhood. Note that `triangle` , with or without an argument, is a
#'   special case of `localtriangle` .
#'
#' @template ergmTerm-general
#'
#' @concept triad-related
#' @concept directed
#' @concept undirected
#' @concept categorical dyadic attribute
InitErgmTerm.localtriangle<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  x<-a$x;attrname<-a$attrname
  if(is.network(x))
    xm<-as.matrix(x, matrix.type="adjacency", attrname)
  else if(is.character(x))
    xm<-as.matrix(nw, matrix.type="adjacency", x)
  else
    xm<-as.matrix(x)
  if(!isSymmetric(xm)){
    ergm_Init_warn("localtriangle requires an undirected neighborhood. Using only mutual ties.")
    xm <- pmin(xm[],(t(xm))[])
  }
  if(!is.null(attrname))
    coef.names <- paste("localtriangle", attrname, sep = ".")
  else
    coef.names <- paste("localtriangle", as.character(sys.call(0)[[3]][2]),
                        sep = ".")
  inputs <- c(NROW(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="localtriangle", coef.names=coef.names, inputs=inputs)
}



#=======================InitErgmTerm functions:  M============================#

################################################################################

#' @templateVar name m2star
#' @title Mixed 2-stars, a.k.a 2-paths
#' @description This term adds one statistic to the model, equal to the number of mixed
#'   2-stars in the network, where a mixed 2-star is a pair of distinct edges
#'   \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)} . A mixed 2-star is
#'   sometimes called a 2-path because it is a directed path of length 2 from
#'   \eqn{i} to \eqn{k} via \eqn{j} . However, in the case of a 2-path the focus
#'   is usually on the endpoints \eqn{i} and \eqn{k} , whereas for a mixed 2-star
#'   the focus is usually on the midpoint \eqn{j} . This term can only be used
#'   with directed networks; for undirected networks see `kstar(2)` . See
#'   also `twopath` .
#'
#' @usage
#' # binary: m2star
#'
#' @template ergmTerm-general
#'
#' @concept directed
InitErgmTerm.m2star<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="m2star", coef.names="m2star",dependence=TRUE, minval = 0) 
}


################################################################################

#' @templateVar name meandeg
#' @title Mean vertex degree
#' @description This term adds one network statistic to the model equal to the
#'   average degree of a node. Note that this term is a constant multiple of
#'   both `edges` and `density` .
#'
#' @usage
#' # binary: meandeg
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmTerm.meandeg<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="meandeg", coef.names="meandeg", dependence=FALSE, minval=0, maxval=if(!is.bipartite(nw)) network.size(nw)-1, conflicts.constraints="edges")
}


################################################################################

#' @templateVar name mm
#' @title Mixing matrix cells and margins
#' @description `attrs` is the rows of the mixing matrix and whose RHS gives
#'   that for its columns. A one-sided formula (e.g.,
#'   `~A` ) is symmetrized (e.g., `A~A` ). A two-sided formula with a dot on one side
#'   calculates the margins of the mixing matrix, analogously to `nodefactor` , with
#'   `A~.` calculating the row/sender/b1 margins and `.~A`
#'   calculating the column/receiver/b2 margins.
#'
#' @usage
#' # binary: mm(attrs, levels=NULL, levels2=-1)
#'
#' @param attrs a two-sided formula whose LHS gives the attribute or
#'   attribute function (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) for the rows of the mixing matrix and whose RHS gives
#'   for its columns. A one-sided formula (e.g., `~A`) is symmetrized (e.g., `A~A`)
#' @templateVar explain subset of rows and columns to be used.
#' @template ergmTerm-levels-doco
#' @param levels2 which specific cells of the matrix to include
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.mm<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.11.0")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrs", "levels", "levels2"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrs", "levels", "levels2"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, -1),
                        required = c(TRUE, FALSE, FALSE))
  }

  # Some preprocessing steps are the same, so run together:
  #' @import purrr
  #' @importFrom utils relist
  spec <-
    list(attrs = a$attrs, levels = a$levels) %>%
    map_if(~!is(., "formula"), ~call("~", .)) %>% # Embed into RHS of formula.
    map_if(~length(.)==2, ~call("~", .[[2]], .[[2]])) %>% # Convert ~X to X~X.
    map(as.list) %>% map(~.[-1]) %>% # Convert to list(X,X).
    map(set_names, c("row", "col")) %>% # Name elements rowspec and colspec.
    transpose() %>%
    unlist(recursive=FALSE) %>% # Convert into a flat list.
    map_if(~is.name(.)&&.==".", ~NULL) %>% # If it's just a dot, convert to NULL.
    map_if(~is.call(.)||(is.name(.)&&.!="."), ~as.formula(call("~", .))) %>% # If it's a call or a symbol, embed in formula.
    relist(skeleton=list(row=c(attrs=NA, levels=NA), col=c(attrs=NA, levels=NA))) %>% # Reconstruct list.
    transpose()

  if(is(a$attrs, "formula"))
    spec[["attrs"]] <- lapply(spec[["attrs"]], function(x){if(is(x,"formula")) environment(x) <- environment(a$attrs); x})
  if(is(a$levels, "formula"))
    spec[["levels"]] <- lapply(spec[["levels"]], function(x){if(is(x,"formula")) environment(x) <- environment(a$levels); x})
  spec <- transpose(spec)
  
  # Extract attribute values.
  attrval <-
    spec %>%
    imap(function(spec, whose){
      if(is.null(spec$attrs)){
        list(valcodes =
               rep(0L,
                   if(!is.bipartite(nw)) network.size(nw)
                   else if(whose=="row") nw%n%"bipartite"
                   else if(whose=="col") network.size(nw) - nw%n%"bipartite"
                   ),
             name = ".",
             levels = NA,
             levelcodes = 0
             )
      }else{
        x <- ergm_get_vattr(spec$attrs, nw, bip = if(is.bipartite(nw)) c(row="b1",col="b2")[whose] else "n")
        name <- attr(x, "name")
        list(name=name, val=x, levels=spec$levels, unique=sort(unique(x)))
      }
    })

  # Undirected unipartite networks with identical attribute
  # specification produce square, symmetric mixing matrices. All
  # others do not.
  symm <- !is.directed(nw) && !is.bipartite(nw) && identical(spec$row$attrs, spec$col$attrs)
  # Are we evaluating the margin?
  marg <- length(attrval$row$unique)==0 || length(attrval$col$unique)==0
  
  # Filter the final level set and encode the attribute values.
  attrval <- attrval %>%
    map_if(~is.null(.$levelcodes), function(v){
      v$levels <- ergm_attr_levels(v$levels, v$val, nw, levels=v$unique)
      v$levelcodes <- seq_along(v$levels)
      v$valcodes <- match(v$val, v$levels, nomatch=0)
      v
    })

  # Construct all pairwise level combinations (table cells) and their numeric codes.
  levels2codes <- expand.grid(row=attrval$row$levelcodes, col=attrval$col$levelcodes) %>% transpose()
  levels2 <- expand.grid(row=attrval$row$levels, col=attrval$col$levels, stringsAsFactors=FALSE) %>% transpose()

  # Drop redundant table cells if symmetrising.
  if(symm){
    levels2keep <- levels2codes %>% map_lgl(with, row <= col)
    levels2codes <- levels2codes[levels2keep]
    levels2 <- levels2[levels2keep]
  }

  # Run the table cell list through the cell filter.
  levels2sel <- ergm_attr_levels(a$levels2, list(row=attrval$row$val, col=attrval$col$val), nw, levels=levels2)
  if(length(levels2sel) == 0) return(NULL)
  levels2codes <- levels2codes[match(levels2sel,levels2, NA)]
  levels2 <- levels2sel; rm(levels2sel)

  # Construct the level names
  levels2names <-
    levels2 %>%
    transpose() %>%
    map(unlist) %>%
    with(paste0(
      "[",
      if(attrval$row$name!=".")
        paste0(attrval$row$name, "=", .$row)
      else ".",
      ",",
      if(attrval$col$name!=".")
        paste0(attrval$col$name, "=", .$col)
      else ".",
      "]"))
  
  coef.names <- paste0("mm",levels2names)

  list(name = "mixmat",
       coef.names = coef.names,
       inputs = c(symm+marg*2, attrval$row$valcodes, attrval$col$valcodes, unlist(levels2codes)),
       dependence = FALSE,
       minval = 0)
}


################################################################################

#' @templateVar name mutual
#' @title Mutuality
#' @description In binary ERGMs, equal to the number of
#'   pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
#'   and \eqn{(j{\rightarrow}i)}{(j,i)} both exist. For valued ERGMs, equal to \eqn{\sum_{i<j} m(y_{i,j},y_{j,i})} ,
#'   where \eqn{m} is determined by `form` argument: `"min"`
#'   for \eqn{\min(y_{i,j},y_{j,i})} , `"nabsdiff"` for
#'   \eqn{-|y_{i,j},y_{j,i}|} , `"product"` for
#'   \eqn{y_{i,j}y_{j,i}} , and `"geometric"` for
#'   \eqn{\sqrt{y_{i,j}}\sqrt{y_{j,i}}} . See Krivitsky (2012) for a
#'   discussion of these statistics. `form="threshold"` simply
#'   computes the binary `mutuality` after
#'   thresholding at `threshold` .
#'   
#'   This term can only be used with directed networks. 
#'   
#' @usage
#' # binary: mutual(same=NULL, by=NULL, diff=FALSE, keep=NULL, levels=NULL)
#'
#' @param same if the optional argument is passed
#'   (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details),
#'   only mutual pairs that match on the attribute are counted;
#'   separate counts for each unique matching value can be obtained by using
#'   `diff=TRUE` with `same`. Only one of `same` or `by` may be used. If both parameters are used, `by` is 
#'   ignored. This paramer is affected by `diff`.
#'
#' @param by if the optional argument is passed (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details),
#'   then each node is counted separately for each mutual pair in which it
#'   occurs and the counts are tabulated by unique values of the attribute.
#'   This means that the sum of the mutual statistics when `by` is used
#'   will equal twice the standard mutual statistic. Only one of `same` or `by` may be used. If both parameters are used, `by` is 
#'   ignored. This paramer is not affected by `diff`.
#' @param keep deprecated
#' @templateVar explain which statistics should be kept whenever the `mutual` term would ordinarily result in multiple statistics.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-keep-dep
#'
#' @concept directed
#' @concept frequently-used
InitErgmTerm.mutual<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                        varnames = c("same", "by", "diff", "keep"),
                        vartypes = c("character", "character", "logical", "numeric"),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, FALSE, "levels"))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                        varnames = c("same", "by", "diff", "keep", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC, "logical", "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, FALSE, NULL, NULL),
                        required = c(FALSE, FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, FALSE, "levels", FALSE))
  }
  
  
  ### Process the arguments
  if (!is.null(a$same) || !is.null(a$by)) {
    if (!is.null(a$same)) {
     attrarg <- a$same
     if (!is.null(a$by)) 
       ergm_Init_warn("Ignoring 'by' argument to mutual because 'same' exists")
    }else{
     attrarg <- a$by
    }
    
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(a$levels, nodecov, nw, levels = sort(unique(nodecov)))
    if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && !is.null(a$keep)) u <- u[a$keep]
    
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }

  ### Construct the list to return
  if (!is.null(a$same) || !is.null(a$by)) {
    if (is.null(a$same)) {
      coef.names <- paste("mutual.by", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
    }else{
     if (a$diff) {
      coef.names <- paste("mutual.same", attrname, u, sep=".")
      inputs <- c(ui, nodecov)
     }else{ 
      coef.names <- paste("mutual", attrname, sep=".")
      inputs <- nodecov
     }
    }
    if (is.null(a$same) && !is.null(a$by)) {
     name <- "mutual_by_attr"
    }else{
     name <- "mutual"
    }
  }else{
     name <- "mutual"
     coef.names <- "mutual"
     inputs <- NULL
  }

  maxval <- network.dyadcount(nw,FALSE)/2

  list(name=name,                      #name: required
       coef.names = coef.names,        #coef.names: required
       inputs=inputs,
       minval = 0,
       maxval = maxval) 
}


#=======================InitErgmTerm functions:  N============================#


################################################################################

#' @templateVar name nearsimmelian
#' @title Near simmelian triads
#' @description This term adds one statistic to the model equal to the number of near
#'   Simmelian triads, as defined by Krackhardt and Handcock (2007). This is a
#'   sub-graph of size three which is exactly one tie short of being complete.
#'
#' @usage
#' # binary: nearsimmelian
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept directed
#' @concept triad-related
InitErgmTerm.nearsimmelian<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nearsimmelian", coef.names="nearsimmelian", minval=0, maxval=network.dyadcount(nw,FALSE)*network.size(nw)*0.5)
}


################################################################################

#' @templateVar name nodecov
#' @title Main effect of a covariate
#' @description This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the sum of
#'   `attr(i)` and `attr(j)` for all edges \eqn{(i,j)} in the
#'   network. For categorical attributes, see `nodefactor` . Note that for
#'   directed networks, `nodecov` equals `nodeicov` plus
#'   `nodeocov` .
#'
#' @usage
#' # binary: nodecov(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept dyad-independent
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept quantitative nodal attribute
InitErgmTerm.nodecov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname","transform","transformname"),
                        vartypes = c("character","function","character"),
                        defaultvalues = list(NULL,function(x)x,""),
                        required = c(TRUE,FALSE,FALSE))
    attrname<-a$attrname
    f<-a$transform
    f.name<-a$transformname
    coef.names <- paste(paste("nodecov",f.name,sep=""),attrname,sep=".")
    nodecov <- f(get.node.attr(nw, attrname, "nodecov", numeric=TRUE))
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                        varnames = c("attr"),
                        vartypes = c(ERGM_VATTR_SPEC),
                        defaultvalues = list(NULL),
                        required = c(TRUE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", multiple="matrix")
    coef.names <- nodecov_names(nodecov, "nodecov")
  }
  list(name="nodecov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}

#' @templateVar name nodecov
#' @template ergmTerm-rdname
#' @aliases nodemain-ergmTerm
#' @usage
#' # binary: nodemain
InitErgmTerm.nodemain<-InitErgmTerm.nodecov

################################################################################

#' @templateVar name nodefactor
#' @title Factor attribute effect
#' @description This term adds multiple network statistics to the
#'   model, one for each of (a subset of) the unique values of the
#'   `attr` attribute (or each combination of the attributes
#'   given). Each of these statistics gives the number of times a node
#'   with that attribute or those attributes appears in an edge in the
#'   network.
#'   
#' @usage
#' # binary: nodefactor(attr, base=1, levels=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-base-dep
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.nodefactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "base", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr                        
    levels <- a$levels    
  }

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

  if (attr(a,"missing")["levels"] && any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
  }

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  } 
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)-1
  ### Construct the list to return
  inputs <- nodepos
  list(name="nodefactor",                                        #required
       coef.names = paste("nodefactor", paste(attrname,collapse="."), u, sep="."), #required
       iinputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################

#' @templateVar name nodeicov
#' @title Main effect of a covariate for in-edges
#' @description This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
#'   value of `attr(j)` for all edges \eqn{(i,j)} in the network. This
#'   term may only be used with directed networks. For categorical attributes,
#'   see `nodeifactor` .
#'
#' @usage
#' # binary: nodeicov(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept directed
#' @concept quantitative nodal attribute
#' @concept frequently-used
InitErgmTerm.nodeicov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attrname","transform","transformname"),
                        vartypes = c("character","function","character"),
                        defaultvalues = list(NULL,function(x)x,""),
                        required = c(TRUE,FALSE,FALSE))
    ### Process the arguments
    attrname<-a$attrname
    f<-a$transform
    f.name<-a$transformname
    coef.names <- paste(paste("nodeicov",f.name,sep=""),attrname,sep=".")
    nodecov <- f(get.node.attr(nw, attrname, "nodeicov", numeric=TRUE))
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attr"),
                        vartypes = c(ERGM_VATTR_SPEC),
                        defaultvalues = list(NULL),
                        required = c(TRUE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", multiple="matrix")
    coef.names <- nodecov_names(nodecov, "nodeicov")
  }
  list(name="nodeicov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}



################################################################################

#' @templateVar name nodeifactor
#' @title Factor attribute effect for in-edges
#' @description This term adds multiple network
#'   statistics to the model, one for each of (a subset of) the unique
#'   values of the `attr` attribute (or each combination of the
#'   attributes given). Each of these statistics gives the number of
#'   times a node with that attribute or those attributes appears as the
#'   terminal node of a directed tie.
#'   
#'   For an analogous term for quantitative vertex attributes, see `nodeicov` .
#'
#' @usage
#' # binary: nodeifactor(attr, base=1, levels=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @concept dyad-independent
#' @concept directed
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.nodeifactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr", "base", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr                        
    levels <- a$levels    
  }

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

  if (attr(a,"missing")["levels"] && any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
  }

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  } 
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)-1
  ### Construct the list to return
  inputs <- nodepos
  list(name="nodeifactor",                                        #required
       coef.names = paste("nodeifactor", paste(attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################

#' @templateVar name nodematch
#' @title Uniform homophily and differential homophily
#' @description When `diff=FALSE` , this term adds one network statistic
#'   to the model, which counts the number of edges \eqn{(i,j)} for which
#'   `attr(i)==attr(j)` . This is also called \dQuote{uniform homophily}, because each group is assumed to have the same propensity for within-group ties. When multiple attribute names are given, the
#'   statistic counts only ties for which all of the attributes
#'   match. When `diff=TRUE` , \eqn{p} network statistics are added
#'   to the model, where \eqn{p} is the number of unique values of the
#'   `attr` attribute. The \eqn{k} th such statistic counts the
#'   number of edges \eqn{(i,j)} for which `attr(i) == attr(j) == value(k)` , where `value(k)` is the \eqn{k} th
#'   smallest unique value of the `attr` attribute. This is also called \dQuote{differential homophily}, because each group is allowed to have a unique propensity for within-group ties. Note that a statistical test of uniform vs. differential homophily should be conducted using the ANOVA function.
#'   
#'   By default, matches on all levels \eqn{k} are
#'   counted. This works for both
#'   `diff=TRUE` and `diff=FALSE` .
#'
#' @usage
#' # binary: nodematch(attr, diff=FALSE, keep=NULL, levels=NULL)
#'
#' @template ergmTerm-attr
#' @param diff specify if the term has uniform or differential homophily
#' @param keep deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-keep-dep
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.nodematch<-InitErgmTerm.match<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, 
                        varnames = c("attrname", "diff", "keep", "levels"),
                        vartypes = c("character", "logical", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, 
                        varnames = c("attr", "diff", "keep", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, FALSE, "levels", FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
                        
  ### Process the arguments
  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
  if(attr(a,"missing")["levels"] && !is.null(a$keep)) u <- u[a$keep]
  
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov==(length(u)+1)
  nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
  ui <- seq(along=u)
  ### Construct the list to return
  if (a$diff) {
    coef.names <- paste("nodematch", paste(attrname,collapse="."), u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", paste(attrname,collapse="."), sep=".")
    inputs <- nodecov
  }
  list(name="nodematch",                                 #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################

#' @templateVar name nodemix
#' @title Nodal attribute mixing
#' @description By default, this term adds one network statistic to
#'   the model for each possible pairing of attribute values. The
#'   statistic equals the number of edges in the network in which the
#'   nodes have that pairing of values. (When multiple attributes are specified, a
#'   statistic is added for each combination of attribute values for
#'   those attributes.) In other words, this term produces one statistic for
#'   every entry in the mixing matrix for the attribute(s). By default, the ordering of
#'   the attribute values is lexicographic: alphabetical (for nominal categories) or
#'   numerical (for ordered categories).
#'   
#' @usage
#' # binary: nodemix(attr, base=NULL, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @param b1levels,b2levels,levels control what statistics are included in the model and the order in which they appear. `levels` applies to unipartite networks; `b1levels` and `b2levels` apply to bipartite networks (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#' @param levels2 similar to the other levels arguments above and applies to all networks. Optionally allows a factor or character matrix to be specified to group certain levels. Level combinations corresponding to `NA` are excluded. Combinations specified by the same character or level will be grouped together and summarised by the same statistic. If an empty string is specified, the level combinations will be ungrouped. Only the upper triangle needs to be specified for undirected networks. For example, `levels2=matrix(c('A', '', NA, 'A'), 2, 2, byrow=TRUE)` on an undirected matrix will group homophilous ties while leaving ties between 1 and 2 ungrouped.
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-base-dep2
#'
#' @concept dyad-independent
#' @concept frequently-used
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.nodemix<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "base", "b1levels", "b2levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels2", FALSE, FALSE))
    attrarg <- a$attrname
    b1levels <- if(!is.null(a$b1levels)) I(a$b1levels) else NULL
    b2levels <- if(!is.null(a$b2levels)) I(a$b2levels) else NULL
  }else if(version <= as.package_version("3.11.0")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "base", "b1levels", "b2levels", "levels", "levels2"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL, NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels2", FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    b1levels <- a$b1levels
    b2levels <- a$b2levels
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "base", "b1levels", "b2levels", "levels", "levels2"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL, NULL, NULL, -1),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels2", FALSE, FALSE, FALSE, FALSE))
    attrarg <- a$attr
    b1levels <- a$b1levels
    b2levels <- a$b2levels
  }

  ### Process the arguments
  if (is.bipartite(nw) && is.directed(nw)) {
    ergm_Init_abort("Directed bipartite networks are not currently possible")
  }

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  
  if (is.bipartite(nw)) {
    #  So undirected network storage but directed mixing
    
    b1nodecov <- ergm_get_vattr(attrarg, nw, bip = "b1")
    b2nodecov <- ergm_get_vattr(attrarg, nw, bip = "b2")
    
    b1namescov <- ergm_attr_levels(b1levels, b1nodecov, nw, sort(unique(b1nodecov)))
    b2namescov <- ergm_attr_levels(b2levels, b2nodecov, nw, sort(unique(b2nodecov)))
    
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    
    levels2.list <- transpose(expand.grid(row = b1namescov, col = b2namescov, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = nr + 1:nc)

    if ((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) {
      levels2.sel <- levels2.list[-a$base]
      has.groups <- FALSE
    } else if (!is.character(a$levels2)) {
      levels2.sel <- ergm_attr_levels(a$levels2, list(row = b1nodecov, col = b2nodecov), nw, levels2.list)
      has.groups <- FALSE
    } else {
      levels2.sel <- ergm_attr_levels(!is.na(a$levels2) & (a$levels2 == ''), list(row = b1nodecov, col = b2nodecov), nw, levels2.list)
      has.groups <- TRUE
    }

    indmat <- matrix(0L, nrow = nr, ncol = nc)
    cn <- c()

    if (has.groups) {
      for (g in sort(unique(as.vector(a$levels2[!is.na(a$levels2) & a$levels2 != ''])))) {
        if (g != '') {
          cn <- c(cn, paste("mix", paste(attrname, collapse="."), g, sep="."))
          indmat[a$levels2 == g] <- length(cn)
          if (!is.directed(nw)) {
            indmat[t(a$levels2 == g)] <- length(cn)
          }
        }
      }
    }

    if (length(levels2.sel) > 0) {
      rows2keep <- match(levels2.sel,levels2.list, NA)
      rows2keep <- rows2keep[!is.na(rows2keep)]

      u <- indices2.grid[rows2keep,]

      ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
      u[,2L] <- u[,2L] - nr
      indmat[as.matrix(u)] <- seq_len(NROW(u)) + length(cn)

      namescov <- c(b1namescov, b2namescov)
      cn <- c(cn, paste("mix", paste(attrname,collapse="."), apply(matrix(namescov[as.matrix(u)],ncol=2),
                                         1,paste,collapse="."), sep="."))
    }

    # Recode to numeric
    b1nodecov <- match(b1nodecov,b1namescov,nomatch=length(b1namescov)+1)
    b2nodecov <- match(b2nodecov,b2namescov,nomatch=length(b2namescov)+1)
    nodecov <- c(b1nodecov, b2nodecov)

    indmat <- cbind(rbind(indmat, 0L), 0L) - 1L

  } else { # So one mode, but could be directed or undirected
    u <- ergm_attr_levels(a$levels, nodecov, nw, sort(unique(nodecov)))
    namescov <- u
    
    nr <- length(u)
    nc <- length(u)

    levels2.list <- transpose(expand.grid(row = u, col = u, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
    uun <- as.vector(outer(u,u,paste,sep="."))
    
    if (!is.directed(nw)) {
        rowleqcol <- indices2.grid$row <= indices2.grid$col
        levels2.list <- levels2.list[rowleqcol]
        indices2.grid <- indices2.grid[rowleqcol,]
        uun <- uun[rowleqcol]
    }

    if ((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) {
      levels2.sel <- levels2.list[-a$base]
      has.groups <- FALSE
    } else if (!is.character(a$levels2)) {
      levels2.sel <- ergm_attr_levels(a$levels2, list(row = nodecov, col = nodecov), nw, levels2.list)
      has.groups <- FALSE
    } else {
      levels2.sel <- ergm_attr_levels(!is.na(a$levels2) & (a$levels2 == ''), list(row = nodecov, col = nodecov), nw, levels2.list)
      has.groups <- TRUE
    }

    indmat <- matrix(0L, nrow=nr, ncol=nc)
    cn <- c()

    if (has.groups) {
      for (g in sort(unique(as.vector(a$levels2[!is.na(a$levels2) & a$levels2 != ''])))) {
        if (g != '') {
          cn <- c(cn, paste("mix", paste(attrname, collapse="."), g, sep="."))
          indmat[a$levels2 == g] <- length(cn)
          if (!is.directed(nw)) {
            indmat[t(a$levels2 == g)] <- length(cn)
          }
        }
      }
    }

    if (length(levels2.sel) > 0) {
      indmat.ungrouped <- matrix(0L, nrow=nr, ncol=nc)
      rows2keep <- match(levels2.sel,levels2.list, NA)
      rows2keep <- rows2keep[!is.na(rows2keep)]

      u <- indices2.grid[rows2keep,]
      uun <- uun[rows2keep]

      ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
      indmat.ungrouped[as.matrix(u)] <- seq_len(NROW(u)) + length(cn)
      if(!is.directed(nw)) indmat.ungrouped <- indmat.ungrouped + t(indmat.ungrouped) - diag(diag(indmat.ungrouped))

      cn <- c(cn, paste("mix", paste(attrname,collapse="."), uun, sep="."))
      indmat <- indmat + indmat.ungrouped
    }

    indmat <- cbind(rbind(indmat, 0), 0) - 1L
    nodecov <- match(nodecov,namescov,nomatch=length(namescov)+1)
  }

  ### Construct the list to return
  list(name = "nodemix", coef.names = cn, # required
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0,
       inputs = NULL, # passed by name below
       nr = as.integer(nr + 1),
       nc = as.integer(nc + 1),
       indmat = as.integer(t(indmat)),
       nodecov = as.integer(c(0L, nodecov) - 1L) # two shifts to make the C code cleaner
       )
}

################################################################################

#' @templateVar name nodeocov
#' @title Main effect of a covariate for out-edges
#' @description This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
#'   value of `attr(i)` for all edges \eqn{(i,j)} in the network. This
#'   term may only be used with directed networks. For categorical attributes,
#'   see `nodeofactor` .
#'   
#' @usage
#' # binary: nodeocov(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-args-3.9.4
#'
#' @concept directed
#' @concept dyad-independent
#' @concept quantitative nodal attribute
InitErgmTerm.nodeocov<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname","transform","transformname"),
                        vartypes = c("character","function","character"),
                        defaultvalues = list(NULL,function(x)x,""),
                        required = c(TRUE,FALSE,FALSE))
    ### Process the arguments
    attrname<-a$attrname
    f<-a$transform
    f.name<-a$transformname
    coef.names <- paste(paste("nodeocov",f.name,sep=""),attrname,sep=".")
    nodecov <- f(get.node.attr(nw, attrname, "nodeocov", numeric=TRUE))
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr"),
                        vartypes = c(ERGM_VATTR_SPEC),
                        defaultvalues = list(NULL),
                        required = c(TRUE))
    ### Process the arguments
    nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", multiple="matrix")
    coef.names <- nodecov_names(nodecov, "nodeocov")
  }
  list(name="nodeocov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}



################################################################################

#' @templateVar name nodeofactor
#' @title Factor attribute effect for out-edges
#' @description This term adds multiple network
#'   statistics to the model, one for each of (a subset of) the unique
#'   values of the `attr` attribute (or each combination of the
#'   attributes given). Each of these statistics gives the number of
#'   times a node with that attribute or those attributes appears as the
#'   node of origin of a directed tie.
#'   
#' @usage
#' # binary: nodeofactor(attr, base=1, levels=-1)
#'
#' @template ergmTerm-attr
#' @param base deprecated
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-levels-not-first
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept dyad-independent
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.nodeofactor<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                        varnames = c("attr", "base", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, LEVELS_BASE1),
                        required = c(TRUE, FALSE, FALSE),
                        dep.inform = list(FALSE, "levels", FALSE))
    attrarg <- a$attr                        
    levels <- a$levels    
  }

  nodecov <- ergm_get_vattr(attrarg, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

  if (attr(a,"missing")["levels"] && any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
  }

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  }  
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)-1

  ### Construct the list to return
  inputs <- nodepos
  list(name="nodeofactor",                                        #required
       coef.names = paste("nodeofactor", paste(attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################

#' @templateVar name nsp
#' @title Nonedgewise shared partners
#' @description This is
#'   just like the `dsp` and `esp` terms, except this term adds
#'   one network statistic to the model for each element in `d`
#'   where the \eqn{i} th such statistic equals the number of
#'   non-edges (that is, dyads that do not have an edge) in the network
#'   with exactly `d[i]` shared partners. This term can be used with
#'   directed and undirected networks.
#'   
#' @usage
#' # binary: nsp(d)
#'
#' @param d a vector of distinct integers
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @templateVar fn nsp
#' @templateVar kind (directed) non-edge `(i,j)`
#' @templateVar see dnsp
#' 
#' @template ergmTerm-sp-to-dsp
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.nsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d<-a$d
  if (any(d==0)) {
    emptynwstats <- rep(0, length(d))
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw, "bipartite")
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw,FALSE)
    }
  }else{
    emptynwstats <- NULL
  }
  ld<-length(d)
  if(ld==0){return(NULL)}
  coef.names <- paste("nsp",d,sep="")
  if(is.directed(nw)){dname <- "tnsp"}else{dname <- "nsp"}

  if (!is.null(emptynwstats)) {
    list(name=dname, coef.names=coef.names, inputs=c(d),
         emptynwstats=emptynwstats, minval=0, auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  } else {
    list(name=dname, coef.names=coef.names, inputs=c(d), minval=0, auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}



#=======================InitErgmTerm functions:  O============================#

################################################################################

#' @templateVar name odegrange
#' @title Out-degree range
#' @description This term adds one
#'   network statistic to the model for each element of `from` (or `to` ); the \eqn{i} th
#'   such statistic equals the number of nodes in the network of out-degree
#'   greater than or equal to
#'   `from[i]` but strictly less than `to[i]` , i.e. with
#'   out-edge count
#'   in semiopen interval `[from,to)` . 
#'
#'   This term can only be used with directed networks; for undirected
#'   networks (bipartite and not)
#'   see `degrange` . For degrees of specific modes of bipartite
#'   networks, see `b1degrange` and `b2degrange` . For
#'   in-degrees, see `idegrange` .
#'
#' @usage
#' # binary: odegrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @template ergmTerm-from-to
#' @template ergmTerm-by
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-attr
#'
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.odegrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degrange_impl("o", TRUE, NULL, nw, arglist, ..., version=version)
}

################################################################################

#' @templateVar name odegree
#' @title Out-degree
#' @description This term adds one network statistic to
#'   the model for each element in `d` ; the \eqn{i} th such statistic equals
#'   the number of nodes in the network of out-degree `d[i]` , i.e. the
#'   number of nodes with exactly `d[i]` out-edges. 
#'   This term can only be used with directed networks; for undirected networks
#'   see `degree` .
#'
#' @usage
#' # binary: odegree(d, by=NULL, homophily=FALSE, levels=NULL)
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-by
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept categorical nodal attribute
#' @concept frequently-used
InitErgmTerm.odegree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  .degree_impl("o", TRUE, NULL, nw, arglist, ..., version=version)
}



################################################################################

#' @templateVar name odegree1.5
#' @title Out-degree to the 3/2 power
#' @description This term adds one network statistic to the model equaling the sum
#'   over the actors of each actor's outdegree taken to the 3/2 power
#'   (or, equivalently, multiplied by its square root). This term is
#'   analogous to the term of Snijders et al. (2010), equation (12). This
#'   term can only be used with directed networks.
#'
#' @usage
#' # binary: odegree1.5
#'
#' @template ergmTerm-general
#'
#' @concept directed
InitErgmTerm.odegree1.5<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="odegreepopularity", coef.names="odegree1.5",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="odegreedist")
}


################################################################################
#' @describeIn ergm-deprecated Use [`odegree1.5`][odegree1.5-ergmTerm] instead.
InitErgmTerm.odegreepopularity<-function (nw, arglist, ...) {
  .Deprecated("odegree1.5")
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="odegreepopularity", coef.names="odegreepopularity",
       minval=0, maxval=network.dyadcount(nw,FALSE)*sqrt(network.size(nw)-1), conflicts.constraints="odegreedist")
}


################################################################################

#' @templateVar name opentriad
#' @title Open triads
#' @description This term
#'   adds one statistic to the model equal to the number of 2-stars minus
#'   three times the number of triangles in the network. It is currently
#'   only implemented for undirected networks.
#'
#' @usage
#' # binary: opentriad
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept triad-related
InitErgmTerm.opentriad<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  list(name="opentriad", coef.names="opentriad", inputs=NULL)
}


################################################################################

#' @templateVar name ostar
#' @title k-Outstars
#' @description This term adds one network statistic to the
#'   model for each element in `k` . The \eqn{i} th such statistic counts the
#'   number of distinct `k[i]` -outstars in the network, where a
#'   \eqn{k} -outstar is defined to be a node \eqn{N} and a set of \eqn{k}
#'   different nodes \eqn{\{O_1, \dots, O_k\}}{{O[1], ..., O[k]}} such that the ties
#'   \eqn{(N{\rightarrow}O_j)}{(N,O_j)} exist for \eqn{j=1, \dots, k} . If `attr` is specified
#'   then the count is the number of \eqn{k} -outstars where all nodes have the
#'   same value of the attribute. This term can only be used with directed
#'   networks; for undirected networks see `kstar` .
#'
#' @usage
#' # binary: ostar(k, attr=NULL, levels=NULL)
#'
#' @param k a vector of distinct integers
#' @template ergmTerm-attr
#' @templateVar explain TODO
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @note `ostar(1)` is equal to both `istar(1)` and `edges` .
#'
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.ostar<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("k", "attrname", "levels"),
                        vartypes = c("numeric", "character", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL    
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("k", "attr", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, NULL),
                        required = c(TRUE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels    
  }
  k<-a$k
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
  }
  lk<-length(k)
  if(lk==0){return(NULL)}

  if(!is.null(attrarg)){
    coef.names <- paste("ostar",k,".",attrname,sep="")
    inputs <- c(k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lk
  }else{
    coef.names <- paste("ostar",k,sep="")
    inputs <- c(k)
  }
  list(name="ostar", coef.names=coef.names, inputs=inputs, minval=0, conflicts.constraints="odegreedist")  
}



#=======================InitErgmTerm functions:  P============================#




#=======================InitErgmTerm functions:  R============================#



################################################################################

#' @templateVar name receiver
#' @title Receiver effect
#' @description This term adds one network statistic for each node equal to the number of
#'   in-ties for that node. This measures the popularity of the node. The term
#'   for the first node is omitted by default because of linear dependence that
#'   arises if this term is used together with `edges` , but its coefficient
#'   can be computed as the negative of the sum of the coefficients of all the
#'   other actors. That is, the average coefficient is zero, following the
#'   Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
#'   1981).  This
#'   term can only be used with directed networks. For undirected networks, see
#'   `sociality` .
#'
#' @usage
#' # binary: receiver(base=1, nodes=-1)
#'
#' @param base deprecated
#' @param nodes specify which nodes' statistics should be included or excluded (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#'
#' @template ergmTerm-base-dep-node
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept dyad-independent
InitErgmTerm.receiver<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base"),
                        vartypes = c("numeric"),
                        defaultvalues = list(1),
                        required = c(FALSE),
                        dep.inform = list("nodes"))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base", "nodes"),
                        vartypes = c("numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(1, LEVELS_BASE1),
                        required = c(FALSE, FALSE),
                        dep.inform = list("nodes", FALSE))
  }
  d <- ergm_attr_levels(a$nodes, 1:network.size(nw), nw, 1:network.size(nw))
  if((!hasName(attr(a,"missing"), "nodes") || attr(a,"missing")["nodes"]) && any(NVL(a$base,0)!=0)) d <- d[-a$base]
  
  ld<-length(d)
  if(ld==0){return(NULL)}
  list(name="receiver", coef.names=paste("receiver",d,sep=""),
       inputs=c(d), emptynwstats=rep(0,length(d)), dependence=FALSE, minval=0, maxval=network.size(nw)-1, conflicts.constraints="idegrees")
}


#=======================InitErgmTerm functions:  S============================#

################################################################################

#' @templateVar name sender
#' @title Sender effect
#' @description This term adds one network statistic for each node equal to the number of
#'   out-ties for that node. This measures the activity of the node. The term for
#'   the first node is omitted by default because of linear dependence that
#'   arises if this term is used together with `edges` , but its coefficient
#'   can be computed as the negative of the sum of the coefficients of all the
#'   other actors. That is, the average coefficient is zero, following the
#'   Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
#'   1981). 
#'   
#'   For undirected networks, see `sociality` .
#'
#' @usage
#' # binary: sender(base=1, nodes=-1)
#'
#' @param base deprecated
#' @param nodes specify which nodes' statistics should be included or excluded (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details)
#'
#' @template ergmTerm-base-dep-node
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept directed
#' @concept dyad-independent
InitErgmTerm.sender<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base"),
                        vartypes = c("numeric"),
                        defaultvalues = list(1),
                        required = c(FALSE),
                        dep.inform = list("nodes"))
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("base", "nodes"),
                        vartypes = c("numeric", ERGM_LEVELS_SPEC),
                        defaultvalues = list(1, LEVELS_BASE1),
                        required = c(FALSE, FALSE),
                        dep.inform = list("nodes", FALSE))
  }
  d <- ergm_attr_levels(a$nodes, 1:network.size(nw), nw, 1:network.size(nw))
  if((!hasName(attr(a,"missing"), "nodes") || attr(a,"missing")["nodes"]) && any(NVL(a$base,0)!=0)) d <- d[-a$base]
  
  ld<-length(d)
  if(ld==0){return(NULL)}
  list(name="sender", coef.names=paste("sender",d,sep=""),
       inputs=c(d), emptynwstats=rep(0,length(d)), dependence=FALSE, minval=0, maxval=network.size(nw)-1, conflicts.constraints="odegrees")
}

################################################################################

#' @templateVar name simmelian
#' @title Simmelian triads
#' @description This term adds one
#'   statistic to the model equal to the number of Simmelian triads, as defined
#'   by Krackhardt and Handcock (2007). This is a complete sub-graph of size
#'   three.
#'
#' @usage
#' # binary: simmelian
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept directed
#' @concept triad-related
InitErgmTerm.simmelian<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="simmelian", coef.names="simmelian", minval=0, maxval=network.edgecount(nw)*network.size(nw)*0.5)
}


################################################################################

#' @templateVar name simmelianties
#' @title Ties in simmelian triads
#' @description This term adds
#'   one statistic to the model equal to the number of ties in the network that
#'   are associated with Simmelian triads, as defined by Krackhardt and Handcock
#'   (2007). Each Simmelian has six ties in it but, because Simmelians can
#'   overlap in terms of nodes (and associated ties), the total number of ties in
#'   these Simmelians is less than six times the number of Simmelians. Hence this
#'   is a measure of the clustering of Simmelians (given the number of
#'   Simmelians).
#'
#' @usage
#' # binary: simmelianties
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept triad-related
#' @concept directed
InitErgmTerm.simmelianties<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="simmelianties", coef.names="simmelianties", minval=0, maxval=network.edgecount(nw)) # TODO: Is this correct?
}



################################################################################

#' @templateVar name smalldiff
#' @title Number of ties between actors with similar attribute values
#' @description This term adds one statistic, having as its
#'   value the number of edges in the network for which the incident
#'   actors' attribute values differ less than `cutoff` ; that is,
#'   number of edges between `i` to `j` such that
#'   `abs(attr[i]-attr[j])<cutoff` .
#'
#' @usage
#' # binary: smalldiff(attr, cutoff)
#'
#' @template ergmTerm-attr
#' @param maximum difference in attribute values for ties to be considered
#'
#' @template ergmTerm-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @concept quantitative nodal attribute
InitErgmTerm.smalldiff<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "cutoff"),
                        vartypes = c("character", "numeric"),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, TRUE))
    attrarg <- a$attrname
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "cutoff"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric"),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, TRUE))
    attrarg <- a$attr
  }
  
  cutoff <- a$cutoff
  if (length(cutoff)>1)
    ergm_Init_abort("cutoff for smalldiff() must be a scalar.")

  nodecov <- ergm_get_vattr(attrarg, nw, accept="numeric")
  attrname <- attr(nodecov, "name")
  
  coef.names <- paste("smalldiff.", attrname, cutoff, sep="")
  inputs <- c(cutoff, nodecov)
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="smalldiff", coef.names=coef.names, inputs=inputs,
       dependence=FALSE)
}


  

################################################################################

#' @templateVar name sociality
#' @title Undirected degree
#' @description This term adds one network statistic for each node equal to the number of
#'   ties of that node. For directed networks, see `sender` and
#'   `receiver` . 
#'
#' @usage
#' # binary: sociality(attr=NULL, base=1, levels=NULL, nodes=-1)
#'
#' @param attr,levels this optional argument is deprecated and will be replaced with a more elegant implementation in a future release. In the meantime, it specifies a categorical vertex attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details). If provided, this term only counts ties between nodes with the same value of the attribute (an actor-specific version of the `nodematch` term), restricted to be one of the values specified by (also deprecated) `levels` if `levels` is not `NULL` .
#' @param base deprecated
#' @param nodes By default, `nodes=-1` means that the statistic for the
#'   first node will be omitted, but this argument may be changed to control
#'   which statistics are included just as for the `nodes` argument of `sender` and
#'   `receiver` terms.
#'
#' @template ergmTerm-base-dep
#'
#' @template ergmTerm-base-dep-node
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-undirected
#'
#' @concept undirected
#' @concept dyad-independent
#' @concept categorical nodal attribute
InitErgmTerm.sociality<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attrname", "base", "levels"),
                        vartypes = c("character", "numeric", "character,numeric,logical"),
                        defaultvalues = list(NULL, 1, NULL),
                        required = c(FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, "nodes", FALSE),
                        dep.warn = list(TRUE, FALSE, TRUE))
    attrarg <- a$attrname    
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attr", "base", "levels", "nodes"),
                        vartypes = c(ERGM_VATTR_SPEC, "numeric", ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, 1, NULL, LEVELS_BASE1),
                        required = c(FALSE, FALSE, FALSE, FALSE),
                        dep.inform = list(FALSE, "nodes", FALSE, FALSE),
                        dep.warn = list(TRUE, FALSE, TRUE, FALSE))
    attrarg <- a$attr
    levels <- a$levels
  }
  
  d <- ergm_attr_levels(a$nodes, 1:network.size(nw), nw, 1:network.size(nw))
  if((!hasName(attr(a,"missing"), "nodes") || attr(a,"missing")["nodes"]) && any(NVL(a$base,0)!=0)) d <- d[-a$base]
  
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
  }
  ld<-length(d)
  if(ld==0){return(NULL)}
  if(!is.null(attrarg)){
    coef.names <- paste("sociality",d,".",attrname,sep="")
    inputs <- c(d, 0, nodecov) # Input requires a "guard" value.
  }else{
    coef.names <- paste("sociality",d,sep="")
    inputs <- c(d,0) # Input requires a "guard" value.
  }
  list(name="sociality", coef.names=coef.names, inputs=inputs, minval=0, maxval=network.size(nw)-1, conflicts.constraints="degrees", dependence=FALSE)
}

#=======================InitErgmTerm functions:  T============================#


################################################################################

#' @templateVar name threetrail
#' @title Three-trails
#' @description For an undirected network, this term adds one statistic equal to the number
#'   of 3-trails, where a 3-trail is defined as a trail of length three that
#'   traverses three distinct edges.
#'   Note that a 3-trail need not
#'   include four distinct nodes; in particular, a triangle counts as three
#'   3-trails. For a directed network, this term adds four statistics
#'   (or some subset of these four),
#'   one for each of the four distinct types of directed three-paths. If the
#'   nodes of the path are written from left to right such that the middle edge
#'   points to the right (R), then the four types are RRR, RRL, LRR, and LRL.
#'   That is, an RRR 3-trail is of the form
#'   \eqn{i\rightarrow j\rightarrow k\rightarrow l}{i-->j-->k-->l} , and RRL
#'   3-trail is of the form
#'   \eqn{i\rightarrow j\rightarrow k\leftarrow l}{i-->j-->k<--l} , etc.
#'   Like in the undirected case, there is no requirement that the nodes be
#'   distinct in a directed 3-trail. However, the three edges must all be
#'   distinct. Thus, a mutual tie \eqn{i\leftrightarrow j}{i<-->j} does not
#'   count as a 3-trail of the form
#'   \eqn{i\rightarrow j\rightarrow i\leftarrow j}{i-->j-->i<--j} ; however,
#'   in the subnetwork \eqn{i\leftrightarrow j \rightarrow k}{i<-->j-->k} ,
#'   there are two directed 3-trails, one LRR
#'   ( \eqn{k\leftarrow j\rightarrow i\leftarrow j}{k<--j-->i-->j} )
#'   and one RRR
#'   ( \eqn{j\rightarrow i\rightarrow j\leftarrow k}{k<--j-->i-->j} ).
#'   
#' @usage
#' # binary: threetrail(keep=NULL, levels=NULL)
#'
#' @param keep deprecated
#' @templateVar explain specify a subset of the four statistics for directed networks.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-keep-dep
#'
#' @template ergmTerm-general
#'
#' @note This term used to be (inaccurately) called `threepath` . That
#'   name has been deprecated and may be removed in a future version.
#'
#' @concept directed
#' @concept undirected
#' @concept triad-related
InitErgmTerm.threetrail <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, 
                         varnames = c("keep"),
                         vartypes = c("numeric"),
                         defaultvalues = list(NULL),
                         required = c(FALSE),
                         dep.inform = list("levels"))
  }else{
    a <- check.ErgmTerm (nw, arglist, 
                         varnames = c("keep", "levels"),
                         vartypes = c("numeric", ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL),
                         required = c(FALSE, FALSE),
                         dep.inform = list("levels", FALSE))
  }  
  vals = c("RRR","RRL","LRR","LRL")
  types <- ergm_attr_levels(a$levels, vals, nw, levels = vals)
  if((!hasName(attr(a,"missing"), "levels") || attr(a,"missing")["levels"]) && !is.null(a$keep)) types <- types[a$keep]
  indices = match(types, vals)
  if (is.directed(nw)) {
    return(list(name = "threetrail", 
                coef.names = paste("threetrail", types, sep="."),
                inputs=indices, minval = 0))
  }
  else {
    return(list(name = "threetrail", coef.names = "threetrail", minval = 0))
  }
}

#' @templateVar name threetrail
#' @template ergmTerm-rdname
#' @aliases threepath-ergmTerm
#' @usage
#' # binary: threepath(keep=NULL, levels=NULL)
InitErgmTerm.threepath <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  ergm_Init_warn(paste("This term is inaccurately named and actually refers to a '3-trail' in that it counts repeated vertices: i-j-k-i is a 3-trail but not a 3-path. See", sQuote("ergmTerm?treepath"), "help for more information. This name has been deprecated and will be removed in a future version: if a 3-trail is what you want, use the term 'threetrail'."))

  f <- InitErgmTerm.threetrail
  f(nw, arglist, ..., version=version)
}

################################################################################

#' @templateVar name transitive
#' @title Transitive triads
#' @description This term adds one statistic to the model, equal to the number of triads in
#'   the network that are transitive. The transitive triads are those of type
#'   `120D` , `030T` , `120U` , or `300` in the categorization
#'   of Davis and Leinhardt (1972). For details on the 16 possible triad types,
#'   see `?triad.classify` in the \CRANpkg{sna} package.
#'   Note the distinction from the `ttriple` term. This term can only be
#'   used with directed networks.
#'
#' @usage
#' # binary: transitive
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept triad-related
InitErgmTerm.transitive<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="transitive", coef.names="transitive", minval = 0)
}

################################################################################

#' @templateVar name triadcensus
#' @title Triad census
#' @description For a directed network, this term adds one network statistic for each of
#'   an arbitrary subset of the 16 possible types of triads categorized by
#'   Davis and Leinhardt (1972) as `003, 012, 102, 021D, 021U, 021C, 111D,
#'   ` `	111U, 030T, 030C, 201, 120D, 120U, 120C, 210,` and `300` . Note that at
#'   least one category should be dropped; otherwise a linear dependency will
#'   exist among the 16 statistics, since they must sum to the total number of
#'   three-node sets. By default, the category `003` , which is the category
#'   of completely empty three-node sets, is dropped. This is considered category
#'   zero, and the others are numbered 1 through 15 in the order given above. Each statistic is the count of the corresponding triad
#'   type in the network. For details on the 16 types, see `?triad.classify`
#'   in the \CRANpkg{sna} package, on which this code is based. For an undirected
#'   network, the triad census is over the four types defined by the number of
#'   ties (i.e., 0, 1, 2, and 3).
#'
#' @usage
#' # binary: triadcensus(levels)
#'
#' @templateVar explain For directed networks, specify a set of terms to add other than the default value of `1:15`. 
#'   For undirected networks, specify which of the four types of ties to include. The default is `1:3` where 0 is dropped.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept triad-related
#' @concept directed
#' @concept undirected
InitErgmTerm.triadcensus<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("d"),
                        vartypes = c("numeric"),
                        defaultvalues = list(NULL),
                        required = c(FALSE))
    d <- a$d                    
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("levels"),
                        vartypes = c(ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL),
                        required = c(FALSE))
    d <- a$levels                      
  }

  emptynwstats<-NULL

  if(is.directed(nw)){
   tcn <- c("003","012", "102", "021D", "021U", "021C", "111D",
            "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
  }else{
#  Undirected
   tcn <- c("0", "1", "2", "3")
  }
  
  if(is.null(d)){
    d <- 1:(length(tcn) - 1)
  }
  
  if(is.numeric(d)){
    d <- d + 1
  }
  
  d <- ergm_attr_levels(d, tcn, nw, levels = tcn)
  
  d <- match(d, tcn) - 1
  
  if (any(d==0)) {
    emptynwstats <- rep(0,length(d))
    nwsize <- network.size(nw)
    # SEARCH_ON_THIS_TO_TRACK_DOWN_TRIADCENSUS_CHANGE
    # to undo triadcensus change, comment out next line:
    emptynwstats[d==0] <- nwsize * (nwsize-1) * (nwsize-2) / 6
  }
  d <- d + 1
  lengthd<-length(d)
  if(lengthd==0){return(NULL)}
  # No covariates here, so "ParamsBeforeCov" unnecessary
  coef.names <- paste("triadcensus",tcn,sep=".")[d]
  if (!is.null(emptynwstats)){
    list(name="triadcensus", coef.names=coef.names, inputs=c(d),
         emptynwstats=emptynwstats, dependence=TRUE)
  }else{
    list(name="triadcensus", coef.names=coef.names, inputs=c(d),
         dependence=TRUE, minval = 0)
  }
}



################################################################################

#' @templateVar name triangle
#' @aliases triangles-ergmTerm
#' @title Triangles
#' @description By default, this term adds one statistic to the model equal to the number of triangles
#'   in the network. For an undirected network, a triangle is defined to be any
#'   set \eqn{\{(i,j), (j,k), (k,i)\}} of three edges. For a directed network, a
#'   triangle is defined as any set of three edges \eqn{(i{\rightarrow}j)}{(i,j)}
#'   and \eqn{(j{\rightarrow}k)}{(j,k)} and either \eqn{(k{\rightarrow}i)}{(k,i)}
#'   or \eqn{(k{\leftarrow}i)}{(i,k)} . The former case is called a "transitive
#'   triple" and the latter is called a "cyclic triple", so in the case of a
#'   directed network, `triangle` equals `ttriple` plus `ctriple`
#'   --- thus at most two of these three terms can be in a model. 
#'
#' @usage
#' # binary: triangle(attr=NULL, diff=FALSE, levels=NULL)
#'
#' # binary: triangles(attr=NULL, diff=FALSE, levels=NULL)
#'
#' @param attr,diff quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If `attr` is specified and `diff` is `FALSE` ,
#'   then the count is restricted to those triples of nodes with
#'   equal values of the vertex attribute specified by `attr` . If `attr` is specified and `diff` is `TRUE` ,
#'   then one statistic is added for each value of `attr` ,
#'   equal to the number of triangles where all three nodes have that value of the attribute.
#' @templateVar explain add one statistic for each value specified if `diff` is `TRUE`.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept frequently-used
#' @concept triad-related
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.triangle<-InitErgmTerm.triangles<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attrname", "diff", "levels"),
                        vartypes = c("character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL    
  }else{
    a <- check.ErgmTerm(nw, arglist,
                        varnames = c("attr", "diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels      
  }

  diff <- a$diff
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (!diff) {
      coef.names <- paste("triangle",attrname,sep=".")
      inputs <- c(nodecov)
    } else {
      coef.names <- paste("triangle",attrname, u, sep=".")
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)
    }
  }else{
    coef.names <- "triangle"
    inputs <- NULL
  }
  list(name="triangle", coef.names=coef.names, inputs=inputs, minval=0)
}



################################################################################

#' @templateVar name tripercent
#' @title Triangle percentage
#' @description By default, this term adds one statistic to the model equal to 100 times the ratio of
#'   the number of triangles in the network to the sum of the number of triangles
#'   and the number of 2-stars not in triangles (the latter is considered a
#'   potential but incomplete triangle). In case the denominator equals zero,
#'   the statistic is defined to be zero. For the definition of triangle, see
#'   `triangle` . This is often called
#'   the mean correlation coefficient. This term can only be
#'   used with undirected networks; for directed networks, it is difficult to
#'   define the numerator and denominator in a consistent and meaningful way.
#'
#' @usage
#' # binary: tripercent(attr=NULL, diff=FALSE, levels=NULL)
#'
#' @param attr,diff quantitative attribute (see Specifying Vertex attributes and Levels (`?nodal_attributes`) for details.) If `attr` is specified and `diff` is `FALSE` ,
#'   then the counts are restricted to those triples of nodes with
#'   equal values of the vertex attribute specified by `attr` . If `attr` is specified and `diff` is `TRUE` ,
#'   then one statistic is added for each value of `attr` ,
#'   equal to the number of triangles where all three nodes have that value of the attribute.
#' @templateVar explain add one statistic for each value specified if `diff` is `TRUE`
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @concept undirected
#' @concept triad-related
#' @concept categorical nodal attribute
InitErgmTerm.tripercent<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attrname", "diff", "levels"),
                        vartypes = c("character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                            
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("attr", "diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels
  }  
  diff <- a$diff
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (!diff) {
      coef.names <- paste("tripercent",attrname,sep=".")
      inputs <- c(1, nodecov)
      attr(inputs, "ParamsBeforeCov") <- 1
    } else {
      coef.names <- paste("tripercent",attrname, u, sep=".")
      inputs <- c(ui, nodecov)
      attr(inputs, "ParamsBeforeCov") <- length(ui)
    }
  }else{
    coef.names <- "tripercent"
    inputs <- NULL
  }
  list(name="tripercent", coef.names=coef.names, inputs=inputs, minval = 0)
}


################################################################################

#' @templateVar name ttriple
#' @title Transitive triples
#' @description By default, this term adds one statistic to the model, equal to the number of transitive
#'   triples in the network, defined as a set of edges \eqn{\{(i{\rightarrow}j),
#'   j{\rightarrow}k), (i{\rightarrow}k)\}}{\{(i,j), (j,k), (i,k)\}} . Note that
#'   `triangle` equals `ttriple+ctriple` for a directed network, so at
#'   most two of the three terms can be in a model. 
#'
#' @usage
#' # binary: ttriple(attr=NULL, diff=FALSE, levels=NULL)
#'
#' @template ergmTerm-attr
#' @param diff If `attr` is specified and `diff` is `FALSE` ,
#'   then the count is over the number of transitive triples where all three nodes have the same value of
#'   the attribute. If `attr` is specified and `diff` is `TRUE` ,
#'   then one statistic is added for each value of `attr` ,
#'   equal to the number of triangles where all three nodes have that value of the attribute.
#' @templateVar explain add one statistic for each value specified if `diff` is `TRUE`.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @concept directed
#' @concept triad-related
#' @concept categorical nodal attribute
InitErgmTerm.ttriple<-function (nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attrname", "diff", "levels"),
                        vartypes = c("character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attrname
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("attr", "diff", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, FALSE, NULL),
                        required = c(FALSE, FALSE, FALSE))
    attrarg <- a$attr
    levels <- a$levels  
  }
  diff <- a$diff
  if(!is.null(attrarg)) {
    nodecov <- ergm_get_vattr(attrarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (!diff) {
      coef.names <- paste("ttriple",attrname,sep=".")
      inputs <- c(nodecov)
     } else {
       coef.names <- paste("ttriple",attrname, u, sep=".")
       inputs <- c(ui, nodecov)
       attr(inputs, "ParamsBeforeCov") <- length(ui)
     }
  }else{
    coef.names <- "ttriple"
    inputs <- NULL
  }
  list(name="ttriple", coef.names=coef.names, inputs=inputs, minval = 0)
}

#' @templateVar name ttriple
#' @template ergmTerm-rdname
#' @aliases ttriad-ergmTerm
#' @usage
#' # binary: ttriad
InitErgmTerm.ttriad<-InitErgmTerm.ttriple

################################################################################

#' @templateVar name twopath
#' @title 2-Paths
#' @description This term adds one statistic to the model, equal to the number of 2-paths in
#'   the network. For a directed network this is defined as a pair of edges
#'   \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)} , where \eqn{i} and
#'   \eqn{j} must be distinct. That is, it is a directed path of length 2 from
#'   \eqn{i} to \eqn{k} via \eqn{j} . For directed networks a 2-path is also a
#'   mixed 2-star but the interpretation is usually different; see `m2star` .
#'   For undirected networks a twopath is defined as a pair of edges
#'   \eqn{\{i,j\}, \{j,k\}} . That is, it is an undirected path of length 2 from
#'   \eqn{i} to \eqn{k} via \eqn{j} , also known as a 2-star.
#'
#' @usage
#' # binary: twopath
#'
#' @template ergmTerm-general
#'
#' @concept directed
#' @concept undirected
InitErgmTerm.twopath<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  if(is.directed(nw)){
   list(name="m2star", coef.names="twopath", dependence=TRUE, minval=0)
  }else{
   k<-2
   lk<-length(k)
   if(lk==0){return(NULL)}
   list(name="kstar", coef.names="twopath", inputs=c(k), minval=0)
  }
}


