#  File R/InitErgmTerm.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#NOTE: a number of undocumented terms have been removed from this file
# the terms still exist on the experimental_terms svn branch


#===========================================================================
# This file contains the following 74 new, easier-to-write ergm-term
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

GWDECAY <- list(
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

# LEVELS_BASE1 is a placeholder for whatever the value of levels= or
# nodes= should be when base==1. For now, it's NULL to prevent the two
# arguments from interfering. Eventually, when base= is removed, it
# will need to be set to -1 either here or by search-and-replace.
LEVELS_BASE1 <- NULL

#=======================InitErgmTerm functions:  A============================#


################################################################################
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
InitErgmTerm.b1degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }

  ### Process the arguments
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) ergm_Init_abort("Term odegrange must have from<to.")

  nb1 <- get.network.attribute(nw, "bipartite")
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw, bip = if(homophily) "n" else "b1")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine b1degrange and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
                         paste("b1deg",from,"+",sep=""),
                         paste("b1deg",from,"to",to,sep=""))
    name <- "b1degrange"
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in d_b1degrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("b1deg",from,"+", ".homophily.",attrname,sep=""),
                         paste("b1deg",from,"to",to, ".homophily.",attrname,sep=""))
    name <- "b1degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_b1degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("b1deg",du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste("b1deg",du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- "b1degrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name,coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval = 0)
  }else{
    list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="b1degreedist")
  }
}

################################################################################
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
InitErgmTerm.b1degree <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("d", "by", "levels"),
                         vartypes = c("numeric", "character", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("d", "by", "levels"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    levels <- a$levels    
  }
  
  byarg <- a$by
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if (!is.null(byarg)) {  # CASE 1:  a$by GIVEN
    nodecov <- ergm_get_vattr(byarg, nw, bip = "b1")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b1degree_by_attr"
    coef.names <- paste("b1deg", du[1,], ".", attrname, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$by NOT GIVEN
    name <- "b1degree"
    coef.names <- paste("b1deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- nb1
    }
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=nb1, conflicts.constraints="odegreedist")
}


################################################################################
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
InitErgmTerm.b2degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }

  ### Process the arguments  
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) ergm_Init_abort("Term odegrange must have from<to.")

  nb1 <- get.network.attribute(nw, "bipartite")
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw, bip = if(homophily) "n" else "b2")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine b2degrange and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
                         paste("b2deg",from,"+",sep=""),
                         paste("b2deg",from,"to",to,sep=""))
    name <- "b2degrange"
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in d_b2degrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("b2deg",from,"+", ".homophily.",attrname,sep=""),
                         paste("b2deg",from,"to",to, ".homophily.",attrname,sep=""))
    name <- "b2degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_b2degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("b2deg",du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste("b2deg",du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- "b2degrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name,coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval = 0)
  }else{
    list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="b2degreedist")
  }
}



################################################################################
InitErgmTerm.b2degree <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("d", "by", "levels"),
                         vartypes = c("numeric", "character", "character,numeric,logical"),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL
  }else{
    ### Check the network and arguments to make sure they are appropriate.
    a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                         varnames = c("d", "by", "levels"),
                         vartypes = c("numeric", ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                         defaultvalues = list(NULL, NULL, NULL),
                         required = c(TRUE, FALSE, FALSE))
    levels <- a$levels    
  }
  
  byarg <- a$by
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if (!is.null(byarg)) {  # CASE 1:  a$by GIVEN
    nodecov <- ergm_get_vattr(byarg, nw, bip = "b2")
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b2degree_by_attr"
    coef.names <- paste("b2deg", du[1,], ".", attrname, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$by NOT GIVEN
    name <- "b2degree"
    coef.names <- paste("b2deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- n-nb1
    }
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=network.size(nw)-nb1, conflicts.constraints="b2degreedist")
}

################################################################################
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
InitErgmTerm.balance<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist)
    
  list(name="balance", coef.names="balance", dependence=TRUE, minval=0)
}


#=======================InitErgmTerm functions:  C============================#


################################################################################
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
InitErgmTerm.ctriple<-InitErgmTerm.ctriad<-function (nw, arglist, ..., version=packageVersion("ergm")) {
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



################################################################################
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
InitErgmTerm.degrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) ergm_Init_abort("Term degrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
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
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
                         paste("deg",from,"+",sep=""),
                         paste("deg",from,"to",to,sep=""))
    name <- "degrange"
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in d_degrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("deg",from,"+", ".homophily.",attrname,sep=""),
                         paste("deg",from,"to",to, ".homophily.",attrname,sep=""))
    name <- "degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("deg",du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste("deg",du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- "degrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name,coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval = 0)
  }else{
    list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="degreedist")
  }
}



################################################################################
InitErgmTerm.degree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
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
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
    coef.names <- paste("degree",d,sep="")
    name <- "degree"
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    # See comment in d_degree_w_homophily function
    coef.names <- paste("deg", d, ".homophily.",attrname, sep="")
    name <- "degree_w_homophily"
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", attrname,u[du[2,]], sep="")
    name <- "degree_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=network.size(nw), conflicts.constraints="degreedist")
}


################################################################################
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
#' @describeIn ergm-deprecated Use [`degree1.5`] instead.
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
InitErgmTerm.density<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="density", coef.names="density", dependence=FALSE, minval = 0, maxval = 1, conflicts.constraints="edges")
}

################################################################################
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

  if(sign.action!="abs" && !is.directed(nw) && !is.bipartite(nw)) ergm_Init_inform("Note that behavior of term diff() on unipartite, undirected networks may be unexpected. See help(\"ergm-terms\") for more information.")
  
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
InitErgmTerm.dyadcov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x","attrname"),
                      vartypes = c("matrix,network,character","character"),
                      defaultvalues = list(NULL,NULL),
                      required = c(TRUE,FALSE))
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix.network(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x)){
    xm<-get.network.attribute(nw,a$x)
    if (is.null(xm)){
      ergm_Init_abort("There is no network attribute named ",a$x)
    }
  }
  else
    xm<-as.matrix(a$x)

#Update the terms list, adding the vectorized adjacency matrix
  if(!is.null(a$attrname))
    cn<-paste("dyadcov", as.character(sys.call(0)[[3]][2]), 
              as.character(a$attrname), sep = ".")
  else
    cn<-paste("dyadcov", as.character(sys.call(0)[[3]][2]), sep = ".")
 
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
InitErgmTerm.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network,character", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix.network(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x)){
    xm<-get.network.attribute(nw,a$x)
    if (is.null(xm)){
      ergm_Init_abort("There is no network attribute named ",a$x)
    }
  }
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov", coef.names = cn, inputs = inputs, dependence=FALSE,
       minval = sum(c(xm)[c(xm)<0]),
       maxval = sum(c(xm)[c(xm)>0])
       )
}

################################################################################
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
  decay<-a$decay; fixed<-a$fixed
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:(network.size(nw) - nb1)
  maxesp <- min(cutoff, network.size(nw)-nb1)

  d <- 1:maxesp
  if (!is.null(attrarg) && !fixed) {
    ergm_Init_abort("The gwb1degree term cannot yet handle a nonfixed decay ",
                    "term with an attribute. Use fixed=TRUE.")
  }
  if(!fixed){# This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwb1degree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="b1degree", coef.names=paste("gwb1degree#",d,sep=""), inputs=c(d),
           conflicts.constraints="b1degreedist", params=list(gwb1degree=NULL,gwb1degree.decay=decay)), GWDECAY)
  } else {
    if(is.null(a$decay)) stop("Term 'gwb1degree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

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
InitErgmTerm.gwb2degree<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ..., version=packageVersion("ergm")) {
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
  decay<-a$decay; fixed<-a$fixed
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:nb1
  maxesp <- min(cutoff,nb1)
  d <- 1:maxesp
  if (!is.null(attrarg) && !fixed) {
      ergm_Init_abort("The gwb2degree term cannot yet handle a nonfixed decay ",
                      "term with an attribute. Use fixed=TRUE.") }
  if(!fixed){# This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwb2degree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="b2degree", coef.names=paste("gwb2degree#",d,sep=""), inputs=c(d),
           conflicts.constraints="b2degreedist", params=list(gwb2degree=NULL,gwb2degree.decay=decay)), GWDECAY)
  } else { 
    if(is.null(a$decay)) stop("Term 'gwb2degree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

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
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrarg) && !fixed) {
    ergm_Init_abort("The gwdegree term cannot yet handle a nonfixed decay ",
                        "term with an attribute. Use fixed=TRUE.")
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwdegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="degree", coef.names=paste("gwdegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="degreedist", params=list(gwdegree=NULL,gwdegree.decay=decay)), GWDECAY)
  } else {
    if(is.null(a$decay)) stop("Term 'gwdegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

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
InitErgmTerm.gwdsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","alpha"),
                      vartypes = c("numeric","logical","numeric","numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    ergm_Init_abort("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".")
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwdsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
#   d <- 1:(network.size(nw)-1)
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
    c(list(name=dname, coef.names=paste("gwdsp#",d,sep=""), 
           inputs=c(d), params=list(gwdsp=NULL,gwdsp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    if(is.null(a$decay)) stop("Term 'gwdsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if (!fixed) # First pass to get MPLE coefficient
      coef.names <- "gwdsp"   # must match params$gwdsp above
    else  # fixed == TRUE
      coef.names <- paste("gwdsp.fixed.",decay,sep="")
  if(is.directed(nw)){dname <- "gwtdsp"}else{dname <- "gwdsp"}
  list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}



################################################################################
InitErgmTerm.gwesp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    ergm_Init_abort("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".")
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwesp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

    #   d <- 1:(network.size(nw)-2)
     maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
    c(list(name=dname, coef.names=paste("esp#",d,sep=""), 
         inputs=c(d), params=list(gwesp=NULL,gwesp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    if(is.null(a$decay)) stop("Term 'gwesp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    coef.names <- paste("gwesp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtesp"}else{dname <- "gwesp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)
  }
}



################################################################################
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
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
  maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrarg) && !fixed ) {
    ergm_Init_abort("The gwidegree term cannot yet handle a nonfixed decay ",
                        "term with an attribute. Use fixed=TRUE.")
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwidegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="idegree", coef.names=paste("gwidegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="idegreedist", params=list(gwidegree=NULL,gwidegree.decay=decay)), GWDECAY)
  } else { 
    if(is.null(a$decay)) stop("Term 'gwidegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

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
InitErgmTerm.gwnsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    ergm_Init_abort("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".")
  }

  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwnsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

#   d <- 1:(network.size(nw)-1)
     maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    if(is.directed(nw)){dname <- "tnsp"}else{dname <- "nsp"}
    c(list(name=dname, coef.names=paste("nsp#",d,sep=""),
           inputs=c(d), params=list(gwnsp=NULL,gwnsp.decay=decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL),
      GWDECAY)
  }else{
    if(is.null(decay)) stop("Term 'gwnsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    coef.names <- paste("gwnsp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtnsp"}else{dname <- "gwnsp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay), auxiliaries=if(cache.sp) .spcache.aux(if(is.directed(nw)) "OTP" else "UTP") else NULL)    
  }
}


################################################################################
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
  decay<-a$decay; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrarg) && !fixed ) {
    ergm_Init_abort("The gwodegree term cannot yet handle a nonfixed decay ",
                        "term with an attribute. Use fixed=TRUE.")
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwodegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(list(minval=0, maxval=network.size(nw), dependence=TRUE, name="odegree", coef.names=paste("gwodegree#",d,sep=""), inputs=c(d),
           conflicts.constraints="odegreedist", params=list(gwodegree=NULL,gwodegree.decay=decay)), GWDECAY)
  } else {
    if(is.null(a$decay)) stop("Term 'gwodegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

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
InitErgmTerm.idegrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term idegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) ergm_Init_abort("Term idegrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine idegrange and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
                         paste("ideg",from,"+",sep=""),
                         paste("ideg",from,"to",to,sep=""))
    name <- "idegrange"
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in d_idegrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("ideg",from,"+", ".homophily.",attrname,sep=""),
                         paste("ideg",from,"to",to, ".homophily.",attrname,sep=""))
    name <- "idegrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_idegrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("ideg",du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste("ideg",du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- "idegrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name,coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval = 0)
  }else{
    list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="idegreedist")
  }
}

################################################################################
InitErgmTerm.idegree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
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
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
    name <- "idegree"
    coef.names <- paste("idegree",d,sep="")
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    name <- "idegree_w_homophily"
    coef.names <- paste("ideg", d, ".homophily.",attrname, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "idegree_by_attr"
    # See comment in d_idegree_by_attr function
    coef.names <- paste("ideg", du[1,], ".", attrname,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=network.size(nw), conflicts.constraints="idegreedist")
}



################################################################################
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
#' @describeIn ergm-deprecated Use [`idegree1.5`] instead.
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
InitErgmTerm.intransitive<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="intransitive", coef.names="intransitive", minval = 0)
}

################################################################################
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
InitErgmTerm.m2star<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="m2star", coef.names="m2star",dependence=TRUE, minval = 0) 
}


################################################################################
InitErgmTerm.meandeg<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="meandeg", coef.names="meandeg", dependence=FALSE, minval=0, maxval=if(!is.bipartite(nw)) network.size(nw)-1, conflicts.constraints="edges")
}


################################################################################
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
InitErgmTerm.nearsimmelian<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nearsimmelian", coef.names="nearsimmelian", minval=0, maxval=network.dyadcount(nw,FALSE)*network.size(nw)*0.5)
}


################################################################################
InitErgmTerm.nodecov<-InitErgmTerm.nodemain<-function (nw, arglist, ..., version=packageVersion("ergm")) {
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



################################################################################
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
   
    levels2.sel <- if((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) levels2.list[-a$base]
                   else ergm_attr_levels(a$levels2, list(row = b1nodecov, col = b2nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
  
    # Recode to numeric
    b1nodecov <- match(b1nodecov,b1namescov,nomatch=length(b1namescov)+1)
    b2nodecov <- match(b2nodecov,b2namescov,nomatch=length(b2namescov)+1)
  
    namescov <- c(b1namescov, b2namescov)
    
    nodecov <- c(b1nodecov, b2nodecov)
    
    cn <- paste("mix", paste(attrname,collapse="."), apply(matrix(namescov[as.matrix(u)],ncol=2),
                                       1,paste,collapse="."), sep=".")
                                       
    ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
    indmat <- matrix(0L, nrow = nr + 1, ncol = nc + 1)
    u[,2L] <- u[,2L] - nr
    indmat[as.matrix(u)] <- seq_len(NROW(u))
    indmat <- indmat - 1L
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
   
    levels2.sel <- if((!hasName(attr(a,"missing"), "levels2") || attr(a,"missing")["levels2"]) && any(NVL(a$base,0)!=0)) levels2.list[-a$base]
                   else ergm_attr_levels(a$levels2, list(row = nodecov, col = nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
    uun <- uun[rows2keep]

    nodecov <- match(nodecov,namescov,nomatch=length(namescov)+1)
    
    cn <- paste("mix", paste(attrname,collapse="."), uun, sep=".")

    ## the +1 for nrow and ncol are needed by the way we code non-included b1 and b2 levels above
    indmat <- matrix(0L, nrow = nr + 1, ncol = nc + 1)
    indmat[as.matrix(u)] <- seq_len(NROW(u))
    if(!is.directed(nw)) indmat <- indmat + t(indmat) - diag(diag(indmat))
    indmat <- indmat - 1L
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
InitErgmTerm.odegrange<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                        
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("from", "to", "by", "homophily", "levels"),
                        vartypes = c("numeric", "numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, Inf, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) ergm_Init_abort("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) ergm_Init_abort("Term odegrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
    attrname <- attr(nodecov, "name")
    u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))
    nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
  }
  if(!is.null(byarg) && !homophily) {
    # Combine odegrange and u into 3xk matrix, where k=length(from)*length(u)
    lu <- length(u)
    du <- rbind(rep(from,lu), rep(to,lu), rep(1:lu, rep(length(from), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[3,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
                         paste("odeg",from,"+",sep=""),
                         paste("odeg",from,"to",to,sep=""))
    name <- "odegrange"
    inputs <- c(rbind(from,to))
  } else if (homophily) {
    if(length(from)==0){return(NULL)}
    # See comment in d_odegrange_w_homophily function
    coef.names <- ifelse(to>=network.size(nw)+1,
                         paste("odeg",from,"+", ".homophily.",attrname,sep=""),
                         paste("odeg",from,"to",to, ".homophily.",attrname,sep=""))
    name <- "odegrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_odegrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("odeg",du[1,],"+.", attrname, u[du[3,]],sep=""),
                         paste("odeg",du[1,],"to",du[2,],".",attrname, u[du[3,]],sep=""))
    name <- "odegrange_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name,coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval = 0)
  }else{
    list(name=name,coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="odegreedist")
  }
}

################################################################################
InitErgmTerm.odegree<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", "character", "logical", "character,numeric,logical"),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- if(!is.null(a$levels)) I(a$levels) else NULL                                                
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                        varnames = c("d", "by", "homophily", "levels"),
                        vartypes = c("numeric", ERGM_VATTR_SPEC, "logical", ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL, FALSE, NULL),
                        required = c(TRUE, FALSE, FALSE, FALSE))
    levels <- a$levels  
  }
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- ergm_get_vattr(byarg, nw)
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
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
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
    name <- "odegree"
    coef.names <- paste("odegree",d,sep="")
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    name <- "odegree_w_homophily"
    # See comment in d_odegree_w_homophily function
    coef.names <- paste("odeg", d, ".homophily.",attrname, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "odegree_by_attr"
    # See comment in d_odegree_by_attr function
    coef.names <- paste("odeg", du[1,], ".", attrname,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }

  list(name = name, coef.names = coef.names, inputs = inputs, emptynwstats = emptynwstats, minval=0, maxval=network.size(nw), dependence=TRUE,
    minval = 0, maxval=network.size(nw), conflicts.constraints="odegreedist")
}



################################################################################
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
#' @describeIn ergm-deprecated Use [`odegree1.5`] instead.
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
InitErgmTerm.opentriad<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  list(name="opentriad", coef.names="opentriad", inputs=NULL)
}


################################################################################
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
InitErgmTerm.simmelian<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="simmelian", coef.names="simmelian", minval=0, maxval=network.edgecount(nw)*network.size(nw)*0.5)
}


################################################################################
InitErgmTerm.simmelianties<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="simmelianties", coef.names="simmelianties", minval=0, maxval=network.edgecount(nw)) # TODO: Is this correct?
}



################################################################################
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
InitErgmTerm.threepath <- function(nw, arglist, ..., version=packageVersion("ergm")) {
  ergm_Init_warn("This term is inaccurately named and actually refers to a '3-trail' in that it counts repeated vertices: i-j-k-i is a 3-trail but not a 3-path. See ergm-terms help for more information. This name has been deprecated and will be removed in a future version: if a 3-trail is what you want, use the term 'threetrail'.")
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

################################################################################
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

################################################################################
InitErgmTerm.transitive<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="transitive", coef.names="transitive", minval = 0)
}

################################################################################
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
InitErgmTerm.ttriple<-InitErgmTerm.ttriad<-function (nw, arglist, ..., version=packageVersion("ergm")) {
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



################################################################################
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


