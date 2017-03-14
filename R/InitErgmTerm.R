#  File R/InitErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

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
#   H:   <hamming>          <hammingmix>
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
# built via <ergm.getmodel>
# 
# --PARAMETERS--
#   nw        : the network given in formula F
#   arglist   : the arguments given with term X in formula F
#   initialfit: whether the parameters for this term have been initially fit
#               (T or F); if FALSE, the ergm belongs to the curved exponential
#               family; if TRUE, the term X does not the ergm to be curved,
#               though other terms may; default=FALSE
#
# --IGNORED PARAMETERS--
#   ... : ignored, but necessary to accomodate other arguments
#         passed by <ergm.getmodel>
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
#    params      : a list of parameter values for curved exponential family model
#                  terms only; each item in the list should be named with the
#                  corresponding parameter name; those that coincide with the
#                  coef.names (used when initialfit=TRUE) will have their 'params'
#                  set by MPLE and their initial values in 'params' are ignored;
#                  otherwise, parameters should be given an initial value in this list
#    map         : a function taking two arguments, theta and length('params'), which
#                  gives the map from the canonical parameters, theta, to the curved
#                  parameters, eta; 'map' is only necessary for curved exponential
#                  family model terms
#   gradient     : a function taking two arguments, theta and length('params'), which
#                  gives the gradient of the eta map above as a p by q matrix, where
#                  p=length(theta), q=length(params); 'gradient' is only necessary
#                  for curved exponential family model terms
#
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




#=======================InitErgmTerm functions:  A============================#


################################################################################
InitErgmTerm.absdiff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list(NULL,1),
                      required = c(TRUE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  ### Construct the list to return
  list(name="absdiff",                                     #name: required
       coef.names = paste(paste("absdiff",if(a$pow!=1) a$pow else "",sep=""), a$attrname, sep="."), #coef.names: required
       inputs = c(a$pow,nodecov),  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}



################################################################################
InitErgmTerm.absdiffcat <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","base"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list(NULL,NULL),
                      required = c(TRUE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=TRUE)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(a$base)) u <- u[-(a$base)]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  u2 <- u[!is.na(u)]
  ### Construct the list to return
  inputs <- c(u2, NAsubstitute, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(u2)+1 # See comment at top of file
  list(name="absdiffcat",                                  #name: required
       coef.names = paste("absdiff", a$attrname, u, sep="."), #coef.names: required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}






################################################################################
InitErgmTerm.altkstar <- function(nw, arglist, initialfit=FALSE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=NULL,
                      varnames = c("lambda","fixed"),
                      vartypes = c("numeric","logical"),
                      defaultvalues = list(1,FALSE),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  if(!initialfit && !a$fixed){ # This is a curved exponential family model
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
       params=list(altkstar=NULL, altkstar.lambda=a$lambda)
       )
  } else {
    if (initialfit) { # coef.names must match "altkstar" from params list above
      coef.names = "altkstar"
    } else { # Not necessary to match; provide more informative label
      coef.names = paste("altkstar", a$lambda, sep=".")
    }
    outlist <- list (name="altkstar",                      #name: required
                     coef.names = coef.names,
                     inputs=a$lambda
                     )
  }
  outlist
}



################################################################################
InitErgmTerm.asymmetric <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attrname", "diff", "keep"),
                      vartypes = c("character", "logical", "numeric"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE))
  ### Process the arguments
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u <- sort(unique(nodecov))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
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
  if (!is.null(a$attrname)) {
    if (a$diff) {
      out$coef.names <- paste("asymmetric", a$attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("asymmetric", a$attrname, sep=".")
      out$inputs <- nodecov
    }
  }

  out
}


#=======================InitErgmTerm functions:  B============================#


################################################################################
InitErgmTerm.b1concurrent<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("by"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  byarg <- a$by
  nb1 <- get.network.attribute(nw, "bipartite")       
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "b1concurrent")[seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b1concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_b1concurrent_by_attr function
    name <- "b1concurrent_by_attr"
    coef.names<-paste("b1concurrent",".", byarg, u, sep="")
    inputs <- c(ui, nodecov)
  }else{
    name <- "b1concurrent"
    coef.names<-paste("b1concurrent",sep="")
    inputs <- NULL
  }
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval=0, maxval=nb1)
}

################################################################################
InitErgmTerm.b1degrange<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                      varnames = c("from", "to", "by", "homophily"),
                      vartypes = c("numeric", "numeric", "character", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term odegrange must have from<to.")

  nb1 <- get.network.attribute(nw, "bipartite")
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "b1degrange")
    if(!homophily) nodecov <- nodecov[seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to b1degrange() has only one value", call.=FALSE)
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
                         paste("b1deg",from,"+", ".homophily.",byarg,sep=""),
                         paste("b1deg",from,"to",to, ".homophily.",byarg,sep=""))
    name <- "b1degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_b1degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("b1deg",du[1,],"+.", byarg, u[du[3,]],sep=""),
                         paste("b1deg",du[1,],"to",du[2,],".",byarg, u[du[3,]],sep=""))
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
InitErgmTerm.b1cov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE, 
                      varnames = c("attrname","transform","transformname"),
                      vartypes = c("character","function","character"),
                      defaultvalues = list(NULL,function(x)x,""),
                      required = c(TRUE,FALSE,FALSE))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste(paste("nodeocov",f.name,sep=""),attrname,sep=".")
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- f(get.node.attr(nw, attrname, "nodeocov", numeric=TRUE)[1:nb1])
  # C implementation is identical
  list(name="nodeocov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}


################################################################################
InitErgmTerm.b1degree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "by"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if (!is.null(a$by)) {  # CASE 1:  a$by GIVEN
    nodecov <- get.node.attr(nw, a$by)[seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
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
    coef.names <- paste("b1deg", du[1,], ".", a$by, u[du[2,]], sep="")
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
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats, minval=0, maxval=nb1)
}


################################################################################
InitErgmTerm.b1factor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))                                    
  attrname<-a$attrname
  base <- a$base
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, attrname, "b1factor")[1:nb1]
  
  if(all(is.na(nodecov)))
	  stop("Argument to b1factor() does not exist", call.=FALSE)
  
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to b1factor() has only one value", call.=FALSE)
  }
  if(base[1]==0){
    coef.names <- paste("b1factor", attrname, paste(u), sep=".")
    inputs <- c(ui, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lu
  }else{
    coef.names <- paste("b1factor",attrname, paste(u[-base]), sep=".")
    inputs <- c(ui[-base], nodecov)
    attr(inputs, "ParamsBeforeCov") <- lu-length(base)
  }
  list(name="b1factor", coef.names=coef.names, inputs=inputs, dependence=FALSE, minval=0)
}



################################################################################
InitErgmTerm.b1star <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "ostar"
    coef.names <- paste("b1star", a$k, ".", a$attrname, sep="")
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
InitErgmTerm.b1starmix <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname", "base", "diff"),
                       vartypes = c("numeric", "character", "numeric", "logical"),
                       defaultvalues = list(NULL, NULL, NULL, TRUE),
                       required = c(TRUE, TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, a$attrname)
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { stop("Only a single scalar k may be used with each b1starmix term") }
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
    coef.names <- paste("b1starmix", a$k, a$attrname,
                        apply(matrix(namescov[u],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- 1:nr
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base] }
    name <- "b1starmixhomophily"
    coef.names <- paste("b1starmix", a$k, a$attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval = 0)
}

################################################################################
InitErgmTerm.b1twostar <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("b1attrname", "b2attrname", "base"),
                       vartypes = c("character", "character", "numeric"),
                       defaultvalues = list(NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  b1nodecov <- get.node.attr(nw, a$b1attrname)[1:nb1]
  b1u<-sort(unique(b1nodecov))
  if(any(is.na(b1nodecov))){ b1u<-c(b1u,NA) }
  if(is.null(a$b2attrname)) { a$b2attrname <- a$b1attrname }
  b2nodecov <- get.node.attr(nw, a$b2attrname)[(1+nb1):n]
  b2u<-sort(unique(b2nodecov))
  if(any(is.na(b2nodecov))){b2u<-c(b2u,NA)}
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  nr <- length(b1u)
  nc <- length(b2u)
  u <- cbind(rep(1:nr, nc*nc), rep(rep(1:nc, each=nr), nc), rep(1:nc, each=nc*nr))
  u <- u[u[,2] <= u[,3],]  
  if (any(NVL(a$base,0)!=0)) { u <- u[-a$base,] }
  coef.names <- paste("b1twostar", a$b1attrname, b1u[u[,1]],  a$b2attrname,
                      apply(matrix(b2u[u[,2:3]],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b1twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]), minval = 0)
}

################################################################################
InitErgmTerm.b2concurrent<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("by"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  byarg <- a$by
  nb1 <- get.network.attribute(nw, "bipartite")
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "b2concurrent")[-seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b2concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_b2concurrent_by_attr function
    coef.names <- paste("b2concurrent",".", byarg,u, sep="")
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
InitErgmTerm.b2cov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("attrname","transform","transformname"),
                      vartypes = c("character","function","character"),
                      defaultvalues = list(NULL,function(x)x,""),
                      required = c(TRUE,FALSE,FALSE))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste(paste("nodeicov",f.name,sep=""),attrname,sep=".")
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- f(get.node.attr(nw, attrname, "nodeicov", numeric=TRUE)[(nb1+1):network.size(nw)])
  list(name="b2cov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}


################################################################################
InitErgmTerm.b2degrange<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE,
                      varnames = c("from", "to", "by", "homophily"),
                      vartypes = c("numeric", "numeric", "character", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term odegrange must have from<to.")

  nb1 <- get.network.attribute(nw, "bipartite")
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "b2degrange")
    if(!homophily) nodecov <- nodecov[-seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to b2degrange() has only one value", call.=FALSE)
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
                         paste("b2deg",from,"+", ".homophily.",byarg,sep=""),
                         paste("b2deg",from,"to",to, ".homophily.",byarg,sep=""))
    name <- "b2degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_b2degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("b2deg",du[1,],"+.", byarg, u[du[3,]],sep=""),
                         paste("b2deg",du[1,],"to",du[2,],".",byarg, u[du[3,]],sep=""))
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
InitErgmTerm.b2degree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "by"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if (!is.null(a$by)) {  # CASE 1:  a$by GIVEN
    nodecov <- get.node.attr(nw, a$by)[-seq_len(nb1)]
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
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
    coef.names <- paste("b2deg", du[1,], ".", a$by, u[du[2,]], sep="")
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
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats, minval=0, maxval=network.size(nw)-nb1)
}

################################################################################
InitErgmTerm.b2factor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  attrname<-a$attrname
  base <- a$base
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, attrname, "b2factor")[(nb1+1):network.size(nw)]
  
  if(all(is.na(nodecov)))
	  stop("Argument to b2factor() does not exist", call.=FALSE)
  
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to b2factor() has only one value", call.=FALSE)
  }
  if(base[1]==0){
    coef.names <- paste("b2factor", attrname, paste(u), sep=".")
    inputs <- c(ui, nodecov)
    attr(inputs, "ParamsBeforeCov") <- lu
  }else{
    coef.names <- paste("b2factor",attrname, paste(u[-base]), sep=".")
    inputs <- c(ui[-base], nodecov)
    attr(inputs, "ParamsBeforeCov") <- lu-length(base)
  }
  list(name="b2factor", coef.names=coef.names, inputs=inputs, dependence=FALSE, minval=0)
}



################################################################################
InitErgmTerm.b2star <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "istar"
    coef.names <- paste("b2star", a$k, ".", a$attrname, sep="")
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
InitErgmTerm.b2starmix <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname", "base", "diff"),
                       vartypes = c("numeric", "character", "numeric", "logical"),
                       defaultvalues = list(NULL, NULL, NULL, TRUE),
                       required = c(TRUE, TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, a$attrname)
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { stop("Only a single scalar k may be used with each b2starmix term") }
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
    coef.names <- paste("b2starmix", a$k, a$attrname,
                        apply(matrix(namescov[u[,2:1]],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- nr+(1:nc)
    if (any(NVL(a$base,0)!=0)) { u <- u[-a$base] }
    name <- "b2starmixhomophily"
    coef.names <- paste("b2starmix", a$k, a$attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs, minval=0)
}

################################################################################
InitErgmTerm.b2twostar <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("b1attrname", "b2attrname", "base"),
                       vartypes = c("character", "character", "numeric"),
                       defaultvalues = list(NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  b1nodecov <- get.node.attr(nw, a$b1attrname)[1:nb1]
  b1u<-sort(unique(b1nodecov))
  if(any(is.na(b1nodecov))){ b1u<-c(b1u,NA) }
  if(is.null(a$b2attrname)) { a$b2attrname <- a$b1attrname }
  b2nodecov <- get.node.attr(nw, a$b2attrname)[(1+nb1):n]
  b2u<-sort(unique(b2nodecov))
  if(any(is.na(b2nodecov))){b2u<-c(b2u,NA)}
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  nr <- length(b1u)
  nc <- length(b2u)
  u <- cbind(rep(1:nc, nr*nr), rep(rep(1:nr, each=nc), nr), rep(1:nr, each=nc*nr))
  u <- u[u[,2] <= u[,3],]  
  if (any(NVL(a$base,0)!=0)) { u <- u[-a$base,] }
  coef.names <- paste("b2twostar", a$b2attrname, b2u[u[,1]],  a$b1attrname,
                      apply(matrix(b1u[u[,2:3]],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b2twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]), minval=0)
}

################################################################################
InitErgmTerm.balance<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)){
   nodecov <- get.node.attr(nw, attrname, "balance")
   u<-sort(unique(nodecov))
   if(any(is.na(nodecov))){u<-c(u,NA)}
#
#  Recode to numeric if necessary
#
   nodecov <- match(nodecov,u,nomatch=length(u)+1)
   ui <- seq(along=u)

   if (length(u)==1)
         stop ("Attribute given to balance() has only one value", call.=FALSE)
   if (!diff) {
     #     No parameters before covariates here, so no need for "ParamsBeforeCov"
     coef.names <- paste("balance",attrname,sep=".")
     inputs <- c(nodecov)
   } else {
     #  Number of input parameters before covariates equals number of
     #  unique elements in nodecov, namely length(u)
     coef.names <- paste("balance",attrname, u, sep=".")
     inputs <- c(ui, nodecov)
     attr(inputs, "ParamsBeforeCov") <- length(ui)
   }
  }else{
    coef.names <- "balance"
    inputs <- NULL
  }
  list(name="balance", coef.names=coef.names, inputs=inputs, dependence=TRUE, minval=0)
}



#=======================InitErgmTerm functions:  C============================#


################################################################################
InitErgmTerm.concurrent<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("by"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  byarg <- a$by
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "concurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)

  }
  if(!is.null(byarg)) {
    if(length(u)==0) {return(NULL)}
    # See comment in d_concurrent_by_attr function
    coef.names <- paste("concurrent",".", byarg,u, sep="")
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
InitErgmTerm.ctriple<-InitErgmTerm.ctriad<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname","diff"),
                      vartypes = c("character","logical"),
                      defaultvalues = list(NULL,FALSE),
                      required = c(FALSE,FALSE))
  attrname <- a$attrname; diff <- a$diff;
  if(!is.null(attrname)){
    nodecov <- get.node.attr(nw, attrname, "ctriple")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ctriple() has only one value", call.=FALSE)
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
                     varnames = c("k"),
                     vartypes = c("numeric"),
                     defaultvalues = list(NULL),
                     required = c(TRUE))
  ### Process the arguments
  if (length(a$k)==0) return(NULL)
  ### Construct the list to return
  list(name="cycle",                            #name: required
       coef.names = paste("cycle", a$k, sep=""),  #coef.names: required
       inputs = c(max(a$k), (2:max(a$k)) %in% a$k),
       minval = 0)
}



#=======================InitErgmTerm functions:  D============================#


################################################################################
InitErgmTerm.degcor<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
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
InitErgmTerm.degrange<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("from", "to", "by", "homophily"),
                      vartypes = c("numeric", "numeric", "character", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degrange")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degrange() has only one value", call.=FALSE)
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
                         paste("deg",from,"+", ".homophily.",byarg,sep=""),
                         paste("deg",from,"to",to, ".homophily.",byarg,sep=""))
    name <- "degrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("deg",du[1,],"+.", byarg, u[du[3,]],sep=""),
                         paste("deg",du[1,],"to",du[2,],".",byarg, u[du[3,]],sep=""))
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
InitErgmTerm.degree<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("d", "by", "homophily"),
                      vartypes = c("numeric", "character", "logical"),
                      defaultvalues = list(NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree() has only one value", call.=FALSE)
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
    coef.names <- paste("deg", d, ".homophily.",byarg, sep="")
    name <- "degree_w_homophily"
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]], sep="")
    name <- "degree_by_attr"
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
InitErgmTerm.degreepopularity<-function (nw, arglist, ...) {
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
InitErgmTerm.diff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow", "dir", "sign.action"),
                      vartypes = c("character","numeric", "character", "character"),
                      defaultvalues = list(NULL,1, "t-h", "identity"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  DIRS <- c("t-h", "tail-head", "b1-b2",
            "h-t", "head-tail", "b2-b1")
  dir <- match.arg(tolower(a$dir), DIRS)
  dir.mul <- if(match(dir, DIRS)<=3) +1 else -1
  
  SIGN.ACTIONS <- c("identity", "abs", "posonly", "negonly")
  sign.action <- match.arg(tolower(a$sign.action), SIGN.ACTIONS)
  sign.code <- match(sign.action, SIGN.ACTIONS)

  if(sign.action!="abs" && !is.directed(nw)) message("Note that behavior of term diff() on undirected networks may be unexpected. See help(\"ergm-terms\") for more information.")
  
  # 1 and 4 are sign codes that allow negative differences.
  if(sign.code %in% c(1, 4) &&  a$pow!=round(a$pow)) stop("In term diff(attrname, pow, sign=",a$sign,"), pow must be an integer.")
  
  ### Construct the list to return
  list(name="diff",                                     #name: required
       coef.names = paste0("diff", if(a$pow!=1) a$pow else "", if(sign.action!="identity") paste0(".", sign.action), if(sign.action!="abs") paste0(".", dir), ".", a$attrname), #coef.names: required
       inputs = c(a$pow, dir.mul, sign.code, nodecov),  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}


################################################################################
InitErgmTerm.dsp<-function(nw, arglist, ...) {
# the following line was commented out in <InitErgm.dsp>:  
#   ergm.checkdirected("dsp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
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
         inputs=c(d), emptynwstats=emptynwstats, minval = 0)
  }else{
    list(name=dname, coef.names=paste("dsp",d,sep=""),inputs=c(d), minval = 0)
  }
}


################################################################################
InitErgmTerm.dyadcov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x","attrname"),
                      vartypes = c("matrix,network","character"),
                      defaultvalues = list(NULL,NULL),
                      required = c(TRUE,FALSE))
  x<-a$x;attrname<-a$attrname
  #Coerce x to an adjacency matrix
  if(is.network(x))
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
  else if(is.character(x))
#   xm<-as.matrix.network(nw,matrix.type="adjacency",x)
    xm<-get.network.attribute(nw,x)
  else
    xm<-as.matrix(x)


#Update the terms list, adding the vectorized adjacency matrix
  if(!is.null(attrname))
    cn<-paste("dyadcov", as.character(sys.call(0)[[3]][2]), 
              as.character(attrname), sep = ".")
  else
    cn<-paste("dyadcov", as.character(sys.call(0)[[3]][2]), sep = ".")
 
  if(is.directed(nw)){
   #Check for symmetry
   # DH:  Since nw is directed, why are we testing for symmetry here?  
   if (any(xm[upper.tri(xm)]!=t(xm)[upper.tri(xm)])){
     xm[lower.tri(xm)]<-t(xm)[lower.tri(xm)]
     warning("asymmetric covariate in dyadcov; using upper triangle only")
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
      stop("There is no network attribute named ",a$x,call.=FALSE)
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
InitErgmTerm.esp<-function(nw, arglist, ...) {
# the following line was commented out in <InitErgm.esp>:  
#    ergm.checkdirected("esp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}

  list(name=dname, coef.names=paste("esp",d,sep=""), inputs=c(d), minval=0)
}



#=======================InitErgmTerm functions:  G============================#

################################################################################
InitErgmTerm.gwb1degree<-function(nw, arglist, initialfit=FALSE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
  # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                      varnames = c("decay", "fixed", "attrname","cutoff"),
                      vartypes = c("numeric", "logical", "character","numeric"),
                      defaultvalues = list(0, TRUE, NULL, 30),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:(network.size(nw) - nb1)
  maxesp <- min(cutoff, network.size(nw)-nb1)

  d <- 1:maxesp
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwb1degree term is not yet able to handle a ",
           "non-fixed decay term.", call.=FALSE) # with an attribute.")
#    }
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^
                  {i-1}*(1+i-exp(-x[2])))
            )
    }
    list(name="b1degree", coef.names=paste("gwb1degree#",d,sep=""),
         inputs=c(d), params=list(gwb1degree=NULL,gwb1degree.decay=decay),
         map=map, gradient=gradient)
  } else {
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb1degree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb1degree() has only one value", call.=FALSE)
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
      coef.names <- paste("gwb1deg",decay,sep="")
      inputs <- c(decay)
    }
    list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
  }
}



################################################################################
InitErgmTerm.gwb2degree<-function(nw, arglist, initialfit=FALSE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
  # default for 'fixed' should be made 'FALSE' when the function can handle it!                    
                      varnames = c("decay", "fixed", "attrname","cutoff"),
                      vartypes = c("numeric", "logical", "character", "numeric"),
                      defaultvalues = list(0, TRUE, NULL, 30),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  cutoff<-a$cutoff
  nb1 <- get.network.attribute(nw,"bipartite")
# d <- 1:nb1
  maxesp <- min(cutoff,nb1)
  d <- 1:maxesp
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwb2degree term is not yet able to handle a ",
           "non-fixed decay term.", call.=FALSE) # with an attribute.")
#    }
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^
                  {i-1}*(1+i-exp(-x[2])))
            )
    }
    list(name="b2degree", coef.names=paste("gwb2degree#",d,sep=""),
         inputs=c(d), params=list(gwb2degree=NULL,gwb2degree.decay=decay),
         map=map, gradient=gradient)
  } else { 
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb2degree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb2degree() has only one value", call.=FALSE)
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
      coef.names <- paste("gwb2deg",decay,sep="")
      inputs <- c(decay)
    }
    list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
  }
}



################################################################################
InitErgmTerm.gwdegree<-function(nw, arglist, initialfit=FALSE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("decay", "fixed", "attrname","cutoff"),
                      vartypes = c("numeric", "logical", "character", "numeric"),
                      defaultvalues = list(0, FALSE, NULL, 30),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed && !initialfit) {
    warning("The gwdegree term cannot yet handle a nonfixed decay ",
            "term with an attribute.  Switching to fixed=TRUE.", call.=FALSE)
    fixed <- TRUE
  }
  if(!initialfit && !fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    list(name="degree", coef.names=paste("gwdegree#",d,sep=""), 
         inputs=c(d), params=list(gwdegree=NULL,gwdegree.decay=decay),
         map=map, gradient=gradient, conflicts.constraints="degreedist")
  } else {
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwdegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwdegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwdegree_by_attr"
      coef.names <- paste("gwdeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
      list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
    }else{
      name <- "gwdegree"
      coef.names <- "gwdegree"
      inputs <- c(decay)
      list(name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="degreedist")
    }
  }
}



################################################################################
InitErgmTerm.gwdsp<-function(nw, arglist, initialfit=FALSE, ...) {
# the following line was commented out in <InitErgm.gwdsp>:  
#   ergm.checkdirected("gwdsp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","alpha"),
                      vartypes = c("numeric","logical","numeric","numeric"),
                      defaultvalues = list(0, FALSE, 30, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  if(!initialfit && !fixed){ # This is a curved exponential family model
#   d <- 1:(network.size(nw)-1)
    maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    map<- function(x,n,...) {
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
    list(name=dname, coef.names=paste("gwdsp#",d,sep=""), 
         inputs=c(d), params=list(gwdsp=NULL,gwdsp.decay=decay),
         map=map, gradient=gradient)
  }else{
    if (initialfit && !fixed) # First pass to get MPLE coefficient
      coef.names <- "gwdsp"   # must match params$gwdsp above
    else  # fixed == TRUE
      coef.names <- paste("gwdsp.fixed.",decay,sep="")
  if(is.directed(nw)){dname <- "gwtdsp"}else{dname <- "gwdsp"}
  list(name=dname, coef.names=coef.names, inputs=c(decay))
  }
}



################################################################################
InitErgmTerm.gwesp<-function(nw, arglist, initialfit=FALSE, ...) {
# the following line was commented out in <InitErgm.gwesp>:
#   ergm.checkdirected("gwesp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(0, FALSE, 30, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!initialfit && !fixed){ # This is a curved exponential family model
#   d <- 1:(network.size(nw)-2)
     maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...){
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...){
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
    list(name=dname, coef.names=paste("esp#",d,sep=""), 
         inputs=c(d), params=list(gwesp=NULL,gwesp.decay=decay),
         map=map, gradient=gradient)
  }else{
    if (initialfit && !fixed)  # First pass to get MPLE coefficient
      coef.names <- "gwesp"
    else # fixed == TRUE
      coef.names <- paste("gwesp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtesp"}else{dname <- "gwesp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay))
  }
}



################################################################################
InitErgmTerm.gwidegree<-function(nw, arglist, initialfit=FALSE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff"),
                      vartypes = c("numeric", "logical", "character","numeric"),
                      defaultvalues = list(0, FALSE, NULL, 30),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
  maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed && !initialfit) {
    warning("The gwidegree term cannot yet handle a nonfixed decay ",
            "term with an attribute.  Switching to fixed=TRUE.", call.=FALSE)
    fixed <- TRUE
  }
  if(!initialfit && !fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    list(name="idegree", coef.names=paste("gwidegree#",d,sep=""), 
         inputs=c(d), params=list(gwidegree=NULL,gwidegree.decay=decay),
         map=map, gradient=gradient, conflicts.constraints="idegreedist")
  } else { 
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwidegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwidegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwidegree_by_attr"
      coef.names <- paste("gwideg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
      list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
    }else{
      name <- "gwidegree"
      coef.names <- "gwidegree"
      inputs <- c(decay)
      list(name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="idegreedist")
    }
  }
}


################################################################################
InitErgmTerm.gwnsp<-function(nw, arglist, initialfit=FALSE, ...) {
# the following line was commented out in <InitErgm.gwnsp>:
#    ergm.checkdirected("gwnsp", is.directed(nw), requirement=FALSE)
# so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff", "alpha"),
                      vartypes = c("numeric","logical","numeric", "numeric"),
                      defaultvalues = list(0, FALSE, 30, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }

  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  if(!initialfit && !fixed){ # This is a curved exponential family model
#   d <- 1:(network.size(nw)-1)
     maxesp <- min(cutoff,network.size(nw)-2)
    d <- 1:maxesp
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...){
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...){
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    if(is.directed(nw)){dname <- "tnsp"}else{dname <- "nsp"}
    list(name=dname, coef.names=paste("nsp#",d,sep=""),
         inputs=c(d), params=list(gwnsp=NULL,gwnsp.decay=decay),
         map=map, gradient=gradient)
  }else{
    if (initialfit && !fixed)  # First pass to get MPLE coefficient
      coef.names <- "gwnsp"
    else # fixed == TRUE
      coef.names <- paste("gwnsp.fixed.",decay,sep="")
    if(is.directed(nw)){dname <- "gwtnsp"}else{dname <- "gwnsp"}
    list(name=dname, coef.names=coef.names, inputs=c(decay))    
  }
}


################################################################################
InitErgmTerm.gwodegree<-function(nw, arglist, initialfit=FALSE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff"),
                      vartypes = c("numeric", "logical", "character","numeric"),
                      defaultvalues = list(0, FALSE, NULL, 30),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff
# d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed && !initialfit) {
    warning("The gwodegree term cannot yet handle a nonfixed decay ",
            "term with an attribute.  Switching to fixed=TRUE.", call.=FALSE)
    fixed <- TRUE
  }
  if(!initialfit && !fixed){ # This is a curved exponential family model
    ld<-length(d)
    if(ld==0){return(NULL)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    list(name="odegree", coef.names=paste("gwodegree#",d,sep=""),
         inputs=c(d), params=list(gwodegree=NULL,gwodegree.decay=decay),
         map=map, gradient=gradient, conflicts.constraints="odegreedist")
  } else {
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwodegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwodegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwodegree_by_attr"
      coef.names <- paste("gwodeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
      list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
    }else{
      name <- "gwodegree"
      coef.names <- "gwodegree"
      inputs <- c(decay)
      list(name=name, coef.names=coef.names, inputs=inputs, conflicts.constraints="odegreedist")
    }
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
  if(is.network(a$x)){													# Arg to hamming is a network
    # check for attribute existance before creating matrix
  
    if( is.null(a$attrname) || is.null(get.edge.attribute(a$x,a$attrname))){ 
      xm<-as.edgelist(a$x)  # so call the non attribute version
    } else {
      xm<-as.edgelist(a$x,a$attrname)
    }
    
    
  }else if(is.character(a$x)){												# Arg to hamming is the name of an attribute in nw
    xm<-get.network.attribute(nw,a$x)
    xm<-as.edgelist(xm)
  }else if(is.null(a$x)){
    xm<-as.edgelist(nw)								# Arg to hamming does not exist; uses nw
  }else if(is.matrix(a$x) && ncol(a$x)!=2){
    xm<-as.edgelist(network.update(network.copy(nw),a$x,matrix.type="adjacency"))
  }else{
    xm<-as.matrix(a$x)													# Arg to hamming is anything else; attempts to coerce
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
      stop("Improper dyadic covariate passed to hamming()", call.=FALSE)
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
    xm <- ergm.Cprepare.el(xm, prototype=nw)
  }
  if (!is.null(covm)) {
    covm <- ergm.Cprepare.el(covm, prototype=nw)
  }else covm <- 0
  inputs <- c(xm, a$defaultweight, covm)
  list(name="hamming", coef.names=coef.names, #name and coef.names: required 
       inputs = inputs, emptynwstats = emptynwstats, dependence = FALSE,
       minval = minval, maxval = maxval)
}

################################################################################
InitErgmTerm.hammingmix<-function (nw, arglist, ...) {
  # There is no reason hammingmix should be directed-only, but for now
  # the undirected version does not seem to work properly, so:
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname","x","base","contrast"),
                      vartypes = c("character","matrix,network","numeric","logical"),
                      defaultvalues = list(NULL,nw,0,FALSE),
                      required = c(TRUE,FALSE,FALSE,FALSE))
  attrname<-a$attrname
  x<-a$x
  base<-a$base
  if (a$contrast) {
    stop("The 'contrast' argument of the hammingmix term is deprecated.  Use 'base' instead")
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
    stop("hammingmix() requires an edgelist", call.=FALSE)
  }
    nodecov <- get.node.attr(nw, attrname, "hammingmix")
    mixmat <- mixingmatrix(nw,attrname)$mat
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)))
#   if(!is.directed(nw)){
#    u <- u[row(mixmat) >= col(mixmat)]
#   }
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
#
#   Recode to numeric if necessary
#
    namescov <- sort(unique(nodecov))
    nodecov <- match(nodecov,namescov)
    if (length(nodecov)==1)
        stop ("Argument to hammingmix() has only one value", call.=FALSE)

  if (any(NVL(base,0)!=0)) {
    u <- u[-base,]
  }
  coef.names <- paste("hammingmix",attrname,
                      apply(matrix(namescov[u],ncol=2),1,paste,collapse="."), 
                      sep=".")
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui),
  inputs=c(ergm.Cprepare.el(xm, prototype=nw), u[,1], u[,2],nodecov)
  attr(inputs, "ParamsBeforeCov") <- nrow(u)
  # The emptynwstats code below does not work right for
  # undirected networks, mostly since hammingmix doesn't work 
  # in this case anyway.
  nw %v% "_tmp_nodecov" <- nodecov
  emptynwstats <- summary(nw ~ nodemix("_tmp_nodecov", base))
  list(name="hammingmix", coef.names=coef.names, inputs=inputs, 
       emptynwstats=emptynwstats, dependence=FALSE)
}




#=======================InitErgmTerm functions:  I============================#

################################################################################
InitErgmTerm.idegrange<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("from", "to", "by", "homophily"),
                      vartypes = c("numeric", "numeric", "character", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term idegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term idegrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "idegrange")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to idegrange() has only one value", call.=FALSE)
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
                         paste("ideg",from,"+", ".homophily.",byarg,sep=""),
                         paste("ideg",from,"to",to, ".homophily.",byarg,sep=""))
    name <- "idegrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_idegrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("ideg",du[1,],"+.", byarg, u[du[3,]],sep=""),
                         paste("ideg",du[1,],"to",du[2,],".",byarg, u[du[3,]],sep=""))
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
InitErgmTerm.idegree<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "homophily"),
                      vartypes = c("numeric", "character", "logical"),
                      defaultvalues = list(NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "idegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to idegree() has only one value", call.=FALSE)
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
    coef.names <- paste("ideg", d, ".homophily.",byarg, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "idegree_by_attr"
    # See comment in d_idegree_by_attr function
    coef.names <- paste("ideg", du[1,], ".", byarg,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name, coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE)
  }else{
    list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval = 0, maxval=network.size(nw), conflicts.constraints="idegreedist")
  }
}


################################################################################
InitErgmTerm.idegreepopularity<-function (nw, arglist, ...) {
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
InitErgmTerm.istar<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("k", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  k <- a$k
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "istar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
#     Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to istar() has only one value", call.=FALSE)
  }else{
  }
  lk<-length(k)
  if(lk==0){return(NULL)}
  if(!is.null(attrname)){
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
InitErgmTerm.kstar<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("k", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  k<-a$k;attrname<-a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "kstar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
#    Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to kstar() has only one value", call.=FALSE)
  }
  lk<-length(k)
  if(lk==0){return(NULL)}
  if(!is.null(attrname)){
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
    warning("localtriangle requires an undirected neighborhood. Using only mutual ties.", call.=FALSE)
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
InitErgmTerm.mutual<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("same", "by", "diff", "keep"),
                      vartypes = c("character", "character", "logical", "numeric"),
                      defaultvalues = list(NULL, NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE))
  ### Process the arguments
  if (!is.null(a$same) || !is.null(a$by)) {
    if (!is.null(a$same)) {
     attrname <- a$same
     if (!is.null(a$by)) 
       warning("Ignoring 'by' argument to mutual because 'same' exists", call.=FALSE)
    }else{
     attrname <- a$by
    }
    nodecov <- get.node.attr(nw, attrname)
    u <- sort(unique(nodecov))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
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
  list(name=name,                      #name: required
       coef.names = coef.names,        #coef.names: required
       inputs=inputs,
       minval = 0,
       maxval = network.dyadcount(nw,FALSE)/2) 
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
InitErgmTerm.nodecov<-InitErgmTerm.nodemain<-function (nw, arglist, ...) {
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
  list(name="nodecov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}



################################################################################
InitErgmTerm.nodefactor<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments

  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      return()
    }
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodefactor",                                        #required
       coef.names = paste("nodefactor", paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}  

################################################################################
InitErgmTerm.nodeicov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname","transform","transformname"),
                      vartypes = c("character","function","character"),
                      defaultvalues = list(NULL,function(x)x,""),
                      required = c(TRUE,FALSE,FALSE))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste(paste("nodeicov",f.name,sep=""),attrname,sep=".")
  nodecov <- f(get.node.attr(nw, attrname, "nodeicov", numeric=TRUE))
  list(name="nodeicov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}



################################################################################
InitErgmTerm.nodeifactor<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = c("attrname", "base"),          
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments

  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      return()
    }
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeifactor",                                        #required
       coef.names = paste("nodeifactor", paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################
InitErgmTerm.nodematch<-InitErgmTerm.match<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "diff", "keep"),
                      vartypes = c("character", "logical", "numeric"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
  u <- sort(unique(nodecov))
  if (!is.null(a$keep)) {
    u <- u[a$keep]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov==(length(u)+1)
  nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
  ui <- seq(along=u)
  ### Construct the list to return
  if (a$diff) {
    coef.names <- paste("nodematch", paste(a$attrname,collapse="."), u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", paste(a$attrname,collapse="."), sep=".")
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
InitErgmTerm.nodemix<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  if (is.bipartite(nw) && is.directed(nw)) {
    stop("Directed bipartite networks are not currently possible")
  }
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
    
  if (is.bipartite(nw)) {
    #  So undirected network storage but directed mixing
    nb1 <- get.network.attribute(nw, "bipartite")       
    #  Recode nodecov to numeric (but retain original sorted names in "namescov")
    b1namescov <- sort(unique(nodecov[1:nb1]))
    b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
    namescov <- c(b1namescov, b2namescov)
    b1nodecov <- match(nodecov[1:nb1],b1namescov)
    b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    nodecov <- c(b1nodecov, b2nodecov + nr)
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}    
    if (any(NVL(a$base,0)!=0)) {
      u <- u[-a$base,]
    }
    name <- "mix"
    cn <- paste("mix", paste(a$attrname,collapse="."), apply(matrix(namescov[u],ncol=2),
                                       1,paste,collapse="."), sep=".")
    inputs <- c(u[,1], u[,2], nodecov)
    attr(inputs, "ParamsBeforeCov") <- NROW(u)
  } else {# So one mode, but could be directed or undirected
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    #   Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    ucount<-sapply(ui,function(x){sum(nodecov==x,na.rm=TRUE)}) #Count cases
    uui <- matrix(1:length(ui)^2,length(ui),length(ui))  #Create int tables
    urm <- t(sapply(ui,rep,length(ui)))   #This is the reverse of what you'd
    ucm <- sapply(ui,rep,length(ui))      #expect for r/c, but it's correct
    uun <- outer(u,u,paste,sep=".")
    if (!is.directed(nw)) {
      uui <- uui[upper.tri(uui,diag=TRUE)]
      urm <- urm[upper.tri(urm,diag=TRUE)]  
      ucm <- ucm[upper.tri(ucm,diag=TRUE)]
      uun <- uun[upper.tri(uun,diag=TRUE)]
    }
    if (any(NVL(a$base,0)!=0)) {
      urm <- as.vector(urm)[-a$base]
      ucm <- as.vector(ucm)[-a$base]
      uun <- as.vector(uun)[-a$base]
    }
    name <- "nodemix"
    cn <- paste("mix", paste(a$attrname,collapse="."), uun, sep=".")
    inputs <- c(urm, ucm, nodecov)
    #attr(inputs, "ParamsBeforeCov") <- 2*length(uui)
    attr(inputs, "ParamsBeforeCov") <- 2*length(uun)
  }
  ### Construct the list to return
  list(name = name, coef.names = cn, # required
       inputs = inputs, 
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}

################################################################################
InitErgmTerm.nodeocov<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = c("attrname","transform","transformname"),
                      vartypes = c("character","function","character"),
                      defaultvalues = list(NULL,function(x)x,""),
                      required = c(TRUE,FALSE,FALSE))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste(paste("nodeocov",f.name,sep=""),attrname,sep=".")
  nodecov <- f(get.node.attr(nw, attrname, "nodeocov", numeric=TRUE))
  list(name="nodeocov", coef.names=coef.names, inputs=c(nodecov), dependence=FALSE)
}



################################################################################
InitErgmTerm.nodeofactor<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (any(NVL(a$base,0)!=0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      return()
    }
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)

  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeofactor",                                        #required
       coef.names = paste("nodeofactor", paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE, # So we don't use MCMC if not necessary
       minval = 0
       )
}  

################################################################################
InitErgmTerm.nsp<-function(nw, arglist, ...) {
# The following line was commented out in <InitErgm.nsp>
#   ergm.checkdirected("nsp", is.directed(nw), requirement=FALSE)
# so I have not included 'directed=TRUE' in the call to <check.ErgmTerm>
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
         emptynwstats=emptynwstats, minval=0)
  } else {
    list(name=dname, coef.names=coef.names, inputs=c(d), minval=0)
  }
}



#=======================InitErgmTerm functions:  O============================#

################################################################################
InitErgmTerm.odegrange<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("from", "to", "by", "homophily"),
                      vartypes = c("numeric", "numeric", "character", "logical"),
                      defaultvalues = list(NULL, Inf, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  from<-a$from; to<-a$to; byarg <- a$by; homophily <- a$homophily
  to <- ifelse(to==Inf, network.size(nw)+1, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term odegrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term odegrange must have from<to.")

  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "odegrange")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to odegrange() has only one value", call.=FALSE)
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
                         paste("odeg",from,"+", ".homophily.",byarg,sep=""),
                         paste("odeg",from,"to",to, ".homophily.",byarg,sep=""))
    name <- "odegrange_w_homophily"
    inputs <- c(rbind(from,to), nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_odegrange_by_attr function
    coef.names <- ifelse(du[2,]>=network.size(nw)+1,
                         paste("odeg",du[1,],"+.", byarg, u[du[3,]],sep=""),
                         paste("odeg",du[1,],"to",du[2,],".",byarg, u[du[3,]],sep=""))
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
InitErgmTerm.odegree<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "homophily"),
                      vartypes = c("numeric", "character", "logical"),
                      defaultvalues = list(NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "odegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to odegree() has only one value", call.=FALSE)
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
    coef.names <- paste("odeg", d, ".homophily.",byarg, sep="")
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    name <- "odegree_by_attr"
    # See comment in d_odegree_by_attr function
    coef.names <- paste("odeg", du[1,], ".", byarg,u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  }
  if (!is.null(emptynwstats)){
    list(name=name, coef.names=coef.names, inputs=inputs,
         emptynwstats=emptynwstats, dependence=TRUE, minval=0)
  }else{
    list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE, minval=0, maxval=network.size(nw), conflicts.constraints="odegreedist")
  }
}


################################################################################
InitErgmTerm.odegreepopularity<-function (nw, arglist, ...) {
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
InitErgmTerm.ostar<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("k", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  k<-a$k
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "ostar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to ostar() has only one value", call.=FALSE)
  }
  lk<-length(k)
  if(lk==0){return(NULL)}

  if(!is.null(attrname)){
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
InitErgmTerm.receiver<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("base"),
                      vartypes = c("numeric"),
                      defaultvalues = list(1),
                      required = c(FALSE))
  d <- 1:network.size(nw)
  if (any(NVL(a$base,0)!=0)) {
    d <- d[-a$base]
  }
  ld<-length(d)
  if(ld==0){return(NULL)}
  list(name="receiver", coef.names=paste("receiver",d,sep=""),
       inputs=c(d), emptynwstats=rep(0,length(d)), dependence=FALSE, minval=0, maxval=network.size(nw)-1, conflicts.constraints="idegrees")
}


#=======================InitErgmTerm functions:  S============================#

################################################################################
InitErgmTerm.sender<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("base"),
                      vartypes = c("numeric"),
                      defaultvalues = list(1),
                      required = c(FALSE))
  d <- 1:network.size(nw)
  if (any(NVL(a$base,0)!=0)) {
    d <- d[-a$base]
  }
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
InitErgmTerm.smalldiff<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "cutoff"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  cutoff <- a$cutoff
  attrname <- a$attrname
  if (length(cutoff)>1)
    stop("cutoff for smalldiff() must be a scalar.", call.=FALSE)
  coef.names <- paste("smalldiff.", attrname, cutoff, sep="")
  nodecov <- get.node.attr(nw, attrname, "smalldiff", numeric=TRUE)
  inputs <- c(cutoff, nodecov)
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="smalldiff", coef.names=coef.names, inputs=inputs,
       dependence=FALSE)
}

  

################################################################################
InitErgmTerm.sociality<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(FALSE, FALSE))
  attrname<-a$attrname
  d <- 1:network.size(nw)
  if (any(NVL(a$base,0)!=0)) {
    d <- d[-a$base]
  }
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "sociality")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to sociality() has only one value", call.=FALSE)
  }
  ld<-length(d)
  if(ld==0){return(NULL)}
  if(!is.null(attrname)){
    coef.names <- paste("sociality",d,".",attrname,sep="")
    inputs <- c(d, 0, nodecov) # Input requires a "guard" value.
  }else{
    coef.names <- paste("sociality",d,sep="")
    inputs <- c(d,0) # Input requires a "guard" value.
  }
  list(name="sociality", coef.names=coef.names, inputs=inputs, minval=0, maxval=network.size(nw)-1, conflicts.constraints="degrees")
}




#=======================InitErgmTerm functions:  T============================#


################################################################################
InitErgmTerm.threepath <- function(nw, arglist, ...) {
  warning("This term is inaccurately named and actually refers to a '3-trail' in that it counts repeated vertices: i-j-k-i is a 3-trail but not a 3-path. See ergm-terms help for more information. This name has been deprecated and will be removed in a future version: if a 3-trail is what you want, use the term 'threetrail'.")
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, 
                       varnames = c("keep"),
                       vartypes = c("numeric"),
                       defaultvalues = list(1:4),
                       required = c(FALSE))
  types <- c("RRR","RRL","LRR","LRL")[a$keep]
  if (is.directed(nw)) {
    return(list(name = "threetrail", 
                coef.names = paste("threetrail", types, sep="."),
                inputs=a$keep, minval = 0))
  }
  else {
    return(list(name = "threetrail", coef.names = "threetrail", minval = 0))
  }
}

################################################################################
InitErgmTerm.threetrail <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, 
                       varnames = c("keep"),
                       vartypes = c("numeric"),
                       defaultvalues = list(1:4),
                       required = c(FALSE))
  types <- c("RRR","RRL","LRR","LRL")[a$keep]
  if (is.directed(nw)) {
    return(list(name = "threetrail", 
                coef.names = paste("threetrail", types, sep="."),
                inputs=a$keep, minval = 0))
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
InitErgmTerm.triadcensus<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  d<-a$d
  emptynwstats<-NULL

  if(is.directed(nw)){
   tcn <- c("003","012", "102", "021D", "021U", "021C", "111D",
            "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
   if(is.null(d)){d <- 1:15}
   if(is.character(d)){d <- match(d, tcn)-1}
  }else{
#  Undirected
   tcn <- c("0", "1", "2", "3")
   if(is.null(d)){d <- 1:3}
   if(is.character(d)){d <- match(d, tcn)-1}
  }
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
InitErgmTerm.triangle<-InitErgmTerm.triangles<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "triangle")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to triangle() has only one value", call.=FALSE)
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
InitErgmTerm.tripercent<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "tripercent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to tripercent() has only one value", call.=FALSE)
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
InitErgmTerm.ttriple<-InitErgmTerm.ttriad<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("attrname", "diff"),
                      vartypes = c("character", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "ttriple")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ttriple() has only one value", call.=FALSE)
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


