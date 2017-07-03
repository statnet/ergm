#  File R/InitErgmTerm.dgw_sp.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#  ------------------------------------------------------------------ 
#   Description of the input and output parameters of the  
#   InitErgmTerm.xxx function, where xxx is the name of your term
#  ------------------------------------------------------------------ 
#
#  INPUTS:
#  Each InitErgmTerm function takes three arguments:
#	  		nw: The network of interest
#      arglist: The list of arguments passed to the term xxx
#         ... : There may be other arguments passed by 
#               ergm.getmodel, so each InitErgmTerm function 
#               must include the ... argument
#  These inputs are automatically supplied by ergm.getmodel.
#
#  OUTPUTS:
#  Each InitErgmTerm function should return a list.  
#     REQUIRED LIST ITEMS:
#          name: This names the C changestats function for term xxx, 
#                but does so by excluding the d_ prefix. The 
#                changestats function is named d_xxxy and 'name' is
#                consequently "xxxy". For example, the b1starmix
#                term has 2 changestats functions based on
#                whether the homophily argument is set. These are
#                d_b1starmix and d_b1starmixhomophily. The 'name' 
#                returned by InitErgmTerm.b1starmix is then one of 
#                "b1starmix" or "b1starmixhomophily" as appropriate.
#    coef.names: Vector of names for the coefficients (parameters)
#                as they will be reported in the output.
#       pkgname: This names the package containing the C changestats
#                function d_[name]. The default is "ergm", which means
#                that if you have code that exists as part of the 
#                (say) "ergm.userterms" package, you MUST specify 
#                pkgname="ergm.userterms"
#
#    OPTIONAL LIST ITEMS:
#        inputs: Vector of (double-precision numeric) inputs that the 
#                changestat function called d_'name' may require.
#                The default is NULL; no inputs are required.  But it
#                MUST be a vector!  Thus, if some of the inputs are,  
#                say, matrices, they must be "flattened" to vectors; if 
#                some are categorical character-valued variables, they
#                must be converted to numbers. Optionally, the inputs 
#                vector may have an attribute named "ParamsBeforeCov",
#                which is the number of input parameters preceding the 
#                covariate vector in 'inputs'.  This is necessary for 
#                compatibility with some of the existing d_xxx changestats 
#                functions in ergm, but is not necessary in general.
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  The default value is TRUE.
#  emptynwstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  For example, the degree0 term 
#                would require 'emptynwstats' since degree0 = number of 
#                nodes for the empty network.
#        params: For curved exponential family model terms only, a list of 
#                (numeric) initial values for the parameters of  
#                curved exponential family model terms. Each item in the  
#                list should be named with the corresponding parameter name 
#                (one or more of these will probably coincide with the 
#                 coef.names).  For example, the gwesp term returns 
#                params=list(gwesp=NULL,gwesp.decay=decay), where decay
#                was specified as an argument to the gwesp term. 
#           map: For curved exponential family model terms only, a function 
#                giving the map from the canonical parameters, theta,
#                associated with the statistics for this term, to eta, 
#                the corresponding curved parameters.  The length of eta 
#                is the same as the length of the 'params' list above.
#                The function takes two arguments:  theta and length(eta).
#      gradient: For curved exponential family model terms only, a function 
#                giving the gradient of the 'map'. If theta has length p 
#                and eta has length q, then gradient should return a
#                p by q matrix. This function takes two arguments:  theta 
#                and length(eta).
#


#  ------------------------------------------------------------------------- 
#   Description of the input parameters to the d_xxxy changestats function, 
#   where xxxy corresponds to the 'name' returned by InitErgmTerm.xxx.
#  -------------------------------------------------------------------------- 
#
#  INPUTS:
#  Each d_xxxy function takes five arguments:
#	    ntoggles: the number of toggles as described in 
#                 "ergm.userterms: A template package"
#          heads: a pointer to the array of the head nodes of the 
#                 proposed edges to be toggled
#          tails: a pointer to the array of the tail nodes of the
#                 proposed edges to be toggled
#            mtp: a pointer to the model, which includes the following:
#                 dstats      : a pointer to the array of changestats,
#                               macro-ed as CHANGE_STAT
#                 nstats      : the length of 'dstats', macro-ed as
#                               N_CHANGE_STATS
#                 inputparams : a pointer to the vector of input 
#                               parameters. This is supplied by the
#                               'inputs' returned by InitErgmTerm.xxx
#                               and is macro-ed as INPUT_PARAM
#                 ninputparams: the length of 'inputparams', macro-ed
#                               as N_INPUT_PARAMS
#            nwp: a pointer to the network.  This includes several 
#                 components and several macros exist for accessing
#                 these. See the changestat.h file for a list of these
#                 components and their macros. 
#  These inputs are automatically supplied to the d_xxxy function by the 
#  network_stats_wrapper function 

.sp.handle_layers <- function(nw, a, type, has_base){
  out <- list()

  if(is(a$Ls.path,"formula")) a$Ls.path <- list(a$Ls.path)
  L.path1 <- a$Ls.path[[1]]
  L.path2 <- a$Ls.path[[2]]
  L.base <- a$L.base

  if(is.null(L.path1) && is.null(L.path2) && is.null(L.base)) return(out)

  if(type=="RTP") stop("Layer-aware shared partner terms do not support reciprocated two-paths at this time.",call.=FALSE)

  
  NVL(L.path1) <- NVL(L.path2, L.base)
  NVL(L.path2) <- NVL(L.path1, L.base)
  if(has_base) NVL(L.base) <- NVL(L.path1, L.path2)
  
  nl <- max(.peek_vattrv(nw, ".LayerID"))

  layer0 <-
    if(has_base)
      as.formula(call("~",call("|",call("|", L.path1[[2]], L.path2[[2]]),L.base[[2]])))
    else
      as.formula(call("~",call("|", L.path1[[2]], L.path2[[2]])))

  out$auxiliaries <- .mk_.layer.net_auxform(layer0, nl)
  aux1 <- .mk_.layer.net_auxform(L.path1, nl)
  out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux1[[2]])
  aux2 <- .mk_.layer.net_auxform(L.path2, nl)
  out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux2[[2]])
  if(has_base){
    aux3 <- .mk_.layer.net_auxform(L.base, nl)
    out$auxiliaries[[2]] <- call("+", out$auxiliaries[[2]], aux3[[2]])
  }
  
  out$any_order <- if(type=="UTP" || (type%in%c("OSP","ISP") && !has_base)) TRUE else !a$L.in_order
  out$coef.names_prefix <- paste0(.lspec_coef.names(list(pth=c(NVL(L.path1,if(L.path2!=L.path1)L.path2)),bse=if(has_base) L.base,inord=a$L.in_order)),":")
  out$name_suffix <- "_ML"
  out$nw1 <- .split_constr_network(nw, ".LayerID", ".LayerName")[[1]] # Needed for emptynwstats.

  out
}




################################################################################
#Term to count ESP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per esp term.  The default, OTP, retains
#the original behavior of esp/gwesp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#
InitErgmTerm.desp<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for esp; valid types are:",paste(type.vec, collapse=","))
  dname<-"esp"
  if(is.directed(nw)){
    conam <- paste("esp",type,sep=".")
    typecode<-which(type==type.vec)
    dname <- "desp"
  }else{
    message("Use the ergm term 'esp' for undirected networks.")
    dname <- "desp"
    conam<-"esp"
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, TRUE)
  
  list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,paste(conam,d,sep="")), auxiliaries=linfo$auxiliaries, inputs=c(linfo$any_order,typecode,d), minval=0)
}


################################################################################
#Geometrically weighted edgewise shared partner term, where shared partners
#are defined by one of various specified types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per esp term.  The default, OTP, retains
#the original behavior of esp/gwesp.  In the undirected case, UTP is
#always used (since it is directedness-safe), and the user's input is
#overridden.  UTP cannot be chosen otherwise, since it won't work.
#
InitErgmTerm.dgwesp<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, 30,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for gwesp; valid types are:",paste(type.vec, collapse=","))
  dname<-"desp"
  if(!is.directed(nw)){  
    message("Use the gwesp term for undirected networks.")
    type <- "UTP"
    typecode<-0
    basenam<-paste("gwesp",sep=".")
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwesp",type,sep=".")
  }
  
  linfo <- .sp.handle_layers(nw, a, type, TRUE)
  
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwesp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

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
    params<-list(gwesp=NULL,gwesp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))

    list(name=paste0(dname,linfo$name_suffix),
         coef.names=paste0(linfo$coef.names_prefix,if(is.directed(nw)) paste("esp.",type,"#",d,sep="") else paste("esp#",d,sep="")),auxiliaries=linfo$auxiliaries, 
         inputs=c(linfo$any_order,typecode,d), params=params, map=map, gradient=gradient)
  }else{
    if(is.null(a$decay)) stop("Term 'dgwesp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwesp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if(is.directed(nw))
      coef.names <- paste(paste("gwesp",type,"fixed.",sep="."),decay, sep="")
    else
      coef.names <- paste("gwesp.fixed.",decay,sep="")

    list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,coef.names), inputs=c(linfo$any_order,decay,typecode,maxesp),auxiliaries=linfo$auxiliaries)
  }
}



#Term to count DSP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per dsp term.  The default, OTP, retains
#the original behavior of dsp/gwdsp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#
InitErgmTerm.ddsp<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for sp; valid types are:",paste(type.vec, collapse=","))
  dname<-"ddsp"
  if(is.directed(nw)){
    conam <- paste("dsp",type,sep=".")
    typecode<-which(type==type.vec)
    dname <- "ddsp"
  }else{
    message("Use the ergm term 'dsp' for undirected networks.")
    conam <- paste("dsp",sep=".")
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, FALSE)
  nw <- linfo$nw1
  
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
  
  list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,paste0(conam,d)), auxiliaries=linfo$auxiliaries, inputs=c(linfo$any_order,typecode,d), minval=0, emptynwstats=emptynwstats)
}



################################################################################
InitErgmTerm.dgwdsp<-function(nw, arglist, ...) {
  # the following line was commented out in <InitErgm.gwdsp>:
  #    ergm.checkdirected("gwdsp", is.directed(nw), requirement=FALSE)
  # so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, 30,"OTP", NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code; valid types are:",paste(type.vec, collapse=","))
  dname<-"ddsp"
  
  if(!is.directed(nw)){  
    message("Use the gwdsp term for undirected networks.")
    type <- "UTP"
    basenam<-"gwdsp"
    typecode<-0
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwdsp",type,sep=".")
  }

  linfo <- .sp.handle_layers(nw, a, type, FALSE)
  
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwdsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

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
    
    params<-list(gwdsp=NULL,gwdsp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))
    
    list(name=paste0(dname,linfo$name_suffix),
         coef.names=paste0(linfo$coef.names_prefix,if(is.directed(nw)) paste("dsp.",type,"#",d,sep="") else paste("dsp#",d,sep="")), 
         inputs=c(linfo$any_order,typecode,d), params=params,
         map=map, gradient=gradient, auxiliaries = linfo$auxiliaries)
  }else{
    if(is.null(a$decay)) stop("Term 'dgwdsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwdsp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if (is.directed(nw)) 
      coef.names <- paste("gwdsp",type,"fixed",decay,sep=".")
    else
      coef.names <- paste("gwdsp.fixed",decay,sep=".")
    
    list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,coef.names), inputs=c(linfo$any_order, decay,typecode,maxesp), auxiliaries=linfo$auxiliaries)
  }
}


#Term to count NSP statistics, where the shared partners may be any of
#several distinct types.
#
#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
#Only one type may be specified per NSP term.  The default, OTP, retains
#the original behavior of nsp/gwnsp.  In the undirected case, the UTP
#routine is used (since it is safe for undirected graphs), irrespective of
#the user's selection.  UTP cannot be chosen otherwise, since it won't work.
#
InitErgmTerm.dnsp<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","character","formula","formula,list","logical"),
                      defaultvalues = list(NULL,"OTP",NULL,NULL,FALSE),
                      required = c(TRUE, FALSE,FALSE,FALSE,FALSE))
  d<-a$d
  ld<-length(d)
  if(ld==0){return(NULL)}
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code for sp; valid types are:",paste(type.vec, collapse=","))
  dname<-"dnsp"
  if(is.directed(nw)){
    conam <- paste("nsp",type,sep=".")
    typecode<-which(type==type.vec)
  }else{
    message("Use the ergm term 'nsp' for undirected networks.")
    conam<-"nsp"
    type<-"UTP"
    typecode<-0
  }

  linfo <- .sp.handle_layers(nw, a, type, TRUE)
  nw <- linfo$nw1

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
  list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,paste0(conam,d)), auxiliaries=linfo$auxiliaries, inputs=c(linfo$any_order,typecode,d), minval=0, emptynwstats=emptynwstats)
}


################################################################################
InitErgmTerm.dgwnsp<-function(nw, arglist, ...) {
  # the following line was commented out in <InitErgm.gwnsp>:
  #    ergm.checkdirected("gwnsp", is.directed(nw), requirement=FALSE)
  # so, I've not passed 'directed=FALSE' to <check.ErgmTerm>  
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha","L.base","Ls.path","L.in_order"),
                      vartypes = c("numeric","logical","numeric","character", "numeric","formula","formula,list","logical"),
                      defaultvalues = list(NULL, FALSE, 30,"OTP", NULL,NULL,NULL,FALSE),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE,FALSE,FALSE,FALSE))
  if(!is.null(a$alpha)){
    stop("For consistency with gw*degree terms, in all gw*sp and dgw*sp terms the argument ", sQuote("alpha"), " has been renamed to " ,sQuote("decay"), ".", call.=FALSE)
  }
  
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...
  
  type<-toupper(a$type[1])
  type.vec<-c("OTP","ITP","RTP","OSP","ISP")
  if(!(type%in%type.vec))
    stop("Illegal type code; valid types are:",paste(type.vec, collapse=","))
  dname<-"dnsp"
  
  if(!is.directed(nw)){  
    message("Use the gwnsp term for undirected networks.")
    type <- "UTP"
    basenam<-"gwdsp"
    typecode<-0
  }else{
    typecode<-which(type==type.vec)
    basenam<-paste("gwnsp",type,sep=".")
  }
  
  linfo <- .sp.handle_layers(nw, a, type, TRUE)

  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'dgwnsp': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)

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
    
    params<-list(gwnsp=NULL,gwnsp.decay=decay)
    names(params)<-c(basenam,paste(basenam,"decay",sep="."))
    
    list(name=paste0(dname,linfo$name_suffix),
         coef.names=paste0(linfo$coef.names_prefix,if(is.directed(nw)) paste("nsp.",type,"#",d,sep="") else paste("nsp#",d,sep="")), 
         inputs=c(linfo$any_order,typecode,d), params=params,
         map=map, gradient=gradient, auxiliaries = linfo$auxiliaries)
  }else{
    if(is.null(a$decay)) stop("Term 'dgwnsp' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    dname<-"dgwnsp"
    maxesp <- min(cutoff,network.size(nw)-2)
    if (is.directed(nw)) 
      coef.names <- paste("gwnsp",type,"fixed",decay,sep=".")
    else
      coef.names <- paste("gwnsp.fixed",decay,sep=".")
    
    list(name=paste0(dname,linfo$name_suffix), coef.names=paste0(linfo$coef.names_prefix,coef.names), inputs=c(linfo$any_order, decay,typecode,maxesp), auxiliaries=linfo$auxiliaries)
  }
}
