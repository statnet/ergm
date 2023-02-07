#  File R/InitErgmTerm.dgw_sp.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

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
#               ergm_model, so each InitErgmTerm function 
#               must include the ... argument
#  These inputs are automatically supplied by ergm_model.
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

################################################################################

#Type codes are as follows (where (i,j) is the focal edge):
#
#  UTP - Undirected two-path (undirected graphs only)
#  OTP - Outgoing two-path (i->k->j)
#  ITP - Incoming two-path (i<-k<-j)
#  RTP - Reciprocated two-path (i<->k<->j)
#  OSP - Outgoing shared partner (i->k<-j)
#  ISP - Incoming shared partner (i<-k->j)
#
# The default, OTP, retains the original behavior of dsp/gwdsp.  If
# and only if the network is undirected, the UTP routine is used
# (since it is safe for undirected graphs), irrespective of the user's
# selection.
#

SPTYPE_DESC <- c(UTP = "undirected two-path",
                 OTP = "outgoing two-path",
                 ITP = "incoming two-path",
                 RTP = "reciprocated two-path",
                 OSP = "outgoing shared partner",
                 ISP = "incoming shared partner")

SPTYPE_CODE <- c(UTP = 0L, OTP = 1L, ITP = 2L, RTP = 3L, OSP = 4L, ISP = 5L)

.d_sp_resolve_type <- function(nw, a, undname, bip){
  type <- toupper(a$type[1])

  ## Deal with cases for which network type determines term SP type:
  if(is.bipartite(nw) && nchar(bip)) type <- c(b1="OSP", b2="ISP")[bip]
  else if(!is.directed(nw)) type <- "UTP"

  ## Check that the ultimate result is valid:
  if(! type%in%names(SPTYPE_CODE))
    ergm_Init_abort("Illegal value for ", sQuote("type"), "; valid types are:", paste.and(sQuote(names(SPTYPE_CODE))))

  type
}

.d_sp_impl <- function(sp, nw, arglist, cache.sp, emptynwstats = function(...)NULL, ...){
  if(sp %in% c("b1","b2")){
    bip <- sp
    sp <- "d"
  }else bip <- ""

  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("d","type"),
                      vartypes = c("numeric","character"),
                      defaultvalues = list(NULL,"OTP"),
                      required = c(TRUE, FALSE))
  utermname <- paste0(bip, sp, "sp")

  d<-a$d
  ld<-length(d)
  if(ld==0) return(NULL)

  type <- .d_sp_resolve_type(nw, a, utermname, bip)
  conam <- if(type=="UTP" || nchar(bip)) utermname else paste0(utermname,".",type)

  list(name=if(nchar(bip)) "ddspbwrap" else paste0("d",utermname), coef.names=paste(conam,d,sep=""), iinputs=c(SPTYPE_CODE[type],d),
       minval=0, emptynwstats=emptynwstats(nw=nw, a=a, d=d, type=type),
       auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
}

.dgw_sp_impl <- function(sp, nw, arglist, cache.sp, gw.cutoff, ...){
  if(sp %in% c("b1","b2")){
    bip <- sp
    sp <- "d"
  }else bip <- ""

  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("decay","fixed","cutoff","type", "alpha"),
                      vartypes = c("numeric","logical","numeric","character", "numeric"),
                      defaultvalues = list(NULL, FALSE, gw.cutoff,"OTP", NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE))
  utermname <- paste0("gw",bip,sp,"sp")
  termname <- paste0("dgw",bip,sp,"sp")

  decay_vs_fixed(a, termname)
  decay<-a$decay;fixed<-a$fixed
  cutoff<-a$cutoff
  decay=decay[1] # Not sure why anyone would enter a vector here, but...

  type <- .d_sp_resolve_type(nw, a, utermname, bip)
  statname <- if(type=="UTP" || nchar(bip)) utermname else paste(utermname,type,sep=".")

  if(!fixed){ # This is a curved exponential family model
    maxsp <- min(cutoff, if(bip=="b1") network.size(nw)-nw%n%"bipartite"
                          else if(bip=="b2") nw%n%"bipartite"
                          else network.size(nw)-2)
    if(maxsp==0) return(NULL)
    d <- seq_len(maxsp)

    params <- setNames(list(NULL,decay), c(statname,paste(statname,"decay",sep=".")))

    c(list(name=if(nchar(bip)) "ddspdistbwrap" else paste0("d",sp,"spdist"),
           coef.names=if(type=="UTP" || nchar(bip)) paste0(bip,sp,"sp#",d) else paste0(bip,sp,"sp.",type,"#",d),
           cutoff.message = ergm_cutoff_message(maxsp, termname, sprintf("number of %ss on some %s", SPTYPE_DESC[type], c(e="edge", d="dyad", n="nonedge")[sp]), "cutoff", "gw.cutoff"),
           iinputs=SPTYPE_CODE[type], params=params, auxiliaries=if(cache.sp) .spcache.aux(type) else NULL), GWDECAY)
  }else{
    coef.names <- paste(statname,"fixed",decay,sep=".")
    list(name=if(nchar(bip)) "dgwdspbwrap" else termname, coef.names=coef.names, inputs=decay, iinputs=SPTYPE_CODE[type], auxiliaries=if(cache.sp) .spcache.aux(type) else NULL)
  }
}

.dsp_emptynwstats <- function(nw, d, ...){
  if(any(d==0)){
    emptynwstats <- numeric(length(d))
    if(is.bipartite(nw)){
      nb1 <- nw %n% "bipartite"
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw,FALSE)
    }
    emptynwstats
  }
}

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

#' @templateVar name esp
#' @title Directed edgewise shared partners
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of edges in the network with exactly `d[i]` shared partners.
#'   
#' @usage
#' # binary: desp(d, type="OTP")
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @template ergmTerm-sp-types
#'
#' @concept directed
#' @aliases desp-ergmTerm
InitErgmTerm.desp<-function(nw, arglist, cache.sp=TRUE, ...) {
  .d_sp_impl("e", nw, arglist, cache.sp, ...)
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

#' @templateVar name gwesp
#' @title Geometrically weighted edgewise shared partner distribution
#' @description This term adds a statistic equal to the geometrically weighted edgewise (not dyadwise) shared partner distribution with decay parameter `decay` parameter.
#'   
#' @usage
#' # binary: dgwesp(decay, fixed=FALSE, cutoff=30, type="OTP")
#'
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying ESP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-sp-types
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept directed
#' @aliases dgwesp-ergmTerm
InitErgmTerm.dgwesp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  .dgw_sp_impl("e", nw, arglist, cache.sp, gw.cutoff, ...)
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
#Only one type may be specified per dsp term.

#' @templateVar name dsp
#' @title Directed dyadwise shared partners
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of dyads in the network with exactly `d[i]` shared partners.
#'   
#' @usage
#' # binary: ddsp(d, type="OTP")
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @template ergmTerm-sp-types
#'
#' @concept directed
#' @aliases ddsp-ergmTerm
InitErgmTerm.ddsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  .d_sp_impl("d", nw, arglist, cache.sp, emptynwstats=.dsp_emptynwstats, ...)
}



################################################################################

#' @templateVar name gwdsp
#' @title Geometrically weighted dyadwise shared partner distribution
#' @description This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared partner distribution with decay parameter `decay` parameter.
#'   
#' @usage
#' # binary: dgwdsp(decay, fixed=FALSE, cutoff=30, type="OTP")
#'
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying DSP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-sp-types
#'
#' @note The GWDSP statistic is equal to the sum of GWNSP plus GWESP.
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept directed
#' @aliases dgwdsp-ergmTerm
InitErgmTerm.dgwdsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  .dgw_sp_impl("d", nw, arglist, cache.sp, gw.cutoff, ...)
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

#' @templateVar name nsp
#' @title Directed non-edgewise shared partners
#' @description This term adds one network statistic to the model for each element in `d` where the \eqn{i} th such statistic equals the number of non-edges in the network with exactly `d[i]` shared partners.
#'   
#' @usage
#' # binary: dnsp(d, type="OTP")
#'
#' @param d a vector of distinct integers
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-directed
#'
#' @template ergmTerm-sp-types
#'
#' @concept directed
#' @aliases dnsp-ergmTerm
InitErgmTerm.dnsp<-function(nw, arglist, cache.sp=TRUE, ...) {
  .d_sp_impl("n", nw, arglist, cache.sp, emptynwstats=.dsp_emptynwstats, ...)
}


################################################################################

#' @templateVar name gwnsp
#' @title Geometrically weighted non-edgewise shared partner distribution
#' @description This term is just like gwesp and gwdsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with decay parameter `decay` parameter.
#'   
#' @usage
#' # binary: dgwnsp(decay, fixed=FALSE, cutoff=30, type="OTP")
#'
#' @templateVar multiplicand shared partner or selected directed analogue count
#' @template ergmTerm-gw-decay-fixed
#' @templateVar underlying NSP
#' @template ergmTerm-gw-cutoff
#' @template ergmTerm-sp-type
#'
#' @template ergmTerm-sp-types
#' @template ergmTerm-cache-sp
#' @template ergmTerm-general
#'
#' @template ergmTerm-gw-alpha-to-decay
#'
#' @concept directed
#' @aliases dgwnsp-ergmTerm
InitErgmTerm.dgwnsp<-function(nw, arglist, cache.sp=TRUE, gw.cutoff=30, ...) {
  .dgw_sp_impl("n", nw, arglist, cache.sp, gw.cutoff, ...)
}

################################################################################

#' @templateVar name gwnsp
#' @template ergmTerm-rdname
#' @usage
#' # binary: gwnsp(decay, fixed=FALSE, cutoff=30, type="OTP")
InitErgmTerm.gwnsp <- InitErgmTerm.dgwnsp

#' @templateVar name gwdsp
#' @template ergmTerm-rdname
#' @usage
#' # binary: gwdsp(decay, fixed=FALSE, cutoff=30, type="OTP")
InitErgmTerm.gwdsp <- InitErgmTerm.dgwdsp

#' @templateVar name gwesp
#' @template ergmTerm-rdname
#' @usage
#' # binary: gwesp(decay, fixed=FALSE, cutoff=30, type="OTP")
InitErgmTerm.gwesp <- InitErgmTerm.dgwesp

#' @templateVar name dsp
#' @template ergmTerm-rdname
#' @usage
#' # binary: dsp(d, type="OTP")
InitErgmTerm.dsp <- InitErgmTerm.ddsp

#' @templateVar name esp
#' @template ergmTerm-rdname
#' @usage
#' # binary: esp(d, type="OTP")
InitErgmTerm.esp <- InitErgmTerm.desp

#' @templateVar name nsp
#' @template ergmTerm-rdname
#' @usage
#' # binary: nsp(d, type="OTP")
InitErgmTerm.nsp <- InitErgmTerm.dnsp
