#  File R/InitErgmTerm.spcache.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

.spcache.aux <- function(type){
  type <- toupper(type)
  trim_env(as.formula(as.call(list(as.name('~'), as.call(list(as.name('.spcache.net'),type=if(type=='ITP')'OTP' else type))))))
}

#' @templateVar name .spcache.net
#' @title Shared Partner Cache
#' @description This auxiliary maintains a hash table mapping dyads to
#'   shared partner counts.
#'
#' @usage
#' # binary: .spcache(type)
#'
#' @template ergmTerm-sp-type
#'
#' @details This auxiliary maintains a single \code{StrictDyadMapUInt}
#'   defined in \file{"ergm_dyad_hashmap.h"}. It is, internally, a
#'   \code{khash} table mapping a pair of \code{Vertex}es (in a
#'   \code{TailHead} structure) onto an \code{unsigned int}.
#'
#' To make use of it in your change statistic:
#'
#' 1. Add `auxiliaries = ~.spcache.net(type)` to the `InitErgmTerm`
#'    output list. You may request it multiple times and/or alongside
#'    other auxiliaries.
#'
#' 2. Add `#include "ergm_dyad_hashmap.h"` to the top of your change
#'    statistic implementation file.
#'
#' 3. In the function (e.g., `c_` function) implementing the change
#'    statistic, use \code{GET_AUX_STORAGE(StoreStrictDyadMapUInt,
#'    spcache);} to obtain it if this is your first or only auxiliary,
#'    or \code{GET_AUX_STORAGE(ind, StoreStrictDyadMapUInt, spcache);} if it is not.
#'
#' 4. Use one of the following macros to access the shared partner count: \describe{
#' \item{\code{GETUDMUI(\var{TAIL}, \var{HEAD}, spcache)}}{if undirected;}
#' \item{\code{GETDDMUI(\var{HEAD}, \var{TAIL}, spcache)}}{if `type = "ITP"` (since ITP is OTP with direction reversed);}
#' \item{\code{GETDDMUI(\var{TAIL}, \var{HEAD}, spcache)}}{in all other cases.}
#' }
#'
#' @template ergmTerm-sp-types
#'
#' @template ergmAuxiliary-general
#'
#' @references \insertAllCited{}
#'
#' @concept directed
#' @concept undirected
#' @concept triad-related
InitErgmTerm..spcache.net<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("type"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  type <- match.arg(tolower(a$type), c("otp","osp","isp","utp","rtp")) # ITP not included, because it's just OTP with direction reversed.

  if (is.directed(nw) == (type == "utp")
      && !(is.bipartite(nw) && type %in% c("osp", "isp")))
    stop("Type UTP may only be used with undirected networks, OSP and ISP with bipartite or directed, and the rest only with directed.")

  list(name=paste0("_",type,"_wtnet"),
       coef.names=c(), dependence=TRUE)
}
