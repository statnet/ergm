#  File R/InitErgmTerm.spcache.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

InitErgmTerm..spcache.net<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("type"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  type <- match.arg(tolower(a$type), c("otp","osp","isp","utp","rtp")) # ITP not included, because it's just OTP with direction reversed.

  if(is.directed(nw)==(type=="utp") && !(NVL(nw%n%"bipartite",0)>0 && type%in%c("osp","isp"))) stop("Type UTP may only be used with undirected networks, OSP and ISP with bipartite or directed, and the rest only with directed.")
  
  list(name=paste0("_",type,"_wtnet"),
       coef.names=c(), dependence=TRUE)
}
