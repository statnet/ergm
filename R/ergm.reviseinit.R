#  File R/ergm.reviseinit.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
###############################################################################
# The <ergm.reviseinit> function revises 'init' to reflect additional
# parameters introduced by curved model terms
#
# --PARAMETERS--
#   m     :  the model, as returned by <ergm.getmodel>
#   init:  the vector of initial theta parameters
#
#
# --RETURNED--
#   init:  the revised 'init'; it is assumed for each term that the
#            parameters in 'm$terms[[j]]$params' matches the terms in 'init'
#            only in a contiguous region of 'init'; this region of overlap
#            is expanded to include any new terms in 'm$terms[[j]]$params',
#            along with their values; for details about 'm$terms', see the
#            <InitErgm> function header
#
###############################################################################

ergm.reviseinit <- function(m, init) {
  for(i in 1:length(m$terms)) {
    if (!is.null(m$terms[[i]]$params)) {
      n1=names(m$terms[[i]]$params)
      n2=names(init)
      overlap = (n2 %in% n1)
      if (length(n1)>sum(overlap)) {
        before = (cumsum(overlap)==0)
        after = (!before & overlap==0)
        newoverlap=c(init[overlap],rep(0,length(n1)-sum(overlap)))
        names(newoverlap)=n1
        for (j in (1+sum(overlap)):length(n1)) {
          newoverlap[j] = m$terms[[i]]$params[[j]]
        }
        init=c(init[before],newoverlap,init[after])
      }
    }
  }
  init
}
      















































































