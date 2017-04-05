#  File ergm/R/ergm.reviseinit.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <ergm.reviseinit> function revises 'init' to reflect additional
# parameters introduced by curved model terms
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
      















































































