#  File ergm/R/ergm.independencemodel.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
################################################################################
# The <ergm.independencemodel> function checks whether the ergm model is a 
# dyadic independence ergm.
###############################################################################

ergm.independencemodel <- function(m) {
  ans <- TRUE
  for (i in 1:length(m$terms)) {
    if(is.null(m$terms[[i]]$dependence) || m$terms[[i]]$dependence)
      ans <- FALSE
  }
  ans
}
