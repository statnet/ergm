#  File ergm/R/ergmMPLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
############################################################################
# The <ergmMPLE> function has different behavior based on whether the given
# formula should be fit or not. If so, <ergm> is called. If not, the elements
# needed for logit regression are computed and returned using <ergm.pl>
###########################################################################

ergmMPLE <- function(formula, fitmodel=FALSE, control=control.ergm(),
                     verbose=FALSE, ...) 
{
  if (fitmodel) {
    return(ergm(formula, estimate="MPLE", control=control, verbose=verbose, ...))
  }
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  MPLEsetup <- ergm.pl(Clist, Clist.miss, model, verbose=verbose, ...)
  list(response = MPLEsetup$zy, predictor = MPLEsetup$xmat, 
       weights = MPLEsetup$wend)
}

