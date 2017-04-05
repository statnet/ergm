############################################################################
# The <ergmMPLE> function has different behavior based on whether the given
# formula should be fit or not. If so, <ergm> is called. If not, the elements
# needed for logit regression are computed and returned using <ergm.pl>
#
# --PARAMETERS--
#   formula : a formula as 'nw ~ term(s)'
#   fitmodel: whether to fit the model given by 'formula' (T or F);
#             default=FALSE
#   control : a list of parameters to control the fitting process; this
#             is ignored if 'fitmodel'=FALSE
#   verbose : whether the <ergm> or <ergm.pl> functions should be verbose
#             (T or F); default=FALSE
#   ..      : additional parameters that will be passed onto <ergm> or
#             <ergm.pl>
#
# --RETURNED--
#   if 'fitmodel'
#     =TRUE  -- an ergm object, as returned by <ergm>
#     =FALSE -- a list with 3 components:
#                response : the vector of dyad values; this is tabulated
#                           according to 'weights' 
#                predictor: the design matrix of change stats;  this is 
#                           tabulated according to 'weights' 
#                weights  : the weights for each entry/row of 'response'/
#                           'predictor'
#
###########################################################################

ergmMPLE <- function(formula, fitmodel=FALSE, control=control.ergm(),
                     verbose=FALSE, ...) 
{
  if (fitmodel) {
    return(ergm(formula, MPLEonly=TRUE, control=control, verbose=verbose, ...))
  }
  nw <- ergm.getnetwork(formula)
  model <- ergm.getmodel(formula, nw, drop=FALSE, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, model)
  Clist.miss <- ergm.design(nw, model, verbose=verbose)
  MPLEsetup <- ergm.pl(Clist, Clist.miss, model, verbose=verbose, ...)
  list(response = MPLEsetup$zy, predictor = MPLEsetup$xmat, 
       weights = MPLEsetup$wend)
}

