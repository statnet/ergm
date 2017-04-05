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

