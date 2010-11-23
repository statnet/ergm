################################################################################
# The <ergm.independencemodel> function checks whether the ergm model is a 
# dyadic independence ergm.
#
# --PARAMETERS--
#   m:  the model object, as returned by <ergm.getmodel>
#
# --RETURNED--
#   ans: whether every term in the model is dyadic independent (T or F)
#
###############################################################################

ergm.independencemodel <- function(m) {
  ans <- TRUE
  for (i in 1:length(m$terms)) {
    if(is.null(m$terms[[i]]$dependence) || m$terms[[i]]$dependence)
      ans <- FALSE
  }
  ans
}
