ergm.independencemodel <- function(m) {
# Return TRUE or FALSE, depending on whether the model object m has
# every term with dependence=FALSE set.

  ans <- TRUE
  for (i in 1:length(m$terms)) {
    if(is.null(m$terms[[i]]$dependence) || m$terms[[i]]$dependence)
      ans <- FALSE
  }
  ans
}
