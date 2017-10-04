#' A convenience container for a list of [`network`] objects, output
#' by \code{\link{simulate.ergm}} among others.
#'
#' @param object,x a list of networks or a `network.list` object.
#' @param ... for `network.list`, additional attributes to be set on
#'   the network list; for others, arguments passed down to
#'   lower-level functions.
#'
#' @export network.list
network.list <- function(object,...){
  if(any(!sapply(list(object), is.network))) stop("network.list() takes a list of networks as its first argument.")
  ddd <- list(...) # FIXME: Use alist here?
  ns <- names(ddd)
  for(i in seq_along(ddd)){
    attr(object, ns[i]) <- ddd[[i]]
  }
  class(object) <- c("network.list","list")
  object
}
