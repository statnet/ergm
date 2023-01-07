#' @describeIn ergm Return `TRUE` if the ERGM was fit to a partially observed network and/or an observational process, such as missing (`NA`) dyads.
#' @export
is.na.ergm <- function(x){
  NVL(x$info$obs, !is.null(x$constrained.obs))
}

#' @describeIn ergm Alias to the `is.na()` method.
#' @export
anyNA.ergm <- function(x, ...){
  is.na(x)
}
