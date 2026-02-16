#  File R/ergm-deprecated.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' @name ergm-deprecated
#' @rdname ergm-deprecated
#' @title Functions that will no longer be supported in future releases of the package
#' @description Functions that have been superceed, were never documented, or will be removed from the package for other reasons
#' @param ...,object,x
#' Arguments to deprecated functions.
#'
#' @keywords misc internal
NULL


#' @describeIn ergm-deprecated extracts the `ergm` parameters; may be
#'   removed in favour of the default method once the number of `ergm`
#'   objects with `$coef` elements in the wild is sufficiently low.
#' @export
coef.ergm <- function(object, ...) {
  if ("coef" %in% names(object)) unclass(object)$coef
  else NextMethod()
}


#' @describeIn ergm-deprecated accesses elements of `ergm` objects;
#'   needed for backwards compatibility when components get renamed.
#' @param name See [Extract].
#' @export
`$.ergm` <- function(x, name) {
  if (name == "coef") {
    .Deprecate_once(msg = "Using x$coef to access the coefficient vector of an ergm is deprecated. Use coef(x) instead.")
    coef(x)
  } else {
    x[[name, exact = FALSE]]
  }
}
