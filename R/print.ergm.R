#  File R/print.ergm.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
###############################################################################
# The <print.ergm> function prints summary information for a given ergm
#
# --PARAMETERS--
#   x     :  an ergm object
#   digits:  the number of significant digits for the coefficients;
#            default=max(3, getOption("digits")-3)
#   ...   :  additional parameters passed from within; these will be ignored
#
# --RETURNED--
#   x
# 
###############################################################################

#' @describeIn ergm Print the call, the estimate, and the method used to obtain it.
#'
#' @param x,digits See [print()].
#'
#' @export
print.ergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
  # The following code is based on stats:::print.lm(), but there really isn't another concise way to do this:
  if(!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse="\n"), "\n", sep="")

  cat("\n",x$estimate.desc," Coefficients:\n",sep="")
  print.default(format(coef(x), digits = digits), print.gap = 2, ...,
                quote = FALSE)

  invisible(x)
}
