#  File R/print.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
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



#' @describeIn ergm
#'
#' @param x,digits See [print()].
#' 
#' Automatically called when an object of class \code{\link{ergm}} is printed.
#' Currently, summarizes the size of the MCMC sample, the \eqn{\theta}
#' vector governing the selection of the sample, and the Monte Carlo MLE. The optional `digits` argument specifies the significant digits for coefficients
#' @export
print.ergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
   if(is.matrix(x$sample)){
#    if(!is.matrix(x$thetasample) && !is.null(x$iterations)){
#     cat("Newton-Raphson iterations: ", x$iterations[1], "\n")
#    }
    cat("MCMC sample of size", nrow(as.matrix(x$sample)), "based on: \n")
    print.default(format(x$MCMCtheta, digits = digits), print.gap = 2, 
        quote = FALSE)
    cat("\nMonte Carlo MLE Coefficients:\n")
    print.default(format(x$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
   }else{
#    if (!is.null(x$iterations)) {
#      cat("Newton-Raphson iterations: ", x$iterations[1], "\n")
#    }

     cat("\n",x$estimate," Coefficients:\n",sep="")
     print.default(format(x$coef, digits = digits), print.gap = 2, 
         quote = FALSE)
   }
  invisible(x)
}
