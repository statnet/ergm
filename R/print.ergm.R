#  File ergm/R/print.ergm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
print.ergm <- function (x, digits = max(3, getOption("digits") - 3), ...) {
#  if(!is.latent(x) || is.null(x$Z.mle)) {                             
  if(is.null(x$Z.mle)) {
   if(is.matrix(x$sample)){
    if(!is.matrix(x$thetasample)){
     cat("Newton-Raphson iterations: ", x$iterations[1], "\n")
    }
    cat("MCMC sample of size", x$samplesize, "based on: \n")
    print.default(format(x$MCMCtheta, digits = digits), print.gap = 2, 
        quote = FALSE)
    cat("\nMonte Carlo MLE Coefficients:\n")
    print.default(format(x$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
   }else{
    cat("Newton-Raphson iterations: ", x$iterations[1], "\n")
    if(!all(x$theta1$independent)){
     cat("\nMPLE Coefficients:\n")
     print.default(format(x$coef, digits = digits), print.gap = 2, 
         quote = FALSE)
    }else{
    cat("\nMLE Coefficients:\n")
    print.default(format(x$coef, digits = digits), print.gap = 2, 
        quote = FALSE)
    }
   }
  }
  else {
    cat("Newton-Raphson iterations: ", x$iterations[1], "\n")
    cat("MLE and Posterior Mean of Coefficients:\n")
    temp <- cbind(x$MCMCtheta,x$coef)
    colnames(temp) <- c("MLE","MAP")
    print.default(format(temp, digits = digits), print.gap = 2,
                  quote = FALSE)
    scale.MLE <-  0
    for(i in 1:ncol(x$Z.mle))
      scale.MLE <- scale.MLE + sum(outer(x$Z.mle[,i],x$Z.mle[,i],function(x,y)abs(x-y)))
    scale.MLE <- scale.MLE / choose(nrow(x$Z.mle),2)
    cat("Scale = ",scale.MLE,"\n")
    n <- network.size(x$newnetwork)
    if(!is.null(x$Z.pmode))
      p <- ncol(x$Z.pmode)
    else
      p <- 0
    r <- length(x$coef) + (n) * p - p*(p+1)/2
    n <- n*(n-1)
    cat("\nBIC = ",2 * x$mle.lik - r*log(n),"\n")
  }
  invisible(x)
}
