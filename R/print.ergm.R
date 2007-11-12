print.ergm <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
  if(!is.latent(x) || is.null(x$Z.mle))
  {
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
  else
  {
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
