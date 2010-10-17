#excerpted from gof.ergm.R:

##FIXME: Add support for curved ERGMs.
##FIXME: Merge with gof.ergm and/or with ergm itself. 
logp.rank.ergm<-function(x,plot=FALSE){
  if(length(x$etamap$curved))
    stop("Curved exponential families are not supported at this time.")
  etasum<-c(x$sample %*% x$coef)
  if(plot){
    plot(density(etasum),main=expression(paste("Density of ",log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y,";",hat(theta)))))),
         xlab=expression(log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y,";",hat(theta))))),zero.line=FALSE)
    abline(v=0)
  }
  mean(etasum<0)
}
