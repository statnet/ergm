#excerpted from gof.ergm.R:

##FIXME: Add support for curved ERGMs.
##FIXME: Merge with gof.ergm and/or with ergm itself.
##FIXME: Publish this somewhere.


######################################################################
# The <logp.rank.ergm> function returns the percent of the eta sums
# which are negative and plots the density of the eta sums if desired
# (eta sums are the sum of the stats * the estimated coefficients)
#
# --PARAMETERS--
#   x   : an ergm object
#   plot: whether to plot the density of the ; default=FALSE
#
# --RETURNED--
#   the percent of eta sums which are negative
#
######################################################################

logp.rank.ergm<-function(x,plot=FALSE){
  if(length(x$etamap$curved))
    stop("Curved exponential families are not supported at this time.")
  etasum<-c(x$sample %*% x$coef)
  if(plot){
    plot(density(etasum),main=expression(paste("Density of ",log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y["obs"],";",hat(theta)))))),
         xlab=expression(log(Pr(paste(Y,";",hat(theta)))/Pr(paste(y["obs"],";",hat(theta))))),zero.line=FALSE)
    abline(v=0)
  }
  mean(etasum<0)
}
