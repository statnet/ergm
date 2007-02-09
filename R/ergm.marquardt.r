
ergm.marquardt<-function() {
  outofbounds <- (apply(statsmatrix,2,max)*
                  apply(statsmatrix,2,min)) > -epsilon
#    if(length(geod) >0){outofbounds <- outofbounds & !geodf}
#    if(length(geosd)>0){outofbounds <- outofbounds & !geosdf}
#    if(length(wdeg) > 0){outofbounds <- outofbounds & !wdegf}
#    if(length(gwesp)> 0){outofbounds <- outofbounds & !gwespf}

  if (any(outofbounds)){
#
#   Some statistics out-of-bounds so print out summary and
#   use Marquardt adjustment
#

#
#     Print out the summary
#
#     Reconstruct the statsmatrix
#
      sm <- ergm.curved.statsmatrix(statsmatrix,theta0,parms.curved)$sm
##      cat("Summary of the simulation:\n")
##      print(apply(sm,2,summary.statsmatrix.ergm),scipen=6)

    cat(paste("Warning: No MCMC-MLE exists!\n",
              "Observed network statistic outside or on boundary of\n",
              "range found in MCMC sample for these statistics:\n",
              paste(names(outofbounds)[outofbounds],collapse=", ")))
    cat(paste("\nApplying stochastic approximation with Marquardt adjustment.\n"))
    
#
#   RM approximation
#
#   Marquardt check on the merit function
#
    E <- apply(statsmatrix0,2,mean)
    H <- cov(statsmatrix0)
    gradsq <- sum(E*E/diag(H))
#
#   Marquardt inflation
#
##     cat("gradsq: old",marquardt$gradsq0," new ",gradsq,"\n")
##     if(gradsq < marquardt$gradsq0){
#
#   Accept proposed value
#
##  marquardt$lambda <- marquardt$lambda/10
    diag(H) <- diag(H)*(1+marquardt$lambda)
    cat("lambda",marquardt$lambda,"\n")
    covarH <- robust.inverse(H)
    H <- cov(statsmatrix0)
    diffcoef <- theta0-theta0
    if(!is.null(offset.fix)){
      diffcoef[!offset.fix] <- as.vector(covarH %*% E)
      coef <- theta0
      coef[!offset.fix] <- coef[!offset.fix] - diffcoef
      cat("New theta",coef,"\n")
      cat("diff",diffcoef,"\n")
    }else{
      diffcoef <- as.vector(covarH %*% E)
      coef <- theta0 - diffcoef
      cat("New theta",coef,"\n")
      cat("diff",diffcoef,"\n")
    }
    if(length(theta0)==1){
      covar <- array(0,dim=c(1,1))
    }else{
      covar <- diag(rep(0,length(theta0)))
    }
    if(!is.null(offset.fix)){
      covar[!offset.fix, !offset.fix] <- covarH
    }
#
#   return(list(failure=TRUE, 
#               coef=guess,
#               MCMCtheta=guess,
#               mc.se=rep(NA,length(guess)),
#               mcmcloglik=0,
#               hessian=diag(guess),
#               acf=ergm.MCMCacf(statsmatrix0)))
    return(list(failure=FALSE, 
                coef=coef,
                MCMCtheta=guess,
                mc.se=rep(NA,length(theta0)),
                samplesize=NA,
                sample=statsmatrix0,
#               sample=matrix(0,ncol=length(guess),nrow=1),
                loglikelihood=NA,
                mcmcloglik=NA,
                covar=covar,
                acf=ergm.MCMCacf(statsmatrix0)))
  }
}


ergm.marquardt2 = function() {
#      Marquardt check on the merit function
#
  E <- apply(statsmatrix0,2,mean)
  H <- cov(statsmatrix0)
  gradsq <- sum(E*E/diag(H))
#
#      Marquardt inflation
#
  diag(H) <- diag(H)*(1+marquardt$lambda)
  cat("lambda",marquardt$lambda,"\n")
  diffcoef <- theta0-theta0
  covarH <- robust.inverse(H)
  diffcoef <- as.vector(covarH %*% E)
  coef <- theta0 - diffcoef
  cat("New theta",coef,"\n")
  cat("diff",diffcoef,"\n")
  
  if(length(theta0)==1){
    covar <- array(0,dim=c(1,1))
  }else{
    covar <- diag(rep(0,length(theta0)))
  }
  covar <- covarH
  return(list(failure=TRUE, 
              coef=coef,
              MCMCtheta=theta0,
              mc.se=rep(NA,length(theta0)),
              mcmcloglik=0,
              samplesize=0,
              sample=statsmatrix0,
#               sample=matrix(0,ncol=length(guess),nrow=1),
              loglikelihood=NA,
              covar=covar,
              acf=ergm.MCMCacf(statsmatrix0)))
#      return(list(failure=TRUE))
}
