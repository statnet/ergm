#  File ergm/R/ergm.maple.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
####################################################################
# The <ergm.maple> function finds a maximizer to the pseudo-
# likelihood function
##################################################################

ergm.maple<-function(pl, m, init=NULL,
                    MPLEtype="glm", family="binomial",
                    save.glm=TRUE,
                    theta1=NULL, verbose=FALSE, ...) {
  if(MPLEtype=="penalized"){
   if(verbose) cat("Using penalized MPLE.\n")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend,
                           start=init)
#  mple$deviance <- 2 * (mplefit$loglik-mplefit$loglik[1])[-1]
   mplefit$deviance <- -2*mplefit$loglik
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEtype=="logitreg"){
    mplefit <- model.matrix(terms(pl$zy ~ .-1,data=data.frame(pl$xmat)),
                           data=data.frame(pl$xmat))
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, offset=pl$foffset, wt=pl$wend,
                             start=init)
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
    mplefit <- try(
          glm(pl$zy ~ .-1 + offset(pl$foffset), data=data.frame(pl$xmat),
                    weights=pl$wend, family=family, start=init),
                    silent = TRUE)
    if (inherits(mplefit, "try-error")) {
      mplefit <- list(coef=pl$theta.offset, deviance=0,
                      cov.unscaled=diag(length(pl$theta.offset)))
      mplefit.summary <- list(cov.unscaled=diag(length(pl$theta.offset)))
    }else{
      mplefit.summary <- summary(mplefit)
    }
   }
#
#  Determine the independence theta and MLE
#  Note that the term "match" is deprecated.
#
   if(is.null(theta1)){
    independent.terms <- 
       c("edges","match","nodecov","nodefactor","nodematch","absdiff",
         "edgecov","dyadcov","sender","receiver","sociality", 
         "nodeofactor","nodeifactor","nodemain",
         "nodemix","mix",
         "b1","b2",
         "testme")
    independent <- rep(0,ncol(pl$xmat))
    names(independent) <- colnames(pl$xmat)
    theta.ind <- independent
    for(i in seq(along=independent.terms)){
     independent[grep(independent.terms[i], colnames(pl$xmat))] <- i
    }
    independent <- independent>0
    if(any(independent)){
     mindfit <- try(glm(pl$zy ~ .-1 + offset(pl$foffset), 
                    data=data.frame(pl$xmat[,independent,drop=FALSE]),
                    weights=pl$wend, family=family),
                    silent = TRUE)
     if (inherits(mindfit, "try-error")) {
      theta1 <- list(coef=NULL, 
                    theta=rep(0,ncol(pl$xmat)),
                    independent=independent,
                    loglikelihood=-pl$numobs*log(2))
     }else{
      mindfit.summary <- summary(mindfit)
      theta.ind[independent] <- mindfit$coef
      theta1 <- list(coef=mindfit$coef, 
                    theta=theta.ind,
                    independent=independent,
                    loglikelihood=-mindfit$deviance/2)
     }
    }else{
     theta1 <- list(coef=NULL, 
                    theta=rep(0,ncol(pl$xmat)),
                    independent=independent,
                    loglikelihood=-pl$numobs*log(2))
    }
   }
#
   options(warn=0)
#  options(warn=2)
   if(nrow(pl$xmat) > pl$maxMPLEsamplesize){
#
#   fix aic and deviance for sampled data
#
    mplefit$deviance <- ergm.logisticdeviance(beta=mplefit$coef,
     X=model.matrix(terms(pl$zy.full ~ .-1,data=data.frame(pl$xmat.full)),
                           data=data.frame(pl$xmat.full)),
     y=pl$zy.full, offset=pl$foffset.full)
    mplefit$aic <- mplefit$deviance + 2*mplefit$rank
   }
  }
  theta <- pl$theta.offset
  real.coef <- mplefit$coef
  real.cov <- mplefit.summary$cov.unscaled
  theta[!m$etamap$offsettheta] <- real.coef
# theta[is.na(theta)] <- 0
  names(theta) <- m$coef.names

#
# Old end
#
  gradient <- rep(NA, length(theta))
#
# Calculate the (global) log-likelihood
#
  loglik <- -mplefit$deviance/2
#
  mc.se <- gradient <- rep(NA, length(theta))
  if(length(theta)==1){
   covar <- array(0,dim=c(1,1))
  }else{
   covar <- diag(rep(0,length(theta)))
  }
# covar <- as.matrix(covar[!m$etamap$offsettheta,!m$etamap$offsettheta])
# covar[!is.na(real.coef),!is.na(real.coef)] <- real.cov
  covar[!is.na(theta)&!m$etamap$offsettheta,!is.na(theta)&!m$etamap$offsettheta] <- real.cov
#
  iteration <-  mplefit$iter 
  samplesize <- NA

# mplefit <- call(MPLEtype, pl$zy ~ 1, family=binomial)
#
  if(MPLEtype=="penalized"){
   mplefit.null <- ergm.pen.glm(pl$zy ~ 1, weights=pl$wend)
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEtype=="logitreg"){
    mplefit.null <- ergm.logitreg(x=matrix(1,ncol=1,nrow=length(pl$zy)),
                                  y=pl$zy, offset=pl$foffset, wt=pl$wend)
   }else{
    mplefit.null <- try(glm(pl$zy ~ 1, family=family, weights=pl$wend),
                        silent = TRUE)
    if (inherits(mplefit.null, "try-error")) {
      mplefit.null <- list(coef=0, deviance=0,
                      cov.unscaled=diag(1))
    }
   }
   options(warn=0)
#  options(warn=2)
  }

  null.deviance <- mplefit$null.deviance
  aic <- mplefit$aic

  if(save.glm){
    glm <- mplefit
    glm.null <- mplefit.null
  }else{
    glm <- NULL
    glm.null <- NULL
  }

# Output results as ergm-class object
  structure(list(coef=theta, sample=NA,
      iterations=iteration, mle.lik=loglik,
      MCMCtheta=theta, loglikelihoodratio=loglik, gradient=gradient,
      hessian=NULL, covar=covar, samplesize=samplesize, failure=FALSE,
      mc.se=mc.se, glm = glm, glm.null = glm.null,
      null.deviance=null.deviance, aic=aic,
      theta1=theta1),
     class="ergm")
}


