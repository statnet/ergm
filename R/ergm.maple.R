#  File R/ergm.maple.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
####################################################################
# The <ergm.maple> function finds a maximizer to the pseudo-
# likelihood function
#
# --PARAMETERS--
#   pl      : a list of pseudo likelihood components, as returned
#             <ergm.pl>. 
#   m       : the model, as returned by <ergm.getmodel>
#   init  : the vector of initial theta coefficients
#   MPLEtype: the method for MPL estimation as "penalized", "glm"
#             or "logitreg"; default="glm"
#   family  : the family to use in the R native routine <glm>; 
#             only applicable if "glm" is the 'MPLEtype';
#             default="binomial"
#   save.glm: whether the mple fit and the null mple fit should be
#             returned (T or F); if false, NULL is returned for both;
#             default==TRUE
#   thetal  : the independence theta; if specified and non-NULL, this
#             is ignored except to return its value in the returned 
#             ergm; default=NULL, in which case 'theta1' is computed           
#   verbose : whether this and the C routines should be verbose
#             (T or F); default=FALSE
#   ...     : additional parameters passed from within; all will
#             be ignored                                                           
#
# --RETURNED--
#   an ergm object as a list containing several items; for details
#   see the return list in the <ergm> function header
#   (<ergm.maple>= ~)
#
##################################################################

ergm.maple<-function(pl, m, init=NULL,
                    MPLEtype="glm", family="binomial",
                    save.glm=TRUE,
                    theta1=NULL, verbose=FALSE, ...) {
  if(MPLEtype=="penalized"){
   if(verbose) message("Using penalized MPLE.")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend,
                           start=init)
#  mple$deviance <- 2 * (mplefit$loglik-mplefit$loglik[1])[-1]
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
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
                    independent=independent)
     }else{
      mindfit.summary <- summary(mindfit)
      theta.ind[independent] <- mindfit$coef
      theta1 <- list(coef=mindfit$coef, 
                    theta=theta.ind,
                    independent=independent)
     }
    }else{
     theta1 <- list(coef=NULL, 
                    theta=rep(0,ncol(pl$xmat)),
                    independent=independent)
    }
   }
#

   if(nrow(pl$xmat) > pl$maxMPLEsamplesize){
#
#   fix deviance for sampled data
#
    mplefit$deviance <- ergm.logisticdeviance(beta=mplefit$coef,
     X=model.matrix(terms(pl$zy.full ~ .-1,data=data.frame(pl$xmat.full)),
                           data=data.frame(pl$xmat.full)),
     y=pl$zy.full, offset=pl$foffset.full)
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

# mplefit <- call(MPLEtype, pl$zy ~ 1, family=binomial)
#
  if(MPLEtype=="penalized"){
   mplefit.null <- ergm.pen.glm(pl$zy ~ 1, weights=pl$wend)
  }else{
   if(MPLEtype=="logitreg"){
    mplefit.null <- ergm.logitreg(x=matrix(1,ncol=1,nrow=length(pl$zy)),
                                  y=pl$zy, offset=pl$foffset, wt=pl$wend)
   }else{
    mplefit.null <- try(glm(pl$zy ~ 1, family=family, weights=pl$wend),
                        silent = TRUE)
    if (inherits(mplefit.null, "try-error")) {
      mplefit.null <- list(coef=0, deviance=0, null.deviance=0,
                      cov.unscaled=diag(1))
    }
   }
  }

  if(save.glm){
    glm <- mplefit
    glm.null <- mplefit.null
  }else{
    glm <- NULL
    glm.null <- NULL
  }

# Output results as ergm-class object
  structure(list(coef=theta,
      iterations=iteration,
      MCMCtheta=theta, gradient=gradient,
      hessian=NULL, covar=covar, failure=FALSE,
      mc.se=mc.se, glm = glm, glm.null = glm.null,
      theta1=theta1),
     class="ergm")
}


