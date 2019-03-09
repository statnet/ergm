#  File R/ergm.mple.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################

#' Find a maximizer to the psuedolikelihood function
#' 
#' The \code{ergm.mple} function finds a maximizer to the psuedolikelihood
#' function (MPLE). It is the default method for finding the ERGM starting
#' coefficient values. It is normally called internally the ergm process and
#' not directly by the user. Generally \code{\link{ergmMPLE}} would be called
#' by users instead.
#' 
#' According to Hunter et al. (2008): "The maximizer of the pseudolikelihood
#' may thus easily be found (at least in principle) by using logistic
#' regression as a computational device." In order for this to work, the
#' predictors of the logistic regression model must be calculated.  These are
#' the change statistics as described in Section 3.2 of Hunter et al. (2008),
#' put into matrix form so that each pair of nodes is one row whose values are
#' the vector of change statistics for that node pair.  The ergm.pl function
#' computes these change statistics and the ergm.mple function implements the
#' logistic regression using R's glm function.  Generally, neither ergm.mple
#' nor ergm.pl should be called by users if the logistic regression output is
#' desired; instead, use the \code{\link{ergmMPLE}} function.
#' 
#' In the case where the ERGM is a dyadic independence model, the MPLE is the
#' same as the MLE.  However, in general this is not the case and, as van Duijn
#' et al. (2009) warn, the statistical properties of MPLEs in general are
#' somewhat mysterious.
#' 
#' MPLE values are used even in the case of dyadic dependence models as
#' starting points for the MCMC algorithm.
#' 
#' @param nw response network.
#' @param fd An \code{\link{rlebdm}} with informative dyads.
#' @param m the model, as returned by \code{\link{ergm_model}}
#' @param init a vector a vector of initial theta coefficients
#' @param MPLEtype the method for MPL estimation as "penalized", "glm"
#'   or "logitreg"; default="glm"
#' @param family the family to use in the R native routine
#'   \code{\link{glm}}; only applicable if "glm" is the 'MPLEtype';
#'   default="binomial"
#' @param maxMPLEsamplesize the sample size to use for endogenous
#'   sampling in the psuedo-likelihood computation; default=1e6
#' @param save.glm whether the mple fit and the null mple fit should
#'   be returned (T or F); if false, NULL is returned for both;
#'   default==TRUE
#' @param theta1 the independence theta; if specified and non-NULL,
#'   this is ignored except to return its value in the returned ergm;
#'   default=NULL, in which case 'theta1' is computed
#' @param control a list of MCMC related parameters; recognized
#'   components include: samplesize : the number of networks to sample
#'   Clist.miss : see 'Clist.miss' above; some of the code uses this
#'   Clist.miss,
#' @param proposal an [ergm_proposal()] object.
#' @param verbose whether this and the C routines should be verbose (T
#'   or F); default=FALSE
#' @param \dots additional parameters passed from within; all will be
#'   ignored
#' @return \code{ergm.mple} returns an ergm object as a list
#'   containing several items; for details see the return list in the
#'   \code{\link{ergm}}
#' 
#' @seealso \code{\link{ergmMPLE}},
#' \code{\link{ergm}},\code{\link{control.ergm}}
#' @references Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris and
#' Martina (2008).  "ergm: A Package to Fit, Simulate and Diagnose
#' Exponential-Family Models for Networks." _Journal of Statistical Software_,
#' *24*(3), pp. 1-29. \url{https://www.jstatsoft.org/article/view/v024i03}
#' 
#' van Duijn MAJ, Gile K, Handcock MS (2009).  "Comparison of Maximum Pseudo
#' Likelihood and Maximum Likelihood Estimation of Exponential Family Random
#' Graph Models." _Social Networks_, *31*, pp. 52-62.
ergm.mple<-function(nw, fd, m, init=NULL,
                    MPLEtype="glm", family="binomial",
                    maxMPLEsamplesize=1e+6,
                    save.glm=TRUE,
                    theta1=NULL, 
		    control=NULL, proposal=NULL,
                    verbose=FALSE,
                    ...) {
  message("Starting maximum pseudolikelihood estimation (MPLE):")
  message("Evaluating the predictor and response matrix.")
  pl <- ergm.pl(nw=nw, fd=fd, m=m,
                theta.offset=init,
                maxMPLEsamplesize=maxMPLEsamplesize,
		control=control,
                verbose=verbose)

  message("Maximizing the pseudolikelihood.")
  if(MPLEtype=="penalized"){
   if(verbose) message("Using penalized MPLE.")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend,
                  start=init[!m$etamap$offsettheta])
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
   if(MPLEtype=="logitreg"){
    mplefit <- model.matrix(terms(pl$zy ~ .-1,data=data.frame(pl$xmat)),
                           data=data.frame(pl$xmat))
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, offset=pl$foffset, wt=pl$wend,
                             start=init[!m$etamap$offsettheta])
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
#     mplefit <- suppressWarnings(try(
#           glm(pl$zy ~ .-1 + offset(pl$foffset), data=data.frame(pl$xmat),
#                weights=pl$wend, family=family),
# # Note:  It appears that specifying a starting vector can lead to problems!
# #               start=init[!m$etamap$offsettheta]),
#                     silent = TRUE))
    glm.result <- .catchToList(glm(pl$zy ~ .-1 + offset(pl$foffset), 
                                  data=data.frame(pl$xmat),
                                  weights=pl$wend, family=family))
    
    # error handling for glm results
    if (!is.null(glm.result$error)) {
      warning(glm.result$error)
      mplefit <- list(coef=pl$theta.offset, deviance=0,
                      cov.unscaled=diag(length(pl$theta.offset)))
      mplefit.summary <- list(cov.unscaled=diag(length(pl$theta.offset)))
    } else if (!is.null(glm.result$warnings)) {
      # if the glm results are crazy, redo it with 0 starting values
      if (max(abs(glm.result$value$coef), na.rm=T) > 1e6) {
        warning("GLM model may be separable; restarting glm with zeros.\n")
        mplefit <- glm(pl$zy ~ .-1 + offset(pl$foffset), 
                       data=data.frame(pl$xmat),
                       weights=pl$wend, family=family, 
                       start=rep.int(0, length(init[!m$etamap$offsettheta])))
        mplefit.summary <- summary(mplefit)
      } else {
        # unknown warning, just report it
        warning(glm.result$warnings)
        mplefit <- glm.result$value
        mplefit.summary <- summary(mplefit)
      }
    } else {
      # no errors or warnings
      mplefit <- glm.result$value
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
         "nodeofactor","nodeifactor","nodemain",
         "edgecov","dyadcov","sender","receiver","sociality", 
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
      
      glm.result <- .catchToList(glm(pl$zy ~ .-1 + offset(pl$foffset), 
                                    data=data.frame(pl$xmat[,independent,drop=FALSE]),
                                    weights=pl$wend, family=family))
      
      # error handling for glm results
      if (!is.null(glm.result$error)) {
        warning(glm.result$error)
        theta1 <- list(coef=NULL, 
                       theta=rep(0,ncol(pl$xmat)),
                       independent=independent)
      } else if (!is.null(glm.result$warnings)) {
        # if the glm results are crazy, redo it with 0 starting values
        if (max(abs(glm.result$value$coef), na.rm=T) > 1e6) {
          warning("GLM model may be separable; restarting glm with zeros.\n")
          mindfit <- glm(pl$zy ~ .-1 + offset(pl$foffset), 
                         data=data.frame(pl$xmat[,independent,drop=FALSE]),
                         weights=pl$wend, family=family,
                         start=rep.int(0, sum(independent)))
          mindfit.summary <- summary(mindfit)
          theta.ind[independent] <- mindfit$coef
          theta1 <- list(coef=mindfit$coef, 
                         theta=theta.ind,
                         independent=independent)
        } else {
          # unknown warning, just report it
          warning(glm.result$warnings)
          mindfit <- glm.result$value
          mindfit.summary <- summary(mindfit)
          theta.ind[independent] <- mindfit$coef
          theta1 <- list(coef=mindfit$coef, 
                         theta=theta.ind,
                         independent=independent)
        }
      } else {
        # no errors or warnings
        mindfit <- glm.result$value
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
  if(ncol(real.cov)==1){real.cov <- as.vector(real.cov)}
  theta[!m$etamap$offsettheta] <- real.coef
# theta[is.na(theta)] <- 0
  names(theta) <- param_names(m,canonical=TRUE)

#
# Old end
#
  gradient <- rep(NA, length(theta))

  # FIXME: Actually, if case-control sampling was used, this should be positive.
  est.cov <- matrix(0, length(theta),length(theta))
  
  if(length(theta)==1){
   covar <- array(0,dim=c(1,1))
   hess <- array(0,dim=c(1,1))
  }else{
   covar <- diag(rep(0,length(theta)))
   hess <- diag(rep(0,length(theta)))
  }
# covar <- as.matrix(covar[!m$etamap$offsettheta,!m$etamap$offsettheta])
# covar[!is.na(real.coef),!is.na(real.coef)] <- real.cov
  covar[!is.na(theta)&!m$etamap$offsettheta,
        !is.na(theta)&!m$etamap$offsettheta] <- real.cov
  hess[!is.na(theta)&!m$etamap$offsettheta,
        !is.na(theta)&!m$etamap$offsettheta] <- -ginv(real.cov)
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
  message("Finished MPLE.")
  # Output results as ergm-class object
  structure(list(coef=theta,
      iterations=iteration, 
      MCMCtheta=theta, gradient=gradient,
      hessian=hess, covar=covar, failure=FALSE,
      est.cov=est.cov, glm = glm, glm.null = glm.null,
      theta1=theta1),
     class="ergm")
}

