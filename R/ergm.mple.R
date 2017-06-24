#  File R/ergm.mple.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
##########################################################################
# The <ergm.mple> function finds a maximizer to the psuedolikelihood
# function
#
# --PARAMETERS--
#   Clist            : a list of parameters used for fitting and returned
#                      by <ergm.Cprepare>
#   Clist.miss       : the corresponding 'Clist' for the network of missing
#                      edges returned by <ergm.design>            
#   m                : the model, as returned by <ergm.getmodel>
#   init           : a vector a vector of initial theta coefficients
#   theta.offset     : a logical vector specifying which of the model
#                      coefficients are offset, i.e. fixed  
#   MPLEtype         : the method for MPL estimation as "penalized", "glm" or
#                      "logitreg"; default="glm"    
#   family           : the family to use in the R native routine <glm>; only
#                      applicable if "glm" is the 'MPLEtype'; default="binomial"
#   maxMPLEsamplesize: the sample size to use for endogenous sampling in the psuedo-
#                      likelihood computation; default=1e6
#   save.glm         : whether the mple fit and the null mple fit should be
#                      returned (T or F); if false, NULL is returned for both;
#                      default==TRUE
#   control       : a list of MCMC related parameters; recognized components
#                      include:
#         samplesize : the number of networks to sample
#         Clist.miss : see 'Clist.miss' above; some of the code uses this Clist.miss,
#                      some uses the one above, does this seem right?
#   MHproposal       : an MHproposal object, as returned by <ergm.getMHproposal>
#   verbose          : whether this and the C routines should be verbose (T or F);
#                      default=FALSE
#   ...              : additional parameters passed from within; all will be
#                      ignored
#
# --RETURNED--
#   an ergm object as a list containing several items; for details see
#   the return list in the <ergm> function header (<ergm.mple>=!);
#
######################################################################################

ergm.mple<-function(Clist, Clist.miss, m, init=NULL,
                    MPLEtype="glm", family="binomial",
                    maxMPLEsamplesize=1e+6,
                    save.glm=TRUE,
		    control=NULL, MHproposal=NULL,
                    verbose=FALSE,
                    ...) {
  pl <- ergm.pl(Clist=Clist, Clist.miss=Clist.miss, m=m,
                theta.offset=init,
                maxMPLEsamplesize=maxMPLEsamplesize,
		control=control, MHproposal=MHproposal,
                ignore.offset=MPLEtype=="logitreg",
                verbose=verbose)

  if(MPLEtype=="penalized"){
   if(verbose) cat("Using penalized MPLE.\n")
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
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, m=m, wt=pl$wend,
                             start=init, maxit=control$MPLE.maxit)
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
##########TODO How to handle offsets???????
   if(nrow(pl$xmat) > pl$maxMPLEsamplesize){
#
#   fix deviance for sampled data
#
    mplefit$deviance <- ergm.logisticdeviance(theta=mplefit$coef,
     X=model.matrix(terms(pl$zy.full ~ .-1,data=data.frame(pl$xmat.full)),
                           data=data.frame(pl$xmat.full)),
     y=pl$zy.full, offset=pl$foffset.full)
   }
  }
  theta <- init
  real.coef <- mplefit$coef
  real.cov <- mplefit.summary$cov.unscaled
  if(ncol(real.cov)==1){real.cov <- as.vector(real.cov)}
  theta[!m$etamap$offsettheta] <- real.coef
# theta[is.na(theta)] <- 0
  names(theta) <- coef.names.model(m, FALSE)

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
    mplefit.null <- ergm.pen.glm(pl$zy ~ -1 + offset(pl$foffset), weights=pl$wend)
  }else if(MPLEtype=="logitreg"){
    mplefit.null <- ergm.logitreg(x=matrix(0,ncol=1,nrow=length(pl$zy)),
                                  y=pl$zy, offset=pl$foffset, wt=pl$wend)
  }else{
    mplefit.null <- try(glm(pl$zy ~ -1 + offset(pl$foffset), family=family, weights=pl$wend),
                        silent = TRUE)
    if (inherits(mplefit.null, "try-error")) {
      mplefit.null <- list(coef=0, deviance=0, null.deviance=0,
                           cov.unscaled=diag(1))
    }
  }

  if(save.glm){
    glm <- mplefit
    glm.null <- mplefit.null
  }else{
    glm <- list(deviance=mplefit$deviance)
    glm.null <- list(deviance=mplefit.null$deviance)
  }

  # Output results as ergm-class object
  structure(list(coef=theta,
      iterations=iteration, 
      MCMCtheta=theta, gradient=gradient,
      hessian=NULL, covar=covar, failure=FALSE,
      glm = glm, glm.null = glm.null),
     class="ergm")
}

