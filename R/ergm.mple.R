#  File R/ergm.mple.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

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
#' @param nw response [`network`] or [`ergm_state`].
#' @param fd An \code{\link{rlebdm}} with informative dyads.
#' @param m the model, as returned by \code{\link{ergm_model}}
#' @param init a vector a vector of initial theta coefficients
#' @param MPLEtype the method for MPL estimation as "penalized", "glm"
#'   or "logitreg"; default="glm"
#' @param family the family to use in the R native routine
#'   \code{\link{glm}}; only applicable if "glm" is the 'MPLEtype';
#'   default="binomial"
#'
#' @templateVar mycontrol control.ergm
#' @template control
#' @template verbose
#'
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
#' *24*(3), pp. 1-29. \doi{10.18637/jss.v024.i03}
#' 
#' van Duijn MAJ, Gile K, Handcock MS (2009).  "Comparison of Maximum Pseudo
#' Likelihood and Maximum Likelihood Estimation of Exponential Family Random
#' Graph Models." _Social Networks_, *31*, pp. 52-62.
ergm.mple<-function(nw, fd, m, init=NULL,
                    MPLEtype="glm", family="binomial",
                    save.xmat=TRUE,
                    control=NULL,
                    verbose=FALSE,
                    ...) {
  message("Starting maximum pseudolikelihood estimation (MPLE):")
  message("Evaluating the predictor and response matrix.")
  pl <- ergm.pl(nw=nw, fd=fd, m=m,
                theta.offset=init,
		control=control,
                ignore.offset=MPLEtype=="logitreg",
                verbose=verbose)

  # test whether the MPLE actually exists
  # FIXME: Figure out how to test for MPLE's existence in penalised and curved MPLEs.
  if(! MPLEtype%in%c("penalized","logitreg"))  mple.existence(pl)

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
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, m=m, wt=pl$wend,
                             start=init, maxit=control$MPLE.maxit, verbose=max(verbose-2,0))
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
      stop(glm.result$error)
    } else if (!is.null(glm.result$warnings)) {
      # if the glm results are crazy, redo it with 0 starting values
      if (max(abs(coef(glm.result$value)), na.rm=T) > 1e6) {
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
  }
  real.coef <- coef(mplefit)
  real.cov <- mplefit.summary$cov.unscaled

  theta <- NVL(init, real.coef)
  theta[!m$etamap$offsettheta] <- real.coef
  names(theta) <- param_names(m, FALSE)

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
        !is.na(theta)&!m$etamap$offsettheta] <- if(length(real.cov)) -sginv(real.cov) else matrix(0,0,0)
#
  iteration <-  mplefit$iter 

# mplefit <- call(MPLEtype, pl$zy ~ 1, family=binomial)
#
  if(MPLEtype=="penalized"){
    mplefit.null <- ergm.pen.glm(pl$zy ~ -1 + offset(pl$foffset), weights=pl$wend)
  }else if(MPLEtype=="logitreg"){
    mplefit.null <- ergm.logitreg(x=matrix(0,ncol=1,nrow=length(pl$zy)),
                                  y=pl$zy, offset=pl$foffset, wt=pl$wend, verbose=max(verbose-2,0))
  }else{
    mplefit.null <- try(glm(pl$zy ~ -1 + offset(pl$foffset), family=family, weights=pl$wend),
                        silent = TRUE)
    if (inherits(mplefit.null, "try-error")) {
      mplefit.null <- list(coefficients=0, deviance=0, null.deviance=0,
                           cov.unscaled=diag(1))
    }
  }

  nobs <- sum(pl$wend)
  df <- nparam(m, offset=FALSE)

  message("Finished MPLE.")
  # Output results as ergm-class object
  structure(list(coefficients=theta,
      iterations=iteration, 
      MCMCtheta=theta, gradient=gradient,
      hessian=hess, covar=covar, failure=FALSE,
      mple.lik = structure(
        ERRVL(try(logLik(mplefit), silent=TRUE), -mplefit$deviance/2),
        nobs = nobs, df = df, class="logLik"),
      mple.lik.null = structure(
        ERRVL(try(logLik(mplefit.null), silent=TRUE), -mplefit.null$deviance/2),
        nobs = nobs, df = df, class="logLik"),
      xmat.full = if(save.xmat) pl$xmat.full
      ),
      class="ergm")
}

#' Test whether the MPLE exists
#'
#' The \code{mple.existence} function tests whether the MPLE actually exists. The code
#' applies the approach introduced by Konis (2007).
#'
#' Konis shows that the MPLE doesn't exist if data may be separated in the sense that
#' there exists a vector beta such that
#'    beta > (T(A^+_{ij})-T(A^-_{ij})) <0  when Aij= 0, and
#'    beta > (T(A^+_{ij})-T(A^-_{ij})) >0  when Aij= 1.
#' Here T(A^+_{ij})-T(A^-_{ij}) is the change
#' statistic of an adjacency matrix A. He derives that finding such beta
#' can be posed as a linear programming problem. In particular,
#' maximize (e' X)beta,
#' subject to X beta >= 0   (1)
#' where e is a vector of ones, and X is the design matrix (T(A^+_ij)-T(A^-_ij)),
#' where each element in a row that corresponds to a dyad with no tie, i.e.,Aij= 0,
#' is being multiplied by -1. If there exist a beta such that (1) has a solution,
#' then the data is separable and the MPLE does not exist.
#'
#' @param pl An ergm.pl-object
#'
#' @references Konis K (2007).  "Linear Programming Algorithms for Detecting Separated
#' Data in Binary LogisticRegression Models (Ph.D. Thesis)." _Worcester College, Oxford University_.
#' \url{https://ora.ox.ac.uk/objects/uuid:8f9ee0d0-d78e-4101-9ab4-f9cbceed2a2a}
#' @noRd
mple.existence <- function(pl) {
#' @importFrom lpSolveAPI make.lp set.column set.objfn set.constr.type set.rhs set.bounds lp.control
  X <- pl$xmat
  y <- pl$zy
  y[y==0] <- -1
  X.bar <- y*X
  e_n <- rep(1, nrow(X.bar))
  obj <- e_n%*%X.bar 
  lprec <- make.lp(nrow=nrow(X.bar), ncol=length(obj)) # set constraint and decision variables
  for(k in 1:length(c(obj))){
    status <- set.column(lprec, k, X.bar[,k])
  }
  status <- set.objfn(lprec, c( obj) )
  status <- set.constr.type(lprec, rep(">=", NROW(X.bar)))
  status <- set.rhs(lprec,  rep(0, NROW(X.bar)))
  status <- set.bounds(lprec, lower = rep(-Inf, length(obj)), upper = rep(Inf, length(obj)))
  control <- lp.control(lprec, pivoting = "firstindex", sense = "max",
                        simplextype = c("primal", "primal"))
  status <- solve(lprec)
  if(status == 3){
    warning("The MPLE does not exist!")
  }
}
