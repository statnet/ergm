#  File R/ergm.mple.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' Find a maximizer to the psuedolikelihood function
#' 
#' The \code{ergm.mple} function finds a maximizer to the psuedolikelihood
#' function (MPLE). It is the default method for finding the ERGM starting
#' coefficient values. It is normally called internally the ergm process and
#' not directly by the user. Generally [ergmMPLE()] would be called
#' by users instead.
#' 
#' According to \insertCite{HuHa08e;textual}{ergm}: "The maximizer of the pseudolikelihood
#' may thus easily be found (at least in principle) by using logistic
#' regression as a computational device." In order for this to work, the
#' predictors of the logistic regression model must be calculated.  These are
#' the change statistics as described in Section 3.2 of \insertCite{HuHa08e;textual}{ergm},
#' put into matrix form so that each pair of nodes is one row whose values are
#' the vector of change statistics for that node pair.  The ergm.pl function
#' computes these change statistics and the ergm.mple function implements the
#' logistic regression using \R's [glm()] function.  Generally, neither ergm.mple
#' nor ergm.pl should be called by users if the logistic regression output is
#' desired; instead, use the [ergmMPLE()] function.
#' 
#' In the case where the ERGM is a dyadic independence model, the MPLE is the
#' same as the MLE.  However, in general this is not the case and, as \insertCite{DuGi09f;textual}{ergm}
#' warn, the statistical properties of MPLEs in general are
#' somewhat mysterious.
#' 
#' MPLE values are used even in the case of dyadic dependence models as
#' starting points for the MCMC algorithm.
#' 
#' @param state,state.obs [`ergm_state`] objects.
#' @param init a vector of initial theta coefficients
#' @param family the family to use in the R native routine
#'   [glm()]; only applicable if "glm" is the 'MPLEtype';
#'   default="binomial"
#'
#' @templateVar mycontrol control.ergm
#' @template control
#' @template verbose
#'
#' @param \dots additional parameters passed from within; all will be
#'   ignored
#' @return \code{ergm.mple} returns an ergm object as a list
#'   containing several items; for details see the return list of
#'   [ergm()]
#' 
#' @seealso [ergmMPLE()],
#' [ergm()],[control.ergm()]
#' @references \insertAllCited{}
ergm.mple<-function(s, s.obs, init=NULL,
                    family="binomial",
                    control=NULL,
                    verbose=FALSE,
                    ...) {
  m <- s$model
  message("Starting maximum pseudolikelihood estimation (MPLE):")
  message("Obtaining the responsible dyads.")
  message("Evaluating the predictor and response matrix.")
  pl <- ergm.pl(s, s.obs,
                theta.offset=init,
		control=control,
                ignore.offset=control$MPLE.type=="logitreg",
                verbose=verbose)

  # test whether the MPLE actually exists
  # FIXME: Figure out how to test for MPLE's existence in penalised and curved MPLEs.
  if(control$MPLE.check && ! control$MPLE.type%in%c("penalized","logitreg"))  mple.existence(pl)

  message("Maximizing the pseudolikelihood.")
  if(control$MPLE.type=="penalized"){
   if(verbose) message("Using penalized MPLE.")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend,
                  start=init[!m$etamap$offsettheta])
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
   if(control$MPLE.type=="logitreg"){
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
     glm.result <- quietly(function() glm(pl$zy ~ .-1 + offset(pl$foffset),
                                          data=data.frame(pl$xmat), weights=pl$wend, family=family))()

   # estimate variability matrix V for Godambe covariance matrix or via bootstrapping, only for dyad dependent models and
   #  init.method="MPLE"
     if(!is.dyad.independent(s$model) && control$MPLE.covariance.method=="Godambe" ||
        control$MPLE.covariance.method=="bootstrap"){
       invHess <- summary(glm.result$result)$cov.unscaled
       mple.cov <- ergm_mplecov(pl=pl, s=s, init=init, theta.mple=coef(glm.result$result), invHess=invHess,
                                verbose=verbose, control=control)
     }
    
    # error handling for glm results
     if (length(glm.result$warnings)) {
       glm.result$warnings <- setdiff(glm.result$warnings, "glm.fit: fitted probabilities numerically 0 or 1 occurred")
       # if the glm results are crazy, redo it with 0 starting values
       if (max(abs(coef(glm.result$result)), na.rm=T) > 1e6) {
         warning("GLM may be separable; restarting glm with zeros.\n")
         glm.result <- list(
           result = glm(pl$zy ~ .-1 + offset(pl$foffset),
                       data=data.frame(pl$xmat),
                       weights=pl$wend, family=family, 
                       start=rep.int(0, length(init[!m$etamap$offsettheta])))
         )
       } else if(length(glm.result$warnings)) {
         # unknown warning, just report it
         for(w in glm.result$warnings) warning(w)
       }
     }

     mplefit <- glm.result$result
     mplefit.summary <- summary(mplefit)
   }
  }
  real.coef <- coef(mplefit)
  if(!is.dyad.independent(s$model) && control$MPLE.covariance.method=="Godambe" ||
     control$MPLE.covariance.method=="bootstrap" ){
    real.cov <- mple.cov
  }else{
    real.cov <- mplefit.summary$cov.unscaled
  }

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
   covar <- diag(rep_along(theta, 0))
   hess <- diag(rep_along(theta, 0))
  }

  covar %[.,.]% (!is.NA(theta) & !m$etamap$offsettheta) <- real.cov
  hess %[.,.]% (!is.NA(theta) & !m$etamap$offsettheta) <-
    EVL3(real.cov, -sginv(., tol = .Machine$double.eps^(3 / 4)), 0)

  iteration <-  mplefit$iter 

# mplefit <- call(control$MPLE.type, pl$zy ~ 1, family=binomial)
#
  if(control$MPLE.type=="penalized"){
    mplefit.null <- ergm.pen.glm(pl$zy ~ -1 + offset(pl$foffset), weights=pl$wend)
  }else if(control$MPLE.type=="logitreg"){
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

  nonident <- check_nonidentifiability(pl$xmat.full, theta, m,
                              tol = control$MPLE.nonident.tol, type="covariates",
                              nonident_action = control$MPLE.nonident,
                              nonvar_action = control$MPLE.nonvar)

  # Output results as ergm-class object
  structure(list(coefficients=theta,
      iterations=iteration, 
      MCMCtheta=theta, gradient=gradient,
      hessian=hess, covar=covar, failure=FALSE,
      mple.lik = structure(
        ERRVL2(logLik(mplefit), -mplefit$deviance/2),
        nobs = nobs, df = df, class="logLik"),
      mple.lik.null = structure(
        ERRVL2(logLik(mplefit.null), -mplefit.null$deviance/2),
        nobs = nobs, df = df, class="logLik"),
      lindep = nonident$lindep
      ),
      class="ergm")
}

#' Test whether the MPLE exists
#'
#' The \code{mple.existence} function tests whether the MPLE actually exists. The code
#' applies the approach introduced by \insertCite{Ko07l;textual}{ergm}.
#'
#' \insertNoCite{Ko07l}{ergm}
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
#' @references \insertAllCited{}
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
  for(k in seq_along(c(obj))){
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
