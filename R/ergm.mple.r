ergm.mple<-function(Clist, m, fix=NULL, theta.offset=NULL,
                    MPLEonly="glm", family="binomial",
                    largestdegree=TRUE, MPLEsamplesize=50000,
                    save.glm=TRUE,
                    maxNumDyadTypes=100000,
                    theta1=NULL, verbose=FALSE, ...)
{
  numobs <- Clist$ndyads
  offset <- rep(0,Clist$ndyads)
  z <- .C("MPLE_wrapper",
           as.double(Clist$heads),    as.double(Clist$tails),
           as.double(Clist$nedges),   as.double(Clist$n), 
           as.integer(Clist$dir),     as.double(Clist$bipartite),
           as.integer(Clist$nterms), 
           as.character(Clist$fnamestring), as.character(Clist$snamestring),
           as.double(Clist$inputs),
           y = integer(maxNumDyadTypes),
           x = double(maxNumDyadTypes*Clist$nparam),
           weightsvector = integer(maxNumDyadTypes),
           as.double(offset), compressedOffset=double(maxNumDyadTypes),
           as.double(maxNumDyadTypes),
           PACKAGE="statnet")
#
  uvals <- z$weightsvector!=0
  zy <- z$y[uvals]
  wend <- z$weightsvector[uvals]
  xmat <- as.matrix(matrix(z$x, ncol=Clist$nparam, byrow=TRUE)[uvals,])
##
## Adjust for meanstats
##
#  if(!is.null(Clist$meanstats)){
#    uobs <- (zy * wend) %*% xmat
#    xmat <- sweep(xmat,2,Clist$meanstats/uobs,"*")
##   xmat <- sweep(sweep(xmat,1,wend,"*"),2,Clist$meanstats/uobs,"*")
##   wend <- wend-wend+mean(wend)
#  }
  rm(z,uvals)
#
# Adjust for the offset
#
  if(is.null(fix) || all(!fix) ){
   fix <- rep(FALSE, length=Clist$nparam)
   theta.offset <- rep(0, length=Clist$nparam)
   names(theta.offset) <- m$coef.names
   foffset <- rep(0, length=nrow(xmat))
   colnames(xmat) <- m$coef.names
  }else{
   foffset <- xmat %*% theta.offset
   xmat <- as.matrix(xmat[,!fix])
   colnames(xmat) <- m$coef.names[!fix]
  }
  
  
#   Note:  Logistic regression model is fit without an intercept term.
#   If an intercept is desired, the 1-star term should be included in
#   the model by the user.
#  cat("number of dyads is", Clist$ndyads, "num parameters", Clist$nparam,"\n")
  
#
#  Warning: This used to force the largest degree to be a foil
#
  if(!is.na(largestdegree)){
   largestdegree <- grep("degree[1-9]",m$coef.names)
   largestdegree <- max(c(0,largestdegree))
   if(largestdegree==0){largestdegree <- NA}
   if(!is.na(largestdegree)){
#   foffset <- foffset - 50*(xmat[,largestdegree]==-1)
    xmat <- cbind(xmat,(xmat[,largestdegree]==-1))
   }
  }

#
# Sample if necessary
#
  if(nrow(xmat) > MPLEsamplesize){
   rsample <- sample((1:nrow(xmat))[zy==1], size=min(MPLEsamplesize,sum(zy)),
                     replace=FALSE)
   rsample <- c(rsample, 
     sample((1:nrow(xmat))[zy==0], size=min(MPLEsamplesize,sum(!zy)),
                     replace=FALSE) )
   tau <- sum(zy*wend)/sum(wend)
   xmat.full <- xmat
   zy.full <- zy
   foffset.full <- foffset
   zy <- zy[rsample]
   wend <- wend[rsample]
   wend <- tau*zy*wend/sum(zy*wend) +
           (1-tau)*(1-zy)*wend/sum((1-zy)*wend)
   wend <- wend*nrow(xmat)/sum(wend)
   xmat <- as.matrix(xmat[rsample,])
   foffset <- foffset[rsample]
  }

  if(MPLEonly=="penalized"){
   mplefit <- ergm.pen.glm(
                  zy ~  xmat -1 + offset(foffset),
                  data=data.frame(xmat), weights=wend)
#  mple$deviance <- 2 * (mplefit$loglik-mplefit$loglik[1])[-1]
   mplefit$deviance <- -2*mplefit$loglik
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEonly=="logitreg"){
    mplefit <- model.matrix(terms(zy ~ .-1,data=data.frame(xmat)),
                           data=data.frame(xmat))
    mplefit <- ergm.logitreg(x=mplefit, y=zy, offset=foffset, wt=wend)
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
    mplefit <- glm(zy ~ .-1 + offset(foffset), data=data.frame(xmat),
                   weights=wend, family=family, ...)
    mplefit.summary <- summary(mplefit)
   }
   if(!is.na(largestdegree)){
    nxmat <- colnames(xmat)
    xmat <- as.matrix(xmat[,-ncol(xmat)])
    colnames(xmat) <- nxmat[-ncol(xmat)]
   }
#
#  Determine the independence theta and MLE
#  Note that the term "match" is depreciated.
#
   if(is.null(theta1)){
    independent.terms <- 
       c("edges","match","nodemain","nodefactor","nodematch","absdiff",
         "edgecov","dyadcov","sender","receiver","sociality", 
         "nodemix","mix",
         "actor.","event.",
         "testme")
    independent <- rep(0,ncol(xmat))
    names(independent) <- colnames(xmat)
    theta.ind <- independent
    for(i in seq(along=independent.terms)){
     independent[grep(independent.terms[i], colnames(xmat))] <- i
    }
    independent <- independent>0
    if(any(independent)){
     mindfit <- glm(zy ~ .-1 + offset(foffset), 
                    data=data.frame(xmat[,independent]),
                    weights=wend, family=family, ...)
     mindfit.summary <- summary(mindfit)
     theta.ind[independent] <- mindfit$coef
     theta1 <- list(coef=mindfit$coef, 
                    theta=theta.ind,
                    independent=independent,
                    loglikelihood=-mindfit$deviance/2)
    }else{
     theta1 <- list(coef=NULL, 
                    theta=rep(0,ncol(xmat)),
                    independent=independent,
                    loglikelihood=-numobs*log(2))
    }
   }
#
   options(warn=0)
#  options(warn=2)
   if(nrow(xmat) > MPLEsamplesize){
#
#   fix aic and deviance for sampled data
#
    mplefit$deviance <- ergm.logisticdeviance(beta=mplefit$coef,
     X=model.matrix(terms(zy.full ~ .-1,data=data.frame(xmat.full)),
                           data=data.frame(xmat.full)),
     y=zy.full, offset=foffset.full)
    mplefit$aic <- mplefit$deviance + 2*mplefit$rank
   }
  }
  theta <- theta.offset
  real.coef <- mplefit$coef
  real.cov <- mplefit.summary$cov.unscaled
  if(!is.na(largestdegree)){
   if(nrow(real.cov)==length(real.coef)){
    real.cov <- real.cov[-nrow(real.cov), -nrow(real.cov)]
   }
   real.coef <- real.coef[-length(real.coef)]
  }
  theta[!fix] <- real.coef
  theta[is.na(theta)] <- 0

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
  covar <- as.matrix(covar[!fix,!fix])
  covar[!is.na(real.coef),!is.na(real.coef)] <- real.cov
#
  iteration <-  mplefit$iter 
  samplesize <- NA

# mplefit <- call(MPLEonly, zy ~ 1, family=binomial)
#
  if(MPLEonly=="penalized"){
   mplefit.null <- ergm.pen.glm(zy ~ 1, weights=wend)
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEonly=="logitreg"){
    mplefit.null <- ergm.logitreg(x=matrix(1,ncol=1,nrow=length(zy)),
                                  y=zy, offset=foffset, wt=wend)
   }else{
    mplefit.null <- glm(zy ~ 1, family=family, weights=wend)
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
ergm.logitreg <- function(x, y, wt = rep(1, length(y)),
                          intercept = FALSE, start = rep(0, p),
                          offset=NULL, maxit=200, ...)
{
  gmin <- function(beta, X, y, w, offset) {
      eta <- (X %*% beta)+offset; p <- plogis(eta)
      -2 * matrix(w *dlogis(eta) * ifelse(y, 1/p, -1/(1-p)), 1) %*% X
  }
  if(is.null(dim(x))) dim(x) <- c(length(x), 1)
  if(is.null(offset)) offset <- rep(0,length(y))
  dn <- dimnames(x)[[2]]
  if(!length(dn)) dn <- paste("Var", 1:ncol(x), sep="")
  p <- ncol(x) + intercept
  if(intercept) {x <- cbind(1, x); dn <- c("(Intercept)", dn)}
  if(is.factor(y)) y <- (unclass(y) != 1)
  fit <- optim(start, ergm.logisticdeviance, gmin,
               X = x, y = y, w = wt, offset=offset,
               method = "BFGS", hessian=TRUE, control=list(maxit=maxit), ...)
  names(fit$par) <- dn
  fit$coef <- fit$par
  fit$deviance <- fit$value
  fit$iter <- fit$counts[1]
  asycov <- try(solve(fit$hessian), silent = TRUE)
  if (inherits(asycov, "try-error")) {
     asycov <- diag(1/diag(-fit$hessian))
  }
  fit$cov.unscaled <- asycov
# cat("\nCoefficients:\n"); print(fit$par)
# # R: use fit$value and fit$convergence
# cat("\nResidual Deviance:", format(fit$value), "\n")
  if(fit$convergence > 0)
      cat("\nConvergence code:", fit$convergence, "\n")
  invisible(fit)
}

ergm.logisticdeviance <- function(beta, X, y,
                            w=rep(1,length(y)), offset=rep(0,length(y))) {
      p <- plogis((X %*% beta)+offset)
      -sum(2 * w * ifelse(y, log(p), log(1-p)))
}
