ergm.mple.ihs<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    MPLEtype="glm", family="binomial",
                    MPLEsamplesize=50000,
                    save.glm=TRUE,
                    maxNumDyadTypes=100000,
                    theta1=NULL, verbose=FALSE, ...)
{
  pl <- ergm.pl.ihs(Clist=Clist, Clist.miss=Clist.miss, m=m,
                    theta.offset=theta.offset,
                    MPLEsamplesize=MPLEsamplesize,
                    maxNumDyadTypes=maxNumDyadTypes,
                    verbose=verbose)

  if(MPLEtype=="penalized"){
   if(verbose) cat("Using penalized MPLE.\n")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend)
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
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, offset=pl$foffset, wt=pl$wend)
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
    mplefit <- try(
          glm(pl$zy ~ .-1 + offset(pl$foffset), data=data.frame(pl$xmat),
                    weights=pl$wend, family=family),
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
#  Note that the term "match" is depreciated.
#
   if(is.null(theta1)){
    independent.terms <- 
       c("edges","match","nodemain","nodefactor","nodematch","absdiff",
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
   if(nrow(pl$xmat) > pl$MPLEsamplesize){
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
  theta[is.na(theta)] <- 0
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

ergm.pl.ihs<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    MPLEsamplesize=50000,
                    maxNumDyadTypes=100000,
                    verbose=FALSE)
{
  if(Clist.miss$nedges>0){
    temp <- matrix(0,ncol=Clist$n,nrow=Clist$n)
    base <- cbind(as.vector(col(temp)), as.vector(row(temp)))
    base <- base[base[, 2] > base[, 1], ]
    if(Clist.miss$dir){
      base <- cbind(base[,c(2,1)],base)
      base <- matrix(t(base),ncol=2,byrow=T)
    }
    ubase <- base[,1] + Clist$n*base[,2]
#   offset <- !is.na(match(ubase, Clist.miss$tails+Clist.miss$heads*Clist$n))
#   Changed by MSH Oct 13. The original (above) seems just wrong!
    offset <- !is.na(match(ubase, Clist.miss$heads+Clist.miss$tails*Clist$n))
    offset <- 1*offset
    numobs <- Clist$ndyads - sum(offset)
  }else{
    offset <- rep(0,Clist$ndyads)
    numobs <- Clist$ndyads
  }
  z <- .C("MPLE_wrapper",
           as.integer(Clist$heads),    as.integer(Clist$tails),
           as.integer(Clist$nedges),   as.integer(Clist$n), 
           as.integer(Clist$dir),     as.integer(Clist$bipartite),
           as.integer(Clist$nterms), 
           as.character(Clist$fnamestring), as.character(Clist$snamestring),
           as.double(Clist$inputs),
           y = integer(maxNumDyadTypes),
           x = double(maxNumDyadTypes*Clist$nparam),
           weightsvector = integer(maxNumDyadTypes),
           as.double(offset), compressedOffset=double(maxNumDyadTypes),
           as.integer(maxNumDyadTypes),
           PACKAGE="ergm")
#
  uvals <- z$weightsvector!=0
  zy <- z$y[uvals]
  wend <- z$weightsvector[uvals]
  xmat <- matrix(z$x, ncol=Clist$nparam, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
##
## Adjust for meanstats
##
#  if(!is.null(Clist$meanstats)){
#    uobs <- (zy * wend) %*% xmat
#    xmat <- sweep(xmat,2,Clist$meanstats/uobs,"*")
##   xmat <- sweep(sweep(xmat,1,wend,"*"),2,Clist$meanstats/uobs,"*")
##   wend <- wend-wend+mean(wend)
#  }
  dmiss <- z$compressedOffset[uvals]
  rm(z,uvals)
#
# Adjust for the offset
#
  if(any(m$etamap$offsettheta)){
   if(is.null(theta.offset)){
    theta.offset <- rep(0, length=Clist$nparam)
    names(theta.offset) <- m$coef.names
    theta.offset[m$etamap$offsettheta] <- -Inf
   }
#  foffset <- xmat %*% theta.offset
#  shouldoffset <- is.infinite(foffset)
   foffset <- xmat[,!m$etamap$offsettheta]%*%theta.offset[!m$etamap$offsettheta]
   shouldoffset <- apply(abs(xmat[,m$etamap$offsettheta])>1e-8,1,any)
   xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE]
   colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
   xmat <- xmat[!shouldoffset,,drop=FALSE]
   zy <- zy[!shouldoffset]
   wend <- wend[!shouldoffset]
   foffset <- foffset[!shouldoffset]
   theta.offset <- theta.offset[!m$etamap$offsettheta]
  }else{
   foffset <- rep(0, length=length(zy))
   theta.offset <- rep(0, length=Clist$nparam)
   if(Clist$nedges>0){
     theta.offset[1] <- log(Clist$nedges/(Clist$ndyads-Clist$nedges))
   }else{
     theta.offset[1] <- log(1/(Clist$ndyads-1))
   }
   names(theta.offset) <- m$coef.names
  }
  
  if(Clist.miss$nedges>0){
    xmat <- xmat[dmiss==0,,drop=FALSE]
    zy <- zy[dmiss==0]
    wend <- wend[dmiss==0]
    foffset <- foffset[dmiss==0]
    colnames(xmat) <- m$coef.names
  }
  
#   Note:  Logistic regression model is fit without an intercept term.
#   If an intercept is desired, the 1-star term should be included in
#   the model by the user.
#  cat("number of dyads is", Clist$ndyads, "num parameters", Clist$nparam,"\n")
  
##
##  Warning: This used to force the largest degree to be a foil
##
#  if(!is.na(largestdegree)){
#   largestdegree <- grep("degree[1-9]",m$coef.names)
#   largestdegree <- max(c(0,largestdegree))
#   if(largestdegree==0){largestdegree <- NA}
#   if(!is.na(largestdegree)){
##   foffset <- foffset - 50*(xmat[,largestdegree]==-1)
#    xmat <- cbind(xmat,(xmat[,largestdegree]==-1))
#   }
#  }

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
   xmat <- xmat[rsample,,drop=FALSE]
   foffset <- foffset[rsample]
  }else{
   xmat.full <- NULL
   zy.full <- NULL
   foffset.full <- NULL
  }
#
  list(xmat=xmat, zy=zy, foffset=foffset, wend=wend, numobs=numobs,
       xmat.full=xmat.full, zy.full=zy.full, foffset.full=foffset.full,
       theta.offset=theta.offset, MPLEsamplesize=MPLEsamplesize)
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
  asycov <- try(robust.inverse(fit$hessian), silent = TRUE)
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

ergm.maple<-function(pl, m,
                    MPLEtype="glm", family="binomial",
                    save.glm=TRUE,
                    theta1=NULL, verbose=FALSE, ...)
{
  if(MPLEtype=="penalized"){
   if(verbose) cat("Using penalized MPLE.\n")
   mplefit <- ergm.pen.glm(
                  pl$zy ~ pl$xmat -1 + offset(pl$foffset),
                  data=data.frame(pl$xmat), weights=pl$wend)
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
    mplefit <- ergm.logitreg(x=mplefit, y=pl$zy, offset=pl$foffset, wt=pl$wend)
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
    mplefit <- try(
          glm(pl$zy ~ .-1 + offset(pl$foffset), data=data.frame(pl$xmat),
                    weights=pl$wend, family=family),
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
#  Note that the term "match" is depreciated.
#
   if(is.null(theta1)){
    independent.terms <- 
       c("edges","match","nodemain","nodefactor","nodematch","absdiff",
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
   if(nrow(pl$xmat) > pl$MPLEsamplesize){
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
  theta[is.na(theta)] <- 0
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

