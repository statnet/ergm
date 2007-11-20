ergm.mple<-function(Clist, Clist.miss, m, theta.offset=NULL,
                    MPLEtype="glm", family="binomial",
#                    largestdegree=TRUE,
                    MPLEsamplesize=50000,
                    save.glm=TRUE,
                    maxNumDyadTypes=100000,
                    theta1=NULL, verbose=FALSE, ...)
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
    theta.offset[m$etamap$offsettheta] <- -1
    foffset <- xmat %*% theta.offset !=0
    theta.offset[m$etamap$offsettheta] <- -Inf
   }else{
    foffset <- xmat %*% theta.offset
   }
   xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE]
   colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
   xmat <- xmat[!foffset,,drop=FALSE]
   zy <- zy[!foffset]
   wend <- wend[!foffset]
   foffset <- foffset[!foffset]
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
  }

  if(MPLEtype=="penalized"){
   if(verbose) cat("Using penalized MPLE.\n")
   mplefit <- ergm.pen.glm(
                  zy ~ xmat -1 + offset(foffset),
                  data=data.frame(xmat), weights=wend)
#  mple$deviance <- 2 * (mplefit$loglik-mplefit$loglik[1])[-1]
   mplefit$deviance <- -2*mplefit$loglik
   mplefit$cov.unscaled <- mplefit$var
   mplefit.summary <- mplefit
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEtype=="logitreg"){
    mplefit <- model.matrix(terms(zy ~ .-1,data=data.frame(xmat)),
                           data=data.frame(xmat))
    mplefit <- ergm.logitreg(x=mplefit, y=zy, offset=foffset, wt=wend)
    mplefit.summary <- list(cov.unscaled=mplefit$cov.unscaled)
   }else{
    mplefit <- try(
         glm(zy ~ .-1 + offset(foffset), data=data.frame(xmat),
                   weights=wend, family=family),
                   silent = TRUE)
    if (inherits(mplefit, "try-error")) {
      mplefit <- list(coef=theta.offset, deviance=0,
                      cov.unscaled=diag(theta.offset))
      mplefit.summary <- list(cov.unscaled=diag(theta.offset))
    }else{
      mplefit.summary <- summary(mplefit)
    }
   }
#   if(!is.na(largestdegree)){
#    nxmat <- colnames(xmat)
#    xmat <- xmat[,-ncol(xmat),drop=FALSE]
#    colnames(xmat) <- nxmat[-ncol(xmat)]
#   }
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
                    data=data.frame(xmat[,independent,drop=FALSE]),
                    weights=wend, family=family)
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
#  if(!is.na(largestdegree)){
#   if(nrow(real.cov)==length(real.coef)){
#    real.cov <- real.cov[-nrow(real.cov), -nrow(real.cov)]
#   }
#   real.coef <- real.coef[-length(real.coef)]
#  }
  theta[!m$etamap$offsettheta] <- real.coef
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
# covar <- as.matrix(covar[!m$etamap$offsettheta,!m$etamap$offsettheta])
# covar[!is.na(real.coef),!is.na(real.coef)] <- real.cov
  covar[!is.na(theta)&!m$etamap$offsettheta,!is.na(theta)&!m$etamap$offsettheta] <- real.cov
#
  iteration <-  mplefit$iter 
  samplesize <- NA

# mplefit <- call(MPLEtype, zy ~ 1, family=binomial)
#
  if(MPLEtype=="penalized"){
   mplefit.null <- ergm.pen.glm(zy ~ 1, weights=wend)
  }else{
   options(warn=-1)
#  options(warn=2)
   if(MPLEtype=="logitreg"){
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

mk.pseudonet<-function(meanstats,f,y,ntoggles=length(meanstats)+1,covS=length(meanstats)^2*8,verbose=FALSE){
  if(verbose) cat("Constructing a fake network with correct (or close) meanstats:\n")
  oldwarn<-options()$warn
  on.exit(options(warn=oldwarn))
  options(warn=-1)
  
  n<-network.size(y)
  dir<-is.directed(y)
  bi<-is.bipartite(y)
  y<-network.copy(y)
  max.edges<-n*(n-1)/(2-dir)
    
  rbern.net<-function(y,p){
    y<-network.copy(y)
    y<-network.update(y,as.matrix.network.edgelist(as.network(n,density=p,directed=dir,bipartite=bi)),matrix.type="edgelist")
    y
  }
  summ.net<-function(y){
    summary(as.formula(paste("y~",f[3])))
  }
  
  if(verbose) cat("n=",n,", dir=",dir,", max.edges=",max.edges,"\n",sep="")
  
  all.edges<-as.matrix.network.edgelist(rbern.net(y,1))
  all.edges<-all.edges[order(runif(dim(all.edges)[1])),,drop=FALSE]
  
  if(verbose) cat("Estimating the covariance matrix of network statistics... ")

  dens.stats<-t(sapply(1:covS,function(i)
                       summ.net(rbern.net(y,(i-1)/(covS-1)))))
  
  wt<-robust.inverse(cov(dens.stats))
  if(verbose) cat("Finished.\n")
  
  if(verbose){
    cat("Weight matrix:\n")
    print(wt)
  }
  
  start.density<-(which.min(mahalanobis(dens.stats,meanstats,wt,inverted=TRUE))-1)/(covS-1)

  if(verbose) cat("Starting density:",start.density,"\n")
  
  decider<-function(target,cur,prop,wt)
    mahalanobis(cur,target,wt,inverted=TRUE)-mahalanobis(prop,target,wt,inverted=TRUE)>=sqrt(.Machine$double.eps)*runif(1,-1,1)

  y<-rbern.net(y,start.density)
  ms<-summ.net(y)
  y.prop<-network.copy(y)
               
  t<-0
  i<-floor(runif(ntoggles,1,max.edges+1))
  balance<-0
  last.acc<-0
  while(t<max.edges && (t<=n || t/2>=(t-last.acc))){
    if(verbose>=2) {print(summary(y))}
    
    if(all(ms==meanstats)){
      if(verbose) cat("\nNetwork with right statistics found.\n")
      return(y)
    }
    
    nt<-floor(runif(1,1,ntoggles+1))
    
    e<-all.edges[unique(i[1:nt]),,drop=FALSE]

    if(dim(e)[1]==1) y.e<-y[e[1],e[2]]
    else y.e<-y[e]

    if(exp(sum((-1)^(1-y.e))*balance/length(y.e))>runif(1)){
      if(verbose) cat(t,"/",max.edges,": i=",paste(i,collapse=","),"
      bal=",balance," ",sep="")
      if(verbose) cat("nt=",length(y.e),": ",y.e,"->",1-y.e,sep="")
      t<-t+1

      if(length(y.e)==1) y.prop[e[1],e[2]]<-1-y.e
      else y.prop[e]<-1-y.e
      balance<-balance+if(length(y.e)==1) y.prop[e[1],e[2]]-y.e else
      sum(y.prop[e]-y.e)
      ms.prop<-summ.net(y.prop)
      
      if(verbose) cat(":",ms.prop-ms)
      
      if(decider(meanstats,ms,ms.prop,wt)){
        y<-network.copy(y.prop)
        ms<-ms.prop
        last.acc<-t
        if(verbose) cat(": Acc")
      }else{
        if(verbose) cat(": Rej")
        y.prop<-network.copy(y)
        ms.prop<-ms
      }

      if(verbose)cat(": stats=",ms,"\n")
    }
    i<-(i+1:ntoggles)%%max.edges+1
  }
  return(y)
}
