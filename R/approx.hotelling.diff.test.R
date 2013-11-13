approx.hotelling.diff.test<-function(x,y=NULL,mu0=NULL){
  if(is.mcmc.list(x)){
    x.spec0 <- lapply(x, .ergm.mvar.spec0)
    x.n <- lapply(x,nrow)
    vcov.xm <- Reduce("+", Map("*", x.spec0, x.n))/Reduce("+", x.n)^2
    x <- as.matrix(x)
  }else{
    x <- as.matrix(x)
    x.spec0 <- .ergm.mvar.spec0(x)
    x.n <- nrow(x)
    vcov.xm <- x.spec0/x.n
  }
  vcov.xm[is.na(vcov.xm)] <- 0

  d <-  colMeans(x)
  vcov.d <- vcov.xm
  
  if(!is.null(y)){
    if(is.mcmc.list(y)){
      y.spec0 <- lapply(y, .ergm.mvar.spec0)
      y.n <- lapply(y,nrow)
      vcov.ym <- Reduce("+", Map("*", y.spec0, y.n))/Reduce("+", y.n)^2
      y <- as.matrix(y)
    }else{
      y <- as.matrix(y)
      y.spec0 <- .ergm.mvar.spec0(y)
      y.n <- nrow(y)
      vcov.ym <- y.spec0/y.n
    }
    vcov.ym[is.na(vcov.ym)] <- 0

    d <- d - colMeans(y)
    vcov.d <- vcov.d + vcov.ym
  }

  novar <- diag(vcov.d)==0
  
  if(is.null(mu0)) mu0 <- rep(0,ncol(x))
  names(mu0)<-colnames(x)

  method <- paste("Chi-squared approximation to the Hotelling's",
                  if(is.null(y)) "One" else "Two",
                  "Sample T^2-test", "with correction for autocorrelation")
  
  # If a statistic doesn't vary and doesn't match, return a 0 p-value:
  if(any((d-mu0)[novar]!=0)){
    warning("Vector(s) ", paste.and(colnames(x)[novar]),
            if(is.null(y)) " do not vary in x or in y and have differences unequal to mu0"
            else " do not vary and do not equal mu0",
            "; P-value has been set to 0.")
        
    chi2 <- +Inf
  }else{
    if(any(novar)){
      warning("Vector(s) ", paste.and(colnames(x)[novar]),
              if(is.null(y)) " do not vary in x or in y but have differences equal to mu0"
              else " do not vary but equal mu0",
              "; they have been ignored for the purposes of testing.")
    }
    chi2 <- t((d-mu0)[!novar])%*%robust.inverse(vcov.d[!novar,!novar,drop=FALSE])%*%(d-mu0)[!novar]
  }
  
  names(chi2) <- "X-squared"
  df <- ncol(x)-sum(novar)
  names(df) <- "df"
  out <- list(statistic=chi2, parameter=df, p.value=pchisq(chi2,df,lower.tail=FALSE),
              method = method,
              null.value=mu0,
              alternative="two.sided",
              estimate = d,
              covariance = vcov.d)
  class(out)<-"htest"
  out
}

## The following function uses small parts of geweke.diag from the
## coda R package. The original code is Copyright (C) 2005-2011 Martyn
## Plummer, Nicky Best, Kate Cowles, Karen Vines
##
## It is incorporated into the ergm package under the terms of the GPL
## v3.
##
## Rather than comparing each mean independently, compares them
## jointly. Note that it returns an htest object, not a geweke.diag
## object.

geweke.diag.mv <- function(x, frac1 = 0.1, frac2 = 0.5){
  if (is.mcmc.list(x)) 
    return(lapply(x, geweke.diag.mv, frac1, frac2))
  x <- as.mcmc(x)
  x1 <- window(x, start=start(x), end=start(x) + frac1 * (end(x) - start(x)))
  x2 <- window(x, start=end(x) - frac2 * (end(x) - start(x)), end=end(x))
  test <- approx.hotelling.diff.test(x1,x2)
  test$method <- paste("Multivariate extension to Geweke's burn-in convergence diagnostic")
  return(test)
}

# Compute the sample estimating equations of an ERGM. If the model is
# linear, all non-offset statistics are passed. If the model is
# curved, the score estimating equations (3.1) by Hunter and
# Handcock (2006) are given instead.
.ergm.esteq <- function(theta, model, statsmatrix){
  esteq <- t(ergm.etagradmult(theta,t(as.matrix(statsmatrix)),model$etamap))[,!model$etamap$offsettheta,drop=FALSE]
  colnames(esteq) <- names(theta)[!model$etamap$offsettheta]
  esteq
}

# This function can be viewed as a multivariate version of coda's
# spectrum0.ar().  Its return value, divided by nrow(cbind(x)), is the
# estimated variance-covariance matrix of the sampling distribution of
# the mean of x if x is a multivatriate time series with AR(p)
# structure, with p determined by AIC.
#
# ar() fails if crossprod(x) is singular, so 
.ergm.mvar.spec0 <- function(x,tol=.Machine$double.eps){
    x <- cbind(x)
    n <- nrow(x)
    p <- ncol(x)

    v <- matrix(NA,p,p)
    novar <- apply(x,2,sd)==0 # FIXME: Add tolerance?
    x <- x[,!novar,drop=FALSE]

    if(ncol(x)){
      arfit <- try(ar(x,aic=TRUE),silent=TRUE)
      if(inherits(arfit,"try-error")){
        warning("Excessive correlation among the statistics. Using a separate effectiveSize approximation.")
        x.factor <- sqrt(n/effectiveSize(x))
        v.var <- t(cov(x)*x.factor)*x.factor
      }else{
        arvar <- arfit$var.pred
        arcoefs <- arfit$ar
        arcoefs <- if(is.null(dim(arcoefs))) sum(arcoefs) else apply(arcoefs,2:3,sum)
        adj <- diag(1,nrow=p-sum(novar)) - arcoefs
        iadj <- solve(adj)
        v.var <- iadj %*% arfit$var.pred %*% t(iadj)
      }
      v[!novar,!novar] <- v.var
    }
    v
}
