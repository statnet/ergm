approx.hotelling.diff.test<-function(x,y=NULL,mu0=NULL){
  # Note that for we want to get the effective sample size before we
  # convert it to a matrix, in case it's an mcmc.list object.
  x.n <- effectiveSize(x)
  x <- as.matrix(x)
  d <-  colMeans(x)
  if(!is.null(y)){
    # y, if it's given, is the constrained sample, so it's OK if it
    # doesn't vary. (E.g, the extreme case --- completely observed
    # network --- is just one configuration of statistics.)
    # Thus, the effective sample size for nonvarying is set to 1.
    y.n <- pmax(effectiveSize(y),1)
    y <- as.matrix(y)
    d <- d - colMeans(y)
  }

  if(!is.null(mu0)) d <- d - mu0

  method <- paste("Chi-squared approximation to the Hotelling's",
                  if(is.null(y)) "One" else "Two",
                  "Sample T^2-test", "with correction for autocorrelation")
  
  # If a statistic doesn't vary and doesn't match, return a 0 p-value:
  if(any(d[x.n==0]!=0)){
    chi2 <- +Inf
    names(chi2) <- "X-squared"
    df <- ncol(x)
    names(df) <- "df"
    nullval <- if(is.null(mu0)) rep(0, ncol(x)) else mu0
    names(nullval) <- colnames(x)
    
    out <- list(statistic=chi2, parameter=df, p.value=0,
                method = method,
                null.value=nullval,
                alternative="two.sided",
                estimate = d)
    class(out)<-"htest"
    return(out)
  }
  
  # If it doesn't vary and matches, ignore it.
  d <- d[x.n!=0]

  # Remove from y first, since we are changing x.n below.
  if(!is.null(y)){
    y <- y[,x.n!=0,drop=FALSE]
    y.n <- y.n[x.n!=0]
  }
  
  x <- x[,x.n!=0,drop=FALSE]
  x.n <- x.n[x.n!=0]

  
  v <- t(cov(x)/sqrt(x.n))/sqrt(x.n)
  if(!is.null(y)) v <- v + t(cov(y)/sqrt(y.n))/sqrt(y.n)
  chi2 <- t(d)%*%robust.inverse(v)%*%d
  names(chi2) <- "X-squared"
  df <- ncol(x)
  names(df) <- "df"
  nullval <- if(is.null(mu0)) rep(0, ncol(x)) else mu0
  names(nullval) <- colnames(x)
  out <- list(statistic=chi2, parameter=df, p.value=pchisq(chi2,df,lower.tail=FALSE),
              method = method,
              null.value=nullval,
              alternative="two.sided",
              estimate = d,
              covariance = v)
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
