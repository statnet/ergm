.dtsq <- function(x, param, df, log = FALSE){
  fx <- x*(df - param + 1)/(param*df)
  p <- df(fx, param, df - param + 1, log=log)
  if(log) p + log((df - param + 1)/(param*df)) else p*((df - param + 1)/(param*df))
}

.ptsq <- function (q, param, df, lower.tail = TRUE, log.p = FALSE){
  fq <- q*(df - param + 1)/(param*df)
  pf(fq, param, df - param + 1, lower.tail=lower.tail, log.p=log.p)
}

.qtsq <- function(p, param, df, lower.tail = TRUE, log.p = FALSE){
  fq <- qf(fq, param, df - param + 1, lower.tail=lower.tail, log.p=log.p)
  fq / ((df - param + 1)/(param*df))
}

approx.hotelling.diff.test<-function(x,y=NULL, mu0=NULL, assume.indep=FALSE, var.equal=FALSE){

  tr <- function(x) sum(diag(as.matrix(x)))

  vars <- list(x=list(v=x))
  if(!is.null(y)) vars$y <- list(v=y)
  
  mywithin <- function(data, ...) within(data, ...) # This is a workaround suggsted by Duncan Murdoch: calling lapply(X, within, {CODE}) would leave CODE unable to see any objects in f.
  vars <- lapply(vars, mywithin, {
    if(!is.mcmc.list(v))
      v <- mcmc.list(mcmc(as.matrix(v)))
    vcovs.indep <- lapply(v, cov)
    if(assume.indep){
      vcovs <- vcovs.indep
    }else{
      vcovs <- lapply(lapply(v, .ergm.mvar.spec0), function(m) matrix(ifelse(is.na(c(m)), 0, c(m)),nrow(m),ncol(m)))
    }
    ms <- lapply(v, base::colMeans)
    m <- colMeans(as.matrix(v))
    ns <- sapply(v,base::nrow)
    n <- sum(ns)

    # These are pooled estimates of the variance-covariance
    # matrix. Note that the outer product of the difference between
    # chain means (times n) is added on as well, because the chains
    # are supposed to all have the same population mean. However, the
    # divisor is then the combined sample size less 1, because we are
    # assuming equal means.
    vcov.indep <- Reduce("+", Map("+", Map("*", vcovs.indep, ns-1), Map("*", Map(outer, lapply(ms,"-",m), lapply(ms,"-",m)), ns) ))/(n-1)
    vcov <- Reduce("+", Map("+", Map("*", vcovs, ns-1), Map("*", Map(outer, lapply(ms,"-",m), lapply(ms,"-",m)), ns) ))/(n-1)

    infl <- tr(vcov) / tr(vcov.indep) # I.e., how much bigger is the trace of the variance-covariance after taking autocorrelation into account than before.
    neff <- n / infl
    
    vcov.m <- vcov/n # Here, vcov already incorporates the inflation due to autocorrelation.
    v <- as.matrix(v)
  })
  rm(mywithin)
  
  x <- vars$x
  y <- vars$y
  
  d <- x$m
  vcov.d <- x$vcov.m
  
  if(!is.null(y)){
    d <- d - y$m
    if(var.equal){
      vcov.pooled <- (x$vcov*(x$n-1) + y$vcov*(y$n-1))/(x$n+y$n-2)
      vcov.d <- vcov.pooled * (1/x$n + 1/y$n)
    }else{
      vcov.d <- vcov.d + y$vcov.m
    }
  }


  p <- ncol(x$v)
  if(is.null(mu0)) mu0 <- rep(0,p)
  names(mu0)<-colnames(x$v)

  novar <- diag(vcov.d)==0
  p <- p-sum(novar)

  ivcov.d <-ginv(vcov.d[!novar,!novar,drop=FALSE])
  
  method <- paste("Hotelling's",
                  if(is.null(y)) "One" else "Two",
                  "-Sample",if(var.equal) "Pooled","T^2-Test", if(!assume.indep) "with correction for autocorrelation")
  
  # If a statistic doesn't vary and doesn't match, return a 0 p-value:
  if(any((d-mu0)[novar]!=0)){
    warning("Vector(s) ", paste.and(colnames(x)[novar]),
            if(is.null(y)) " do not vary in x or in y and have differences unequal to mu0"
            else " do not vary and do not equal mu0",
            "; P-value has been set to 0.")
        
    T2 <- +Inf
  }else{
    if(any(novar)){
      warning("Vector(s) ", paste.and(colnames(x)[novar]),
              if(is.null(y)) " do not vary in x or in y but have differences equal to mu0"
              else " do not vary but equal mu0",
              "; they have been ignored for the purposes of testing.")
    }
    T2 <- c(t((d-mu0)[!novar])%*%ivcov.d%*%(d-mu0)[!novar])
  }

  NANVL <- function(z, ifNAN) ifelse(is.nan(z),ifNAN,z)
  
  names(T2) <- "T^2"
  pars <- c(param = p, df = if(is.null(y)){
    NANVL(x$neff,1)-1
  }else if(var.equal){
    NANVL(x$neff,1)+NANVL(y$neff,1)-2
  }else{
    mywith <- function(data, ...) with(data, ...)
    # This is the Krishnamoorthy and Yu (2004) degrees of freedom formula, courtesy of Wikipedia.
    df <- (p+p^2)/sum(NANVL(sapply(vars, mywith, (tr(vcov.m[!novar,!novar] %*% ivcov.d %*% vcov.m[!novar,!novar] %*% ivcov.d) +
                                            tr(vcov.m[!novar,!novar] %*% ivcov.d)^2)/neff), 0))
    rm(mywith)
    df
  })

  if(pars[1]>=pars[2]) warning("Effective degrees of freedom (",pars[2],") must exceed the number of varying parameters (",pars[1],"). P-value will not be computed.")
  out <- list(statistic=T2, parameter=pars, p.value=if(pars[1]<pars[2]) .ptsq(T2,pars[1],pars[2],lower.tail=FALSE) else NA,
              method = method,
              null.value=mu0,
              alternative="two.sided",
              estimate = d,
              covariance = vcov.d)
  class(out)<-"htest"
  out
}

## The following function's bookkeeping parts (e.g., handling of
## mcmc.list and calculation of windows starts and ends) are loosely
## based on parts of geweke.diag() from the coda R package.
##
## Rather than comparing each mean independently, compares them
## jointly. Note that it returns an htest object, not a geweke.diag
## object.
##
## If approx.hotelling.diff.test returns an error, then assume that
## burn-in is insufficient.
.geweke.diag.mv <- function(x, frac1 = 0.1, frac2 = 0.5){
  if(inherits(x, "mcmc.list"))
    return(lapply(x, .geweke.diag.mv, frac1, frac2))
  x <- as.mcmc(x)
  x.len <- end(x) - start(x)
  x1 <- window(x, start=start(x), end=start(x) + frac1*x.len)
  x2 <- window(x, start=end(x) - frac2*x.len, end=end(x))

  test <- approx.hotelling.diff.test(x1,x2,var.equal=TRUE) # When converged, the chain should have the same variance throughout.
  if(is.na(test$p.value)) test$p.value <- 0 # Interpret too-small a sample size as insufficient burn-in.

  test$method <- paste("Multivariate extension to Geweke's burn-in convergence diagnostic")
  test
}

# Compute the sample estimating equations of an ERGM. If the model is
# linear, all non-offset statistics are passed. If the model is
# curved, the score estimating equations (3.1) by Hunter and
# Handcock (2006) are given instead.
.ergm.esteq <- function(theta, model, statsmatrix){
  etamap <- NVL(model$etamap, model)
  esteq <- t(ergm.etagradmult(theta,t(as.matrix(statsmatrix)),etamap))[,!etamap$offsettheta,drop=FALSE]
  if(is.mcmc(statsmatrix)){
    esteq <- mcmc(esteq, start=start(statsmatrix), end=end(statsmatrix), thin=thin(statsmatrix))
    varnames(esteq) <- NVL(names(theta), .coef.names.model(model, FALSE))[!etamap$offsettheta]
  }else  colnames(esteq) <- NVL(names(theta), .coef.names.model(model, FALSE))[!etamap$offsettheta]

  esteq
}

# This function can be viewed as a multivariate version of coda's
# spectrum0.ar().  Its return value, divided by nrow(cbind(x)), is the
# estimated variance-covariance matrix of the sampling distribution of
# the mean of x if x is a multivatriate time series with AR(p)
# structure, with p determined by AIC.
#
# ar() fails if crossprod(x) is singular, which is remedied by mapping
# the variables onto the principal components of x, dropping redundant
# dementions.
#
# FIXME: Actually, for MCMC with multiple chains, we should be using the pooled mean.
.ergm.mvar.spec0 <- function(x, order.max=NULL, aic=is.null(order.max), tol=.Machine$double.eps^0.5, ...){
    x <- cbind(x)
    n <- nrow(x)
    p <- ncol(x)

    v <- matrix(NA,p,p)
    novar <- abs(apply(x,2,stats::sd))<tol
    x <- x[,!novar,drop=FALSE]

    if(ncol(x)){
      # Map the variables onto their principal components, dropping
      # redundant (linearly-dependent) dimensions. Here, we keep the
      # eigenvectors such that the reciprocal condition number defined
      # as s.min/s.max, where s.min and s.max are the smallest and the
      # biggest singular values, respectively, is greater than the
      # tolerance.
      e <- eigen(cov(x), symmetric=TRUE)
      Q <- e$vec[,sqrt(pmax(e$val,0)/max(e$val))>tol*2,drop=FALSE]
      xr <- x%*%Q # Columns of xr are guaranteed to be linearly independent.

      # Calculate the time-series variance of the mean on the PC scale.

      if(is.null(order.max)){ord <- 10*log10(nrow(xr))}
      arfit <- .catchToList(ar(xr,aic=is.null(order.max), order.max=ord, ...))
      while(!is.null(arfit$error) & ord > 1){
        ord <- ord - 1
        arfit <- .catchToList(ar(xr,aic=is.null(order.max), order.max=ord, ...))
      }
      arfit <- arfit$value
      arvar <- arfit$var.pred
      arcoefs <- arfit$ar
      arcoefs <- if(is.null(dim(arcoefs))) sum(arcoefs) else apply(arcoefs,2:3,base::sum)

      adj <- diag(1,nrow=ncol(xr)) - arcoefs
      iadj <- solve(adj)
      v.var <- iadj %*% arvar %*% t(iadj)

      # Reverse the mapping for the variance estimate.
      v.var <- Q%*%v.var%*%t(Q)
    
      v[!novar,!novar] <- v.var
    }
    v
}

