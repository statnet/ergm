#  File R/approx.hotelling.diff.test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
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



#' Approximate Hotelling T^2-Test for One Sample Means
#' 
#' A multivariate hypothesis test for a single population mean or a
#' difference between them. This version attempts to adjust for
#' multivariate autocorrelation in the samples.
#'
#' @param x a numeric matrix of data values with cases in rows and
#'   variables in columns.
#' @param y an optinal matrix of data values with cases in rows and
#'   variables in columns for a 2-sample test.
#' @param mu0 an optional numeric vector: for a 1-sample test, the
#'   poulation mean under the null hypothesis; and for a 2-sample
#'   test, the difference between population means under the null
#'   hypothesis; defaults to a vector of 0s.
#' @param assume.indep if `TRUE`, performs an ordinary Hotelling's
#'   test without attempting to account for autocorrelation.
#' @param var.equal for a 2-sample test, perform the pooled test:
#'   assume population variance-covariance matrices of the two
#'   variables are equal.
#'
#' @return An object of class `htest` with the following information:
#' \item{statistic}{The \eqn{T^2} statistic.}
#' \item{parameter}{Degrees of freedom.}
#' \item{p.value}{P-value.}
#' \item{method}{Method specifics.}
#' \item{null.value}{Null hypothesis mean or mean difference.}
#' \item{alternative}{Always `"two.sided"`.}
#' \item{estimate}{Sample difference.}
#' \item{covariance}{Estimated variance-covariance matrix of the estimate of the difference.}
#' 
#' It has a print method [print.htest()].
#'
#' @seealso [t.test()]
#' @note For [`mcmc.list`] input, the variance for this test is
#'   estimated with unpooled means. This is not strictly correct.
#' @references
#' 
#' Hotelling, H. (1947). Multivariate Quality Control. In C. Eisenhart, M. W.
#' Hastay, and W. A. Wallis, eds. Techniques of Statistical Analysis. New York:
#' McGraw-Hill.
#'
#' @export approx.hotelling.diff.test
approx.hotelling.diff.test<-function(x,y=NULL, mu0=0, assume.indep=FALSE, var.equal=FALSE){
  if(!is.mcmc.list(x))
    x <- mcmc.list(mcmc(as.matrix(x)))
  if(!is.null(y) && !is.mcmc.list(y))
    y <- mcmc.list(mcmc(as.matrix(y)))

  if(is.null(mu0)) mu0 <- rep(0,nvar(x))
  mu0 <- rep(mu0, length.out = nvar(x))

  tr <- function(x) sum(diag(as.matrix(x)))

  vars <- list(x=list(v=x))
  if(!is.null(y)) vars$y <- list(v=y)

  v <- NULL # Prevent a spurious R CMD check warning.
  mywithin <- function(data, ...) within(data, ...) # This is a workaround suggsted by Duncan Murdoch: calling lapply(X, within, {CODE}) would leave CODE unable to see any objects in f.
  vars <- lapply(vars, mywithin, {
    vcovs.indep <- lapply(v, cov)
    if(assume.indep){
      vcovs <- vcovs.indep
    }else{
      vcovs <- lapply(v, spectrum0.mvar)
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
  })
  rm(mywithin, v)
  
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


  p <- nvar(x$v)
  names(mu0)<-varnames(x$v)

  novar <- diag(vcov.d)==0
  p <- p-sum(novar)

  ivcov.d <-ginv(vcov.d[!novar,!novar,drop=FALSE])
  
  method <- paste("Hotelling's",
                  NVL2(y, "Two", "One"),
                  "-Sample",if(var.equal) "Pooled","T^2-Test", if(!assume.indep) "with correction for autocorrelation")
  
  # If a statistic doesn't vary and doesn't match, return a 0 p-value:
  if(any((d-mu0)[novar]!=0)){
    warning("Vector(s) ", paste.and(colnames(x)[novar]),
            NVL2(y,
                 " do not vary and do not equal mu0",
                 " do not vary in x or in y and have differences unequal to mu0"),
            "; P-value has been set to 0.")
        
    T2 <- +Inf
  }else{
    if(any(novar)){
      warning("Vector(s) ", paste.and(colnames(x)[novar]),
              NVL2(y,
                   " do not vary but equal mu0",
                   " do not vary in x or in y but have differences equal to mu0"),
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

#' Multivariate version of `coda`'s [coda::geweke.diag()].
#' 
#' Rather than comparing each mean independently, compares them
#' jointly. Note that it returns an `htest` object, not a `geweke.diag`
#' object.
#'
#' @param x an [`mcmc`], [`mcmc.list`], or just a matrix with
#'   observations in rows and variables in columns.
#' @param frac1,frac2 the fraction at the start and, respectively, at
#'   the end of the sample to compare.
#' @param split.mcmc.list when given an `mcmc.list`, whether to test
#'   each chain individually.
#' @note If [approx.hotelling.diff.test()] returns an error, then
#'   assume that burn-in is insufficient.
#' @return An object of class `htest`, inheriting from that returned
#'   by [approx.hotelling.diff.test()], but with p-value considered to
#'   be 0 on insufficient sample size.
#'
#' @seealso [coda::geweke.diag()], [approx.hotelling.diff.test()]
#' @export
geweke.diag.mv <- function(x, frac1 = 0.1, frac2 = 0.5, split.mcmc.list = FALSE){
  # The following function's bookkeeping parts (e.g., handling of
  # mcmc.list and calculation of windows starts and ends) are loosely
  # based on parts of geweke.diag() from the coda R package.
  
  if(is.mcmc.list(x)){
    if(split.mcmc.list){
      return(lapply(x, geweke.diag.mv, frac1, frac2))
    }
  }else{
    x <- as.mcmc(x)
  }
  
  x.len <- end(x) - start(x)
  x1 <- window(x, start=start(x), end=start(x) + frac1*x.len)
  x2 <- window(x, start=end(x) - frac2*x.len, end=end(x))

  test <- approx.hotelling.diff.test(x1,x2,var.equal=TRUE) # When converged, the chain should have the same variance throughout.
  if(is.na(test$p.value)) test$p.value <- 0 # Interpret too-small a sample size as insufficient burn-in.

  test$method <- paste("Multivariate extension to Geweke's burn-in convergence diagnostic")
  test
}

#' Multivariate version of `coda`'s [spectrum0.ar()].
#'
#' Its return value, divided by `nrow(cbind(x))`, is the estimated
#' variance-covariance matrix of the sampling distribution of the mean
#' of `x` if `x` is a multivatriate time series with AR(\eqn{p}) structure, with
#' \eqn{p} determined by AIC.
#'
#' @param x a matrix with observations in rows and variables in
#'   columns.
#' @param order.max maximum (or fixed) order for the AR model.
#' @param aic use AIC to select the order (up to `order.max`).
#' @param tol drop components until the reciprocal condition number of
#'   the transformed variance-covariance matrix is greater than this.
#' @param ... additional arguments to [ar()].
#'
#' @note [ar()] fails if `crossprod(x)` is singular,
#' which is remedied by mapping the variables onto the principal
#' components of `x`, dropping redundant dimentions.
#' @export spectrum0.mvar
spectrum0.mvar <- function(x, order.max=NULL, aic=is.null(order.max), tol=.Machine$double.eps^0.5, ...){
  x <- cbind(x)
  n <- nrow(x)
  p <- ncol(x)
  
  v <- matrix(0,p,p)
  novar <- abs(apply(x,2,stats::sd))<tol
  x <- x[,!novar,drop=FALSE]

  # Index of the first local minimum in a sequence.
  first_local_min <- function(x){
    d <- diff(c(Inf,x,Inf))
    min(which(d>=0))-1
  }
  
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

    ord <- NVL(order.max, ceiling(10*log10(nrow(xr))))
    arfit <- .catchToList(ar(xr,aic=is.null(order.max), order.max=ord, ...))
    # If ar() failed or produced a variance matrix estimate that's
    # not positive semidefinite, try with a lower order.
    while((!is.null(arfit$error) || ERRVL(try(any(eigen(arfit$value$var.pred, only.values=TRUE)$values<0), silent=TRUE), TRUE)) && ord > 1){
      ord <- ord - 1
      arfit <- .catchToList(ar(xr,aic=is.null(order.max), order.max=ord, ...))
    }
    
    arfit <- arfit$value
    if(aic && arfit$order>(ord <- first_local_min(arfit$aic)-1)){
      arfit <- ar(xr, aic=ord==0, order.max=max(ord,1)) # Workaround since ar() won't take order.max=0.
    }
    
    arvar <- arfit$var.pred
    arcoefs <- arfit$ar
    arcoefs <- NVL2(dim(arcoefs), apply(arcoefs,2:3,base::sum), sum(arcoefs))
    
    adj <- diag(1,nrow=ncol(xr)) - arcoefs
    iadj <- solve(adj)
    v.var <- iadj %*% arvar %*% t(iadj)
    
    # Reverse the mapping for the variance estimate.
    v.var <- Q%*%v.var%*%t(Q)
    
    v[!novar,!novar] <- v.var
  }
  v
}

