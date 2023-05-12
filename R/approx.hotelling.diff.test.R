#  File R/approx.hotelling.diff.test.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
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
  fq <- qf(p, param, df - param + 1, lower.tail=lower.tail, log.p=log.p)
  fq / ((df - param + 1)/(param*df))
}



#' Approximate Hotelling T^2-Test for One or Two Population Means
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
#' @param ... additional arguments, passed on to [spectrum0.mvar()],
#'   etc.; in particular, `order.max=` can be used to limit the order
#'   of the AR model used to estimate the effective sample size.
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
#' \item{covariance.x}{Estimated variance-covariance matrix of the estimate of the mean of `x`.}
#' \item{covariance.y}{Estimated variance-covariance matrix of the estimate of the mean of `y`.}
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
approx.hotelling.diff.test<-function(x,y=NULL, mu0=0, assume.indep=FALSE, var.equal=FALSE, ...){
  if(!is.mcmc.list(x))
    x <- mcmc.list(mcmc(as.matrix(x)))
  if(!is.null(y) && !is.mcmc.list(y))
    y <- mcmc.list(mcmc(as.matrix(y)))

  if(is.null(mu0)) mu0 <- rep(0,nvar(x))
  mu0 <- rep(mu0, length.out = nvar(x))

  tr <- function(x) sum(diag(as.matrix(x)))

  vars <- list(x=list(v=x))
  if(!is.null(y)) vars$y <- list(v=y)
  else var.equal <- FALSE

  vars <- lapply(if(is.null(y)) list(x=x) else list(x=x,y=y), function(v, ...){
    vm <- as.matrix(v)
    vcov.indep <- cov(vm)
    if(assume.indep){
      vcov <- vcov.indep
      infl <- 1
    }else if(!var.equal){
      vcov <- ERRVL(try(spectrum0.mvar(v, ...), silent=TRUE),
                    stop("Unable to compute autocorrelation-adjusted standard errors."))
      infl <- attr(vcov, "infl")
    }else{
      infl <- vcov <- NULL
    }
    m <- colMeans(vm)
    n <- nrow(vm)
    
    neff <- n / infl
    vcov.m <- vcov/n # Here, vcov already incorporates the inflation due to autocorrelation.

    list(v=v, vm=vm, m=m, n=n, vcov.indep=vcov.indep, vcov=vcov, infl=infl, neff=neff, vcov.m=vcov.m)
  }, ...)
  
  x <- vars$x
  y <- vars$y
  
  d <- x$m
  vcov.d <- x$vcov.m
  
  if(!is.null(y)){
    d <- d - y$m
    if(var.equal){
      # If we are pooling variances *and* estimating autocorrelation, then pool the two variables before calling spectrum0.mvar().
      if(!assume.indep){
        # Center both variables about the same mean.
        xp <- sweep.mcmc.list(x$v, x$m)
        yp <- sweep.mcmc.list(y$v, y$m)

        if(nchain(xp) == nchain(yp)){
          # Equal numbers of chains -> concatenate.
          xyp <- as.mcmc.list(Map(function(x,y) as.mcmc(rbind(x, matrix(NA, ceiling(10*log10(niter(xp) + niter(yp))), nvar(xp)), y)), xp, yp))
        }else{
          # Differing numbers of chains -> pad the shorter of the variables (if any) with NAs.
          padto <- max(niter(x$v), niter(y$v))
          xp <-
            if((xpad <- padto - niter(x$v)) > 0) lapply.mcmc.list(xp, function(z) rbind(z, matrix(NA, xpad, nvar(x$v))))
            else lapply.mcmc.list(xp, function(z) as.matrix(xp)) # To avoid MCMC burnin/interval inconsistency.
          yp <-
            if((ypad <- padto - niter(y$v)) > 0) lapply.mcmc.list(yp, function(z) rbind(z, matrix(NA, xpad, nvar(y$v))))
            else lapply.mcmc.list(yp, function(z) as.matrix(xp))  # To avoid MCMC burnin/interval inconsistency.

          # Now simply treat them as different chains of an MCMC sample.
          xyp <- as.mcmc.list(c(unclass(xp), unclass(yp)))
        }

        vcov <- x$vcov <- y$vcov <- ERRVL(try(spectrum0.mvar(xyp, ...), silent=TRUE),
                                  stop("Unable to compute autocorrelation-adjusted standard errors."))
        infl <- x$infl <- y$infl <- attr(vcov, "infl")
        x$neff <- x$n / infl
        y$neff <- y$n / infl
        x$vcov.m <- vcov / x$n
        y$vcov.m <- vcov / y$n
      }

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

  if(p==0) stop("data are essentially constant")
  
  ivcov.d <-sginv(vcov.d[!novar,!novar,drop=FALSE])
  
  method <- paste0("Hotelling's ",
                  NVL2(y, "Two", "One"),
                  "-Sample",if(var.equal) " Pooled"," T^2-Test", if(!assume.indep) " with correction for autocorrelation")
  
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
    T2 <- xTAx((d-mu0)[!novar], ivcov.d)
  }

  NANVL <- function(z, ifNAN) ifelse(is.nan(z),ifNAN,z)
  
  names(T2) <- "T^2"
  pars <- c(param = p, df = if(is.null(y)){
    NANVL(x$neff,1)-1
  }else if(var.equal){
    NANVL(x$neff,1)+NANVL(y$neff,1)-2
  }else{
    # This is the Krishnamoorthy and Yu (2004) degrees of freedom formula, courtesy of Wikipedia.
    (p+p^2)/(
      NANVL((tr(x$vcov.m[!novar,!novar] %*% ivcov.d %*% x$vcov.m[!novar,!novar] %*% ivcov.d) +
             tr(x$vcov.m[!novar,!novar] %*% ivcov.d)^2)/x$neff,0) +
      NANVL((tr(y$vcov.m[!novar,!novar] %*% ivcov.d %*% y$vcov.m[!novar,!novar] %*% ivcov.d) +
             tr(y$vcov.m[!novar,!novar] %*% ivcov.d)^2)/y$neff,0))
  })

  if(pars[1]>=pars[2]) warning("Effective degrees of freedom (",pars[2],") must exceed the number of varying parameters (",pars[1],"). P-value will not be computed.")
  out <- list(statistic=T2, parameter=pars, p.value=if(pars[1]<pars[2]) .ptsq(T2,pars[1],pars[2],lower.tail=FALSE) else NA,
              method = method,
              null.value=mu0,
              alternative="two.sided",
              estimate = d,
              covariance = vcov.d,
              covariance.x = x$vcov.m,
              covariance.y = y$vcov.m,
              novar = novar)
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
#' @param ... additional arguments, passed on to
#'   [approx.hotelling.diff.test()], which passes them to
#'   [spectrum0.mvar()], etc.; in particular, `order.max=` can be used
#'   to limit the order of the AR model used to estimate the effective
#'   sample size.
#' @note If [approx.hotelling.diff.test()] returns an error, then
#'   assume that burn-in is insufficient.
#' @return An object of class `htest`, inheriting from that returned
#'   by [approx.hotelling.diff.test()], but with p-value considered to
#'   be 0 on insufficient sample size.
#'
#' @seealso [coda::geweke.diag()], [approx.hotelling.diff.test()]
#' @export
geweke.diag.mv <- function(x, frac1 = 0.1, frac2 = 0.5, split.mcmc.list = FALSE, ...){
  # The following function's bookkeeping parts (e.g., handling of
  # mcmc.list and calculation of windows starts and ends) are loosely
  # based on parts of geweke.diag() from the coda R package.
  
  if(is.mcmc.list(x)){
    if(split.mcmc.list){
      return(lapply(x, geweke.diag.mv, frac1, frac2, ...))
    }
  }else{
    x <- as.mcmc(x)
  }
  
  x.len <- end(x) - start(x)
  x1 <- window(x, start=start(x), end=start(x) + frac1*x.len)
  x2 <- window(x, start=end(x) - frac2*x.len, end=end(x))

  test <-
    ERRVL(try(approx.hotelling.diff.test(x1,x2,var.equal=TRUE,...), silent=TRUE),
    {
      warning("Multivariate Geweke diagnostic failed, probably due to insufficient sample size.", call.=FALSE, immediate.=TRUE)
      test <- structure(list(p.value=NA), class="htest")
    })
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
#' @return A square matrix with dimension equalling to the number of
#'   columns of `x`, with an additional attribute `"infl"` giving the
#'   factor by which the effective sample size is reduced due to
#'   autocorrelation, according to the Vats, Flegal, and Jones (2015)
#'   estimate for ESS.
#' 
#' @note [ar()] fails if `crossprod(x)` is singular,
#' which is remedied by mapping the variables onto the principal
#' components of `x`, dropping redundant dimentions.
#' @export spectrum0.mvar
spectrum0.mvar <- function(x, order.max=NULL, aic=is.null(order.max), tol=.Machine$double.eps^0.5, ...){
  breaks <- if(is.mcmc.list(x)) c(0,cumsum(sapply(x, niter))) else NULL
  x.full <- as.matrix(x)
  x <- na.omit(x.full)
  p <- ncol(x)
  
  v <- matrix(0,p,p)

  # Save the scale of each variable, then drop nonvarying and standardise.
  # TODO: What if the variable actually has a tiny magnitude?
  xscl <- apply(x, 2L, stats::sd)
  novar <- xscl < tol
  x <- x[,!novar,drop=FALSE]
  x.full <- x.full[,!novar,drop=FALSE]
  xscl <- xscl[!novar]
  if(ncol(x) == 0) stop("All variables are constant.")
  x <- sweep(x, 2L, xscl, "/", check.margin = FALSE)
  x.full <- sweep(x.full, 2L, xscl, "/", check.margin = FALSE)

  # Index of the first local minimum in a sequence.
  first_local_min <- function(x){
    d <- diff(c(Inf,x,Inf))
    min(which(d>=0))-1
  }
  
  # Map the variables onto their principal components, dropping
  # redundant (linearly-dependent) dimensions. Here, we keep the
  # eigenvectors such that the reciprocal condition number defined as
  # s.min/s.max, where s.min and s.max are the smallest and the
  # biggest eigenvalues, respectively, is greater than the tolerance.
  e <- eigen(cov(x), symmetric=TRUE)
  Q <- e$vectors[,sqrt(pmax(e$values,0)/max(e$values))>tol*2,drop=FALSE]
  xr <- x.full%*%Q # Columns of xr are guaranteed to be linearly independent.

  ind.var <- cov(xr, use="complete.obs") # Get the sample variance of the transformed columns.

  # Convert back into an mcmc.list object.
  xr <-
    if(!is.null(breaks)) do.call(mcmc.list,lapply(lapply(seq_along(breaks[-1]), function(i) xr[(breaks[i]+1):(breaks[i+1]),,drop=FALSE]), mcmc))
    else as.mcmc.list(mcmc(xr))
  
  ord <- NVL(order.max, ceiling(10*log10(niter(xr))))
  xr <- do.call(rbind, c(lapply(unclass(xr)[-nchain(xr)], function(z) rbind(cbind(z), matrix(NA, ord, nvar(z)))), unclass(xr)[nchain(xr)]))

  safe_ar <- possibly(ar, NULL)

  # Calculate the time-series variance of the mean on the PC scale.
  # If ar() failed or produced a variance matrix estimate that's
  # not positive semidefinite, try with a lower order.
  arfit <- NULL; ord <- ord + 1
  while((is.null(arfit) || ERRVL(try(any(eigen(arfit$var.pred, only.values=TRUE)$values<0), silent=TRUE), TRUE)) && ord > 0){
    ord <- ord - 1
    if(ord<=0) stop("Unable to fit ar() even with order 1; this is likely to be due to insufficient sample size or a trend in the data.")
    arfit <- safe_ar(xr,aic=is.null(order.max), order.max=ord, na.action=na.pass, ...)
  }
  
  if(aic && arfit$order>(ord <- first_local_min(arfit$aic)-1)){
    arfit <- ar(xr, aic=ord==0, order.max=max(ord,1), na.action=na.pass) # Workaround since ar() won't take order.max=0.
  }
  
  arvar <- arfit$var.pred
  arcoefs <- arfit$ar
  arcoefs <- NVL2(dim(arcoefs), apply(arcoefs,2:3,base::sum), sum(arcoefs))
  
  adj <- diag(1,nrow=ncol(xr)) - arcoefs
  v.var <- sandwich_ssolve(adj, arvar)

  infl <- exp((determinant(v.var)$modulus-determinant(ind.var)$modulus)/ncol(ind.var))
  
  # Reverse the mapping for the variance estimate.
  v.var <- xAxT(xscl*Q, v.var)
  
  v[!novar,!novar] <- v.var
  
  attr(v, "infl") <- infl
  v
}
