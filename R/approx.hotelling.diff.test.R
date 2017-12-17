#  File R/approx.hotelling.diff.test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
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

  mu0 <- rep(mu0, length.out = nvar(x))

  tr <- function(x) sum(diag(as.matrix(x)))

  vars <- list(x=list(v=x))
  if(!is.null(y)) vars$y <- list(v=y)
  
  mywithin <- function(data, ...) within(data, ...) # This is a workaround suggsted by Duncan Murdoch: calling lapply(X, within, {CODE}) would leave CODE unable to see any objects in f.
  vars <- lapply(vars, mywithin, {
    vm <- as.matrix(v)
    vcov.indep <- cov(vm)
    if(assume.indep){
      vcov <- vcov.indep
    }else{
      vcov <- spectrum0.mvar(v)
      vcov[is.na(vcov)] <- 0
    }
    m <- colMeans(vm)
    n <- nrow(vm)
    
    infl <- tr(vcov) / tr(vcov.indep) # I.e., how much bigger is the trace of the variance-covariance after taking autocorrelation into account than before.
    neff <- n / infl
    
    vcov.m <- vcov/n # Here, vcov already incorporates the inflation due to autocorrelation.
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


  p <- nvar(x$v)
  names(mu0)<-varnames(x$v)

  novar <- diag(vcov.d)==0
  p <- p-sum(novar)

  if(p==0) stop("data are essentially contstant")
  
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
              covariance = vcov.d,
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

# Compute the sample estimating equations of an ERGM. If the model is
# linear, all non-offset statistics are passed. If the model is
# curved, the score estimating equations (3.1) by Hunter and
# Handcock (2006) are given instead.
.ergm.esteq <- function(theta, model, statsmatrix){
  etamap <- NVL(model$etamap, model)
  esteq <- t(ergm.etagradmult(theta,t(as.matrix(statsmatrix)),etamap))[,!etamap$offsettheta,drop=FALSE]
  if(is.mcmc(statsmatrix)){
    esteq <- mcmc(esteq, start=start(statsmatrix), end=end(statsmatrix), thin=thin(statsmatrix))
    varnames(esteq) <- NVL(names(theta), coef.names.model(model, FALSE))[!etamap$offsettheta]
  }else  colnames(esteq) <- NVL(names(theta), coef.names.model(model, FALSE))[!etamap$offsettheta]

  esteq
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
  breaks <- if(is.mcmc.list(x)) c(0,cumsum(sapply(x, niter))) else 0
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  v <- matrix(NA,p,p)
  novar <- abs(apply(x,2,stats::sd))<tol
  x <- x[,!novar,drop=FALSE]

  # Index of the first local minimum in a sequence.
  first_local_min <- function(x){
    d <- diff(c(Inf,x,Inf))
    min(which(d>=0))-1
  }
  
  # Map the variables onto their principal components, dropping
  # redundant (linearly-dependent) dimensions. Here, we keep the
  # eigenvectors such that the reciprocal condition number defined
  # as s.min/s.max, where s.min and s.max are the smallest and the
  # biggest singular values, respectively, is greater than the
  # tolerance.
  e <- eigen(cov(x), symmetric=TRUE)
  Q <- e$vec[,sqrt(pmax(e$val,0)/max(e$val))>tol*2,drop=FALSE]
  xr <- x%*%Q # Columns of xr are guaranteed to be linearly independent.

  # Convert back into an mcmc.list object.
  xr <- do.call(mcmc.list,lapply(lapply(seq_along(breaks[-1]), function(i) xr[(breaks[i]+1):(breaks[i+1]),,drop=FALSE]), mcmc))
  
  # Calculate the time-series variance of the mean on the PC scale.
  ord <- NVL(order.max, ceiling(10*log10(niter(xr))))
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
  
  adj <- diag(1,nrow=nvar(xr)) - arcoefs
  iadj <- solve(adj)
  v.var <- iadj %*% arvar %*% t(iadj)
  
  # Reverse the mapping for the variance estimate.
  v.var <- Q%*%v.var%*%t(Q)
  
  v[!novar,!novar] <- v.var
  
  v
}

## Note: This code is derived from the stats package. It incorporates
## patches that have been submitted to the R bugzilla. It uses the NA
## padding trick to get the ar estimate for the combined MCMC samples.

#' An ar.yw method for [`mcmc.list`] objects.
#'
#' @param x an [`mcmc.list`] object.
#' @param aic,order.max,na.action,demean,series,... See [ar()]. 
#'
#' @method ar.yw mcmc.list
#' @export
ar.yw.mcmc.list <-
    function (x, aic = TRUE, order.max = NULL, na.action = na.fail,
              demean = TRUE, series = NULL, ...)
{

  x <- lapply(x, na.action)
  na.action <- na.pass
  x <- do.call(rbind, lapply(x, function(z) rbind(cbind(z), matrix(NA, niter(z), nvar(z)))))
  ## Note: Assuming that the patch for stats is accepted, the
  ## following code can be replaced with a call to ar.yw().
    if(is.null(series)) series <- deparse(substitute(x))
    xfreq <- frequency(x)
    x <- as.matrix(x)
    if(!is.numeric(x))
        stop("'x' must be numeric")
    if(any(is.na(x)!=is.na(x[,1]))) stop("NAs 'x' must be the the same row-wise")
    nser <- ncol(x)
    if (demean) {
        xm <- colMeans(x, na.rm=TRUE)
        x <- sweep(x, 2L, xm, check.margin=FALSE)
    } else xm <- rep.int(0, nser)
    n.used <- nrow(x)
    n.obs <- sum(!is.na(x[,1])) # number of non-missing rows
    order.max <- if (is.null(order.max))
	min(n.obs - 1L, floor(10 * log10(n.obs))) else round(order.max)
    if (order.max < 1L) stop("'order.max' must be >= 1")
    else if (order.max >= n.obs) stop("'order.max' must be < 'n.obs'")
    xacf <- acf(x, type = "covariance", lag.max = order.max, plot = FALSE,
                demean = demean, na.action=na.pass)$acf
    if(nser > 1L) {
        ## multivariate case
        snames <- colnames(x)
        A <- B <- array(0, dim = c(order.max + 1L, nser, nser))
        A[1L, , ] <- B[1L, , ] <- diag(nser)
        EA <- EB <- xacf[1L, , , drop = TRUE]
        partialacf <- array(dim = c(order.max, nser, nser))
        xaic <- numeric(order.max + 1L)
        solve.yw <- function(m) {
            # Solve Yule-Walker equations with Whittle's
            # generalization of the Levinson(-Durbin) algorithm
            betaA <- betaB <- 0
            for (i in 0L:m) {
                betaA <- betaA + A[i + 1L, , ] %*% xacf[m + 2L - i, , ]
                betaB <- betaB + B[i + 1L, , ] %*% t(xacf[m + 2L - i, , ])
            }
            KA <- -t(qr.solve(t(EB), t(betaA)))
            KB <- -t(qr.solve(t(EA), t(betaB)))
            EB <<- (diag(nser) - KB %*% KA) %*% EB
            EA <<- (diag(nser) - KA %*% KB) %*% EA
            Aold <- A
            Bold <- B
            for (i in seq_len(m + 1L)) {
                A[i + 1L, , ] <<- Aold[i + 1L, , ] + KA %*% Bold[m + 2L - i, , ]
                B[i + 1L, , ] <<- Bold[i + 1L, , ] + KB %*% Aold[m + 2L - i, , ]
            }
        }
        cal.aic <- function() { # omits mean params, that is constant adj
            det <- abs(prod(diag(qr(EA)$qr)))
            return(n.obs * log(det) + 2 * m * nser * nser)
        }
        cal.resid <- function() {
            resid <- array(0, dim = c(n.used - order, nser))
            for (i in 0L:order)
                resid <- resid + x[(order - i + 1L):(n.used - i), , drop = FALSE] %*% t(ar[i + 1L, , ])
            return(rbind(matrix(NA, order, nser), resid))
        }
        order <- 0L
        for (m in 0L:order.max) {
            xaic[m + 1L] <- cal.aic()
            if (!aic || xaic[m + 1L] == min(xaic[seq_len(m + 1L)])) {
                ar <- A
                order <- m
                var.pred <- EA * n.obs/(n.obs - nser * (m + 1L))
            }
            if (m < order.max) {
                solve.yw(m)
                partialacf[m + 1L, , ] <- -A[m + 2L, , ]
            }
        }
        xaic <- setNames(xaic - min(xaic), 0L:order.max)
        resid <- cal.resid()
        if(order) {
            ar <- -ar[2L:(order + 1L), , , drop = FALSE]
            dimnames(ar) <- list(seq_len(order), snames, snames)
        } else ar <- array(0, dim = c(0L, nser, nser),
                           dimnames = list(NULL, snames, snames))
        dimnames(var.pred) <- list(snames, snames)
        dimnames(partialacf) <- list(seq_len(order.max), snames, snames)
        colnames(resid) <- colnames(x)
    } else {
        if (xacf[1L] == 0) stop("zero-variance series")
        ## univariate case
        r <- as.double(drop(xacf))
        z <- .Fortran(C_eureka, as.integer(order.max), r, r,
                      coefs = double(order.max^2),
                      vars = double(order.max),
                      double(order.max))
        coefs <- matrix(z$coefs, order.max, order.max)
        partialacf <- array(diag(coefs), dim = c(order.max, 1L, 1L))
        var.pred <- c(r[1L], z$vars)
        xaic <- n.obs * log(var.pred) + 2 * (0L:order.max) + 2 * demean
        maic <- min(aic)
	xaic <- setNames(if(is.finite(maic)) xaic - min(xaic) else
			 ifelse(xaic == maic, 0, Inf),
			 0L:order.max)
        order <- if (aic) (0L:order.max)[xaic == 0L] else order.max
        ar <- if (order) coefs[order, seq_len(order)] else numeric()
        var.pred <- var.pred[order + 1L]
        ## Splus compatibility fix
        var.pred <- var.pred * n.obs/(n.obs - (order + 1L))
        resid <- if(order) c(rep.int(NA, order), embed(x, order + 1L) %*% c(1, -ar))
        else as.vector(x) # we had as.matrix() above
    }
    res <- list(order = order, ar = ar, var.pred = var.pred, x.mean  =  drop(xm),
                aic  =  xaic, n.used = n.used, n.obs = n.obs, order.max = order.max,
                partialacf = partialacf, resid = resid, method = "Yule-Walker",
                series = series, frequency = xfreq, call = match.call())
    if(nser == 1L && order)
        res$asy.var.coef <-
            solve(toeplitz(drop(xacf)[seq_len(order)]))*var.pred/n.obs
    class(res) <- "ar"
    res
}
