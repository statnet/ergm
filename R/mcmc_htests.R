#  File R/approx.hotelling.diff.test.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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


#' Degrees of freedom for a Hotelling's \eqn{T^2} test
#'
#' One-sample, two-sample homoscedastic, and two-sample
#' heteroscedastic scenarios are handled, the latter using the
#' Krishnamoorthy and Yu (2004, eq. 7) degrees of freedom formula.
#'
#' In case some of the inputs or results are NaN, conservative
#' estimates are used.
#'
#' @param n (effective) sample size for the sample(s), scalar for
#'   one-sample, vector of 2 for two-sample.
#' @param v for unpooled, a list of 2 estimated variance-covariance
#'   matrices.
#'
#' @returns A scalar giving the degrees of freedom.
#' @noRd
hotelling_t2_df <- function(n, v) {
  NANVL <- function(z, ifNAN) ifelse(is.nan(z), ifNAN, z)

  if (length(n) == 1) { # one-sample
    NANVL(n[[1]], 1) - 1
  } else if (is.null(v)) { # 2-sample pooled
    NANVL(n[[1L]], 1) + NANVL(n[[2L]], 1) - 2
  } else { # 2-sample unpooled
    tr <- function(x) sum(diag(as.matrix(x)))

    d1 <- qrssolve(v[[1L]] + v[[2L]], v[[1]])
    d2 <- diag(1, nrow(d1)) - d1

    p <- attr(d1, "rank")

    (p + p^2) / (
      NANVL((tr(d1 %*% d1) + tr(d1)^2) / n[[1L]], 0) +
      NANVL((tr(d2 %*% d2) + tr(d2)^2) / n[[2L]], 0)
    )
  }
}


#' Approximate Hotelling T^2-Test for One or Two Population Means
#'
#' A multivariate hypothesis test for a single population mean or a
#' difference between them. This version attempts to adjust for
#' multivariate autocorrelation in the samples and handles data that
#' are not full rank. \insertNoCite{Ho47m}{ergm}
#' \insertNoCite{KrYu04m}{ergm} \insertNoCite{Lu05n}{ergm}
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
#' @seealso [t.test()], [spectrum0.mvar()]
#' @note For [`mcmc.list`] input, the variance for this test is
#'   estimated with unpooled means. This is not strictly correct.
#' @references \insertAllCited{}
#'
#' @export approx.hotelling.diff.test
approx.hotelling.diff.test<-function(x,y=NULL, mu0=0, assume.indep=FALSE, var.equal=FALSE, ...){
  if(!is.mcmc.list(x))
    x <- mcmc.list(mcmc(as.matrix(x)))
  if(!is.null(y) && !is.mcmc.list(y))
    y <- mcmc.list(mcmc(as.matrix(y)))

  if(is.null(mu0)) mu0 <- rep(0,nvar(x))
  mu0 <- rep(mu0, length.out = nvar(x))

  if(is.null(y)) var.equal <- FALSE

  data <- compact(list(x = x, y = y))

  m <- map(data, colMeans.mcmc.list)
  n <- map_int(data, function(dat) sum(map_int(dat, nrow)))

  NO_AR_ERR <- "Unable to compute autocorrelation-adjusted standard errors."

  v <- if (var.equal) {
         if (assume.indep) {
           map(data, var.mcmc.list) |>
             map2(n - 1L, `*`) |>
             reduce(`+`) |>
             (`/`)(sum(n) - 2) |>
             list() |>
             rep.int(2L) |>
             setNames(c("x", "y"))
         } else {
           # If we are pooling variances *and* estimating
           # autocorrelation.
           map2(data, m, sweep.mcmc.list) |> # Center each around mean.
             unlist(recursive = FALSE) |> # Combine.
             spectrum0.mvar(...) |> # Estimate.
             ERRVL2(stop(NO_AR_ERR)) |> # Handle error.
             list() |> # Enclose in list.
             rep.int(2L) |> # Make 2 copies.
             setNames(c("x", "y"))
         }
       } else {
         if (assume.indep) {
           map(data, var.mcmc.list)
         } else {
           map(data, spectrum0.mvar, ...) |> ERRVL2(stop(NO_AR_ERR))
         }
       }

  infl <- map(v, attr, "infl") |> map_dbl(`%||%`, 1)
  v <- map(v, as.matrix) # Ensure v is always a matrix.

  neff <- map2_dbl(n, infl, `/`)

  # Here, vcov already incorporates the inflation due to
  # autocorrelation, so n, not neff.
  vcov.m <- map2(v, n, `/`)

  d <- m$x - (m$y %||% 0)
  vcov.d <- vcov.m$x + (vcov.m$y %||% 0)

  names(mu0) <- varnames(x)
  nv <- diag(vcov.d) == 0
  
  method <- paste0("Hotelling's ", NVL2(y, "Two", "One"), "-Sample",
                   if (var.equal) " Pooled"," T^2-Test",
                   if (!assume.indep) " with correction for autocorrelation")

  T2 <- try(xTAx_seigen((d - mu0), vcov.d), silent = TRUE)

  # If a statistic doesn't vary and doesn't match, return a 0 p-value:
  if (is(T2, "try-error")) {
    if (grepl("x is not in the span of A", T2, fixed = TRUE)) {
      warning("Observed differences are not in the span of the variation of the variables; not possible under the null, so p-value = 0.")
      T2 <- structure(+Inf, rank = attr(xTAx_seigen((d - d), vcov.d), "rank"))
    } else stop(attr(T2, "condition"))
  }

  p <- attr(T2, "rank")
  if (p == 0) stop("data are essentially constant")
  else if (p < nvar(x)) message("Data are not full-rank (", nvar(x),
                                " variables total, ", p,
                                " linearly independent).")
 
  names(T2) <- "T^2"
  pars <- c(param = p, df = hotelling_t2_df(neff, if (!var.equal) vcov.m))

  if (pars[1] >= pars[2]) warning("Effective degrees of freedom (", pars[2],
                                  ") must exceed the number of varying parameters (",
                                  pars[1], "). P-value will not be computed.")
  structure(
    list(statistic = T2, parameter = pars, method = method, null.value = mu0,
         p.value = if (pars[1] < pars[2]) .ptsq(T2, pars[1], pars[2], lower.tail = FALSE) else NA,
         alternative = "two.sided", estimate = d, covariance = vcov.d,
         covariance.x = vcov.m$x, covariance.y = vcov.m$y),
    class = "htest")
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
    ERRVL2(approx.hotelling.diff.test(x1,x2,var.equal=TRUE,...),
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
#' Its return value, divided by the total sample size, is the
#' estimated variance-covariance matrix of the sampling distribution
#' of the mean of `x` if `x` is one or more multivatriate time series
#' with AR(\eqn{p}) structure (\insertCite{Lu05n;nobrackets}{ergm},
#' Prop. 3.3), with \eqn{p} determined by AIC.
#'
#' @param x a matrix with observations in rows and variables in
#'   columns or a list of such matrices, which are then treated as
#'   independent realizations of the same AR process; [`mcmc.list`]
#'   objects are accepted but treated no differently from ordinary
#'   lists.
#' @param order.max maximum (or fixed) order for the AR model.
#' @param aic use AIC to select the order (up to `order.max`).
#' @param tol tolerance used in detecting multicollinearity. See Note
#'   below.
#' @param ... additional arguments to [ar()].
#'
#' @return A square matrix with dimension equalling to the number of
#'   columns of `x`, with an additional attribute `"infl"` giving the
#'   factor by which the effective sample size is reduced due to
#'   autocorrelation, according to the \insertCite{VaFl15m;textual}{ergm}
#'   estimate for ESS.
#' 
#' @note AR fitting can fail under multicollinearity. This is remedied
#'   as follows:
#'
#' 1. Standardize the variables.
#' 2. Use eigendecomposition to extract principal components.
#' 3. Standardize the principal components.
#' 4. Drop those components whose standard deviation differs from 1 by more than `tol`. This should filter out redundant components or those too numerically unstable.
#' 5. Fit a multivariate AR model (via OLS), and calculate the variance.
#' 6. Reverse the mapping in steps 1-4 to obtain the variance of the original data.
#' @export spectrum0.mvar
spectrum0.mvar <- function(x, order.max=NULL, aic=is.null(order.max), tol=.Machine$double.eps^0.5, ...){
  if (is.list(x)) {
    lens <- map_int(x, nrow)
    x <- do.call(rbind, x)
  } else lens <- NULL

  p <- ncol(x)
  v <- matrix(0,p,p)

  # Save the scale of each variable, then drop nonvarying and standardise.
  xscl <- apply(x, 2L, stats::sd)
  novar <- xscl / max(xscl) < tol
  p.var <- sum(!novar)
  x <- x[, !novar, drop = FALSE]
  xscl <- xscl[!novar]
  if(ncol(x) == 0) stop("All variables are constant.")
  x <- sweep(x, 2L, xscl, `/`, check.margin = FALSE)
  
  # Map the variables onto their principal components, dropping
  # redundant (linearly-dependent) dimensions.
  e <- eigen(cov(x), symmetric=TRUE)
  x <- x %*% e$vectors %*% diag(1 / sqrt(pmax(e$values, 0)), p.var)
  x.sd <- apply(x, 2L, sd, na.rm = TRUE)
  ind <- !is.na(x.sd) & abs(x.sd - 1) < tol
  x <- x[, ind, drop = FALSE]
  ind.var <- var(x)

  # This matrix applies the reverse of the above transformation to the
  # resulting covariance matrix.
  unmapper <- diag(xscl, p.var) %*% e$vectors[, ind, drop=FALSE] %*% diag(sqrt(e$values[ind]), sum(ind))

  # Convert back into a list.
  x <- if (!is.null(lens)) split_len(x, lens, margin = 1) else list(x)
  
  # Calculate the time-series variance of the mean on the PC scale.
  arfit <- fit_var_ols_multi(x, lags = order.max, aic = aic)
  
  arvar <- arfit$var.pred
  arcoefs <- arfit$ar
  arcoefs <- NVL2(dim(arcoefs), apply(arcoefs, 1:2, base::sum), sum(arcoefs))
  
  adj <- diag(1, nrow = NROW(arvar)) - arcoefs
  v.var <- sandwich_qrssolve(adj, arvar)

  infl <- exp((determinant(v.var)$modulus-determinant(ind.var)$modulus)/ncol(ind.var))
  
  # Reverse the mapping for the variance estimate.
  v.var <- xAxT(unmapper, v.var)
  
  v[!novar,!novar] <- v.var

  structure(v, infl = as.vector(infl), rank = sum(ind))
}


#' Use OLS to fit an up to an AR(L) model
#'
#' This is a more focused version of ar.ols() that can also handle
#' multiple independent chains.
#'
#' @param Xl a list of matrices with common number of columns.
#' @param lags maximum lag to calculate (or a list thereof)
#' @param aic whether to return all lag values (`FALSE`) or only the
#'   one with the best aic (`TRUE`)
#' @param intercept whether to include an intercept
#'
#' @return For most uses, the return format is consistent with ar().
#'
#' @noRd
fit_var_ols_multi <- function(Xl, lags = NULL, aic = FALSE, intercept = TRUE) {
  p   <- ncol(Xl[[1]])
  T_c <- sapply(Xl, nrow)
  Tmax <- max(T_c)

  ## Default lag choice
  NVL(lags) <- ceiling(10 * log10(Tmax))

  ## Handle excessive lags
  if (aic) {
    lag_eff <- min(max(lags), Tmax - 1)
    lags_to_fit <- 0:lag_eff
  } else {
    lags_to_fit <- lags
    lag_eff <- max(lags[lags_to_fit <= Tmax - 1], 0, na.rm = TRUE)
  }

  ## Observation counts
  n_obs <- sapply(0:max(lag_eff, 0), function(l)
    sum(pmax(T_c - l, 0))
    )

  ## ------------------------------------------------------------
  ## 1. Accumulate sufficient statistics
  ## ------------------------------------------------------------

  XtX <- matrix(0, p * lag_eff, p * lag_eff)
  XtY <- matrix(0, p * lag_eff, p)
  YtY <- matrix(0, p, p)

  oneZ <- if (intercept) matrix(0, 1, p * lag_eff)
  oneY <- if (intercept) matrix(0, 1, p)
  one1 <- if (intercept) 0

  for (X in Xl) {
    Tc <- nrow(X)
    YtY <- YtY + crossprod(X)

    if (intercept) {
      oneY <- oneY + colSums(X)
      one1 <- one1 + Tc
    }

    for (i in seq_len(lag_eff)) {
      if (i >= Tc) next
      rows <- (i + 1):Tc

      Xi <- X[rows - i, , drop = FALSE]
      Yi <- X[rows, , drop = FALSE]

      XtY[((i - 1) * p + 1):(i * p), ] <-
        XtY[((i - 1) * p + 1):(i * p), ] + crossprod(Xi, Yi)

      if (intercept)
        oneZ[, ((i - 1) * p + 1):(i * p)] <-
          oneZ[, ((i - 1) * p + 1):(i * p)] + colSums(Xi)
    }

    for (i in seq_len(lag_eff)) {
      if (i >= Tc) next
      for (j in seq_len(lag_eff)) {
        if (j >= Tc) next
        k <- max(i, j)
        rows <- (k + 1):Tc
        XtX[((i - 1) * p + 1):(i * p),
        ((j - 1) * p + 1):(j * p)] <-
          XtX[((i - 1) * p + 1):(i * p),
          ((j - 1) * p + 1):(j * p)] +
          crossprod(X[rows - i, , drop = FALSE],
                    X[rows - j, , drop = FALSE])
      }
    }
  }

  ## ------------------------------------------------------------
  ## 2. Fit requested lags
  ## ------------------------------------------------------------

  fits <- vector("list", length(lags_to_fit))

  for (k in seq_along(lags_to_fit)) {
    l <- lags_to_fit[k]

    if (!aic && l > Tmax - 1) {
      fits[[k]] <- list(aic = NaN)
      next
    }

    n <- n_obs[l + 1]
    df <- n - (p * l + if (intercept) p else 0)
    if (df <= 0) {
      fits[[k]] <- list(aic = NaN)
      next
    }

    ## VAR(0)
    if (l == 0) {
      Sigma <- YtY / df
      fits[[k]] <- list(
        ar = array(0, c(p, p, 0)),
        x.intercept = if (intercept) colMeans(do.call(rbind, Xl)),
        var.pred = Sigma,
        aic = n * determinant(Sigma)$modulus
      )
      next
    }

    idx <- 1:(p * l)

    fit <- tryCatch({

      XtX_l <- XtX[idx, idx, drop = FALSE]
      XtY_l <- XtY[idx, , drop = FALSE]

      if (intercept) {
        XtX_l <- rbind(
          cbind(XtX_l, t(oneZ[, idx, drop = FALSE])),
          c(oneZ[, idx, drop = FALSE], one1)
        )
        XtY_l <- rbind(XtY_l, oneY)
      }

      B <- qr.solve(XtX_l, XtY_l)

      RSS <- YtY - t(B) %*% XtY_l
      Sigma <- RSS / df

      aic <- n * log(det(Sigma)) +
        2 * (p^2 * l + if (intercept) p else 0)

      A_arr <- array(0, c(p, p, l))
      for (i in seq_len(l))
        A_arr[, , i] <- t(B[((i - 1) * p + 1):(i * p), ])

        list(
          ar = A_arr,
          x.intercept = if (intercept) B[p * l + 1, ] else NULL,
          var.pred = Sigma,
          aic = aic
        )

    }, error = function(e) list(aic = NaN))

    fits[[k]] <- fit
  }

  ## ------------------------------------------------------------
  ## 3. Return logic
  ## ------------------------------------------------------------

  if (aic)
    return(fits[[which.min(sapply(fits, `[[`, "aic"))]])

  if (length(lags_to_fit) == 1)
    return(fits[[1]])

  fits
}
