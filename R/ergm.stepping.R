#  File R/ergm.stepping.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

#' Transform points represented by rows of `m` such that their
#' centroid is shifted `1-gamma.tr` of the way toward `x` and their
#' spread around their centroid is scaled by a factor of
#' `gamma.scl`. Both `gamma.tr` and `gamma.scl` can be vectors equal
#' in length to the number of columns of `m`.
#' @noRd
.shift_scale_points <- function(m, x, gamma.tr, gamma.scl=gamma.tr){
  m.c <- colMeans(m)
  m <- sweep_cols.matrix(m, m.c, disable_checks=TRUE)
  m.c <- c(m.c*gamma.tr + x*(1-gamma.tr))
  t(t(m) * gamma.scl + m.c)
}

## Given two matrices x1 and x2 with d columns (and any positive
## numbers of rows), find the greatest gamma<=steplength.max s.t., the
## points of x2 shrunk towards the centroid of x1 a factor of gamma,
## are all in the convex hull of x1, as is the centroid of x2 shrunk
## by margin*gamma.
##
## If -1 <= margin < 0, the algorithm "forgives" some amount of either
## centroid or x2 being outside of the convex hull of x1.

## This is a variant of Hummel et al. (2010)'s steplength algorithm
## also usable for missing data MLE.
.Hummel.steplength <- function(x1, x2=NULL, margin=0.05, steplength.max=1, x2.num.max=c(max(ncol(rbind(x1))*2, 30), 10), parallel=c("observational","never"), control=NULL, verbose=FALSE){
  parallel <- match.arg(parallel)
  margin <- 1 + margin
  x1 <- rbind(x1); m1 <- rbind(colMeans(x1)); ; n1 <- nrow(x1)
  if(is.function(x2.num.max)) x2.num.max <- x2.num.max(x1)
  if(is.null(x2)){
    m2 <- rbind(rep(0,ncol(x1)))
    parallel <- FALSE
  }else{                                      
    x2 <- rbind(x2)
    m2 <- rbind(colMeans(x2))
    parallel <- parallel == "observational"
  }
  n2 <- nrow(x2)

  # Drop duplicated elements in x1, *and* those elements in x2 that
  # duplicate those in x1.
  d12 <- duplicated(rbind(x1,x2))
  d1 <- d12[seq_len(n1)]
  d2 <- d12[-seq_len(n1)]
  x1 <- x1[!d1,,drop=FALSE]
  x2 <- x2[!d2,,drop=FALSE]
  if(length(x2)==0) x2 <- NULL

  if(verbose>1) message("Eliminating repeated points: ", sum(d1),"/",length(d1), " from target set, ", sum(d2),"/",length(d2)," from test set.")

  cl <- if(parallel && !is.null(control)) ergm.getCluster(control, verbose)

  ## Use PCA to rotate x1 into something numerically stable and drop
  ## unused dimensions, then apply the same affine transformation to
  ## m1 and x2:
  if(nrow(x1)>1){
    ## Center:
    x1m <- colMeans(x1) # note that colMeans(x)!=m1
    x1c <- sweep(x1, 2, x1m, "-")
    ## Rotate x1 onto its principal components, dropping linearly dependent dimensions:
    e <- eigen(crossprod(x1c), symmetric=TRUE)
    Q <- e$vec[,sqrt(pmax(e$val,0)/max(e$val))>sqrt(.Machine$double.eps)*2,drop=FALSE]
    x1cr <- x1c%*%Q # Columns of x1cr are guaranteed to be linearly independent.

    ## Scale x1:
    x1crsd <- pmax(apply(x1cr, 2, sd), sqrt(.Machine$double.eps))
    x1crs <- sweep(x1cr, 2, x1crsd, "/")

    ## Now, apply these operations to m1 and x2:
    m1crs <- sweep(sweep(m1, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
    if(!is.null(x2)) x2crs <- sweep(sweep(x2, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
    m2crs <- sweep(sweep(m2, 2, x1m, "-")%*%Q, 2, x1crsd, "/")
  }else{
    if(is.null(x2)){
      if(isTRUE(all.equal(m1,m2,check.attributes=FALSE))) return(1) else return(0)
    }else{
      if(apply(x2, 1, all.equal, m1, check.attributes=FALSE) %>% map_lgl(isTRUE) %>% all) return(1) else return(0)
    }
  }

  if(!is.null(x2) && nrow(x2crs) > x2.num.max[1]){
    ## If constrained sample size > x2.num.max
    if(verbose>1){message("Using fast and approximate Hummel et al search.")}
    d <- rowSums(sweep(x2crs, 2, m1crs)^2)
    x2crs <- x2crs[order(-d)[1:x2.num.max[1]],,drop=FALSE]
  }

  max_run <- if(length(x2.num.max) > 1) x2.num.max[2] else Inf

  if(!is.null(cl)){
    # NBs: parRapply() would call shrink_into_CH() for every
    # row. Direct reference to split.data.frame() is necessary here
    # since no matrix method.
    x2crs <- split.data.frame(x2crs, rep_len(seq_along(cl), nrow(x2crs)))
    min(steplength.max, unlist(parallel::parLapply(cl=cl, x2crs, shrink_into_CH, x1crs, verbose=verbose, max_run=max_run))/margin)
  }else{
    min(steplength.max, shrink_into_CH(if(!is.null(x2)) x2crs else m2crs, x1crs, verbose=verbose, max_run=max_run)/margin)
  }
}
