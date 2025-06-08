#  File R/ergm.pl.R in package ergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' @rdname ergm.mple
#' @description \code{ergm.pl} is an even more internal workhorse
#'   function that prepares many of the components needed by
#'   \code{ergm.mple} for the regression routines that are used to
#'   find the MPLE estimated ergm. It should not be called directly by
#'   the user.
#'
#' @param theta.offset a numeric vector of length equal to the number
#'   of statistics of the model, specifying (positionally) the
#'   coefficients of the offset statistics; elements corresponding to
#'   free parameters are ignored.
#' @param ignore.offset If \code{FALSE} (the default), columns
#'   corresponding to terms enclosed in \code{offset()} are not
#'   returned with others but are instead processed by multiplying
#'   them by their corresponding coefficients (which are fixed, by
#'   virtue of being offsets) and the results stored in a separate
#'   column.
#'
#' @return \code{ergm.pl} returns a list containing:
#'
#' \item{xmat}{the compressed and possibly sampled matrix of change
#'   statistics}
#'
#' \item{xmat.full}{as `xmat` but with offset terms}
#'
#' \item{zy}{the corresponding vector of responses, i.e. tie values}
#'
#' \item{foffset}{if `ignore.offset==FALSE`, the combined offset statistics multiplied by their parameter values}
#'
#' \item{wend}{the vector of weights for `xmat`
#'   and `zy`}
#'
#' @keywords internal
#' @export
ergm.pl<-function(state, state.obs, theta.offset=NULL,
                    control, ignore.offset=FALSE,
                    verbose=FALSE) {
  on.exit(ergm_Cstate_clear())
  on.exit(PL_workspace_clear(), add=TRUE)

  m <- state$model
  fd <- if(is(state.obs, "rlebdm")) state.obs
        else as.rlebdm(state$proposal$arguments$constraints, state.obs$proposal$arguments$constraints, which="informative")
  d <- sum(fd)
  el <- as.edgelist(state)
  elfd <- as.rlebdm(el) & fd
  e <- sum(elfd)

  maxDyads <- if(is.function(control$MPLE.samplesize)) control$MPLE.samplesize(d=d, e=e) else control$MPLE.samplesize

  z <- .Call("MPLE_wrapper",
             state,
             # MPLE settings
             as.double(to_ergm_Cdouble(fd)),
             as.integer(maxDyads),
             PACKAGE="ergm")
  y <- z$y
  x <- z$x
  rm(z)
  rownames(x) <- param_names(m,canonical=TRUE)

  if(verbose) message(paste("MPLE covariate matrix has",ncol(y), "rows."))

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(y)<d && sum(y[1,])/sum(y)<1/2){
    if(verbose) message("A sparse network with too many unique dyads encountered. Using case-control MPLE.")
    # Strip out the rows associated with ties.
    tokeep <- y[2,] != 0
    x <- x[, tokeep, drop=FALSE]
    y <- y[, tokeep, drop=FALSE]
    y[1,] <- 0L

    z <- .Call("MPLE_wrapper",
               state,
               # MPLE settings
               as.double(to_ergm_Cdouble(elfd)),
               as.integer(.Machine$integer.max), # maxDyads
               PACKAGE="ergm")

    y.e <- z$y
    x.e <- z$x
    rm(z)
    rownames(x.e) <- param_names(m,canonical=TRUE)

    if(verbose) message(paste("MPLE covariate matrix has", ncol(y), "rows."))

    # Divvy up the sampling weight of the ties:
    y.e <- y.e / sum(y.e) * e

    # Divvy up the sampling weight of the nonties:
    y <- y / sum(y) * (d-e)

    y <- cbind(y, y.e)
    x <- cbind(x, x.e)
  }

  #
  # Adjust for the offset
  #

  x.full <- x

  if(any(m$etamap$offsettheta) && !ignore.offset){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    # Compute the offset's effect.
    foffset <- .multiply.with.inf(t(x[m$etamap$offsetmap,,drop=FALSE]),
                                  cbind(ergm.eta(theta.offset,m$etamap)[m$etamap$offsetmap]))
    
    # Remove offset covariate columns.
    x <- x[!m$etamap$offsettheta,,drop=FALSE]
    # Now, iff a row's offset effect is infinite, then it carries no
    # further information whatsoever, so it should be dropped.
    x <- x[,is.finite(foffset),drop=FALSE]
    y <- y[,is.finite(foffset),drop=FALSE]
    foffset <- foffset[is.finite(foffset)]
  }else{
    foffset <- rep(0, length=ncol(y))
  }

  # Convert the new-style MPLE results into old-style.
  y0 <- which(y[2,] != 0)
  y1 <- which(y[1,] != 0)
  x <- t(x[,c(y0, y1), drop=FALSE])
  x.full <- t(x.full[,c(y0, y1), drop=FALSE])
  foffset <- foffset[c(y0,y1)]
  list(zy = rep(c(0L,1L), c(length(y0), length(y1))),
       xmat = x,
       foffset=foffset,
       xmat.full=x.full,
       wend = c(y[2,y0],y[1,y1])
       )
}

PL_workspace_clear <- function(){
  .Call("MPLE_workspace_free", PACKAGE="ergm")
}
