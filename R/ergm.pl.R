#  File R/ergm.pl.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

#' @rdname ergm.mple
#' @description \code{ergm.pl} is an even more internal workhorse
#'   function that prepares many of the components needed by
#'   \code{ergm.mple} for the regression rountines that are used to
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
ergm.pl<-function(nw, fd, m, theta.offset=NULL,
                    control, ignore.offset=FALSE,
                    verbose=FALSE) {
  on.exit(ergm_Cstate_clear())
  on.exit(PL_workspace_clear(), add=TRUE)

  state <- ergm_state(nw, model=m)
  d <- sum(fd)
  el <- as.edgelist(state)
  elfd <- as.rlebdm(el) & fd
  e <- sum(elfd)

  maxDyads <- if(is.function(control$MPLE.samplesize)) control$MPLE.samplesize(d=d, e=e) else control$MPLE.samplesize

  # Convert the new-style MPLE results into old-style.
  old_PL <- function(z){
    y0 <- which(z$y[1,] != 0)
    y1 <- which(z$y[2,] != 0)
    x <- t(z$x[,c(y0, y1), drop=FALSE])
    list(y = rep(c(0L,1L), c(length(y0), length(y1))),
         x = x,
         weightsvector = c(z$y[1,y0],z$y[2,y1])
         )
  }

  z <- .Call("MPLE_wrapper",
             state,
             # MPLE settings
             as.double(to_ergm_Cdouble(fd)),
             as.integer(maxDyads),
             PACKAGE="ergm")
  z <- old_PL(z)

  if (verbose) {
    message(paste("MPLE covariate matrix has", length(z$y), "rows."))
  }
  zy <- z$y
  wend <- as.numeric(z$weightsvector)
  xmat <- z$x
  colnames(xmat) <- param_names(m,canonical=TRUE)
  rm(z)

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(wend)<d && mean(zy)<1/2){
    if(verbose) message("A sparse network with too many unique dyads encountered. Using case-control MPLE.")
    # Strip out the rows associated with ties.
    wend <- wend[zy==0]
    xmat <- xmat[zy==0,,drop=FALSE]
    zy <- zy[zy==0]

    z <- .Call("MPLE_wrapper",
               state,
               # MPLE settings
               as.double(to_ergm_Cdouble(elfd)),
               as.integer(.Machine$integer.max), # maxDyads
               PACKAGE="ergm")
    z <- old_PL(z)

    if (verbose) {
      message(paste("MPLE covariate matrix has", length(z$y), "rows."))
    }
    zy.e <- z$y
    wend.e <- as.numeric(z$weightsvector)
    xmat.e <- z$x
    colnames(xmat.e) <- param_names(m,canonical=TRUE)
    rm(z)

    # Divvy up the sampling weight of the ties:
    wend.e <- wend.e / sum(wend.e) * e

    # Divvy up the sampling weight of the nonties:
    wend <- wend / sum(wend) * (d-e)

    zy <- c(zy,zy.e)
    wend <- c(wend, wend.e)
    xmat <- rbind(xmat, xmat.e)

    rm(zy.e, wend.e, xmat.e)
  }

  #
  # Adjust for the offset
  #

  xmat.full <- xmat

  if(any(m$etamap$offsettheta) && !ignore.offset){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    # Compute the offset's effect.
    foffset <- .multiply.with.inf(xmat[,m$etamap$offsetmap,drop=FALSE], 
                                  cbind(ergm.eta(theta.offset,m$etamap)[m$etamap$offsetmap]))
    
    # Remove offset covariate columns.
    xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE] 
    colnames(xmat) <- param_names(m,canonical=TRUE)[!m$etamap$offsettheta]
    # Now, iff a row's offset effect is infinite, then it carries no
    # further information whatsoever, so it should be dropped.
    xmat <- xmat[is.finite(foffset),,drop=FALSE]
    zy <- zy[is.finite(foffset)]
    wend <- wend[is.finite(foffset)]
    foffset <- foffset[is.finite(foffset)]
  }else{
    foffset <- rep(0, length=length(zy))
  }
  
  list(xmat=xmat, zy=zy, foffset=foffset, wend=wend,
       xmat.full=xmat.full)
}

PL_workspace_clear <- function(){
  .Call("MPLE_workspace_free", PACKAGE="ergm")
}
