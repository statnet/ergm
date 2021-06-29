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

  state <- ergm_state(nw, model=m)
  d <- sum(fd)
  el <- as.edgelist(state)
  elfd <- as.rlebdm(el) & fd
  e <- sum(elfd)

  maxNumDyadTypes <- as.integer(min(if(is.function(control$MPLE.max.dyad.types)) control$MPLE.max.dyad.types(d=d, e=e) else control$MPLE.max.dyad.types,
                         1.05*d)) # a little larger than d so the hash table doesn't bog down
  maxDyads <- if(is.function(control$MPLE.samplesize)) control$MPLE.samplesize(d=d, e=e) else control$MPLE.samplesize

  if(as.double(maxNumDyadTypes)*nparam(m, canonical=TRUE) > .Machine$integer.max) {
    stop("The maximum number of unique dyad types times the number of statistics exceeds 32 bit limits, so the MPLE cannot proceed; try reducing either MPLE.max.dyad.types or the number of terms in the model.")
  }

  z <- .Call("MPLE_wrapper",
             state,
             # MPLE settings
             as.double(to_ergm_Cdouble(fd)),
             as.integer(maxDyads),
             as.integer(maxNumDyadTypes),
             PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  if (verbose) {
    message(paste("MPLE covariate matrix has", sum(uvals), "rows."))
  }
  zy <- z$y[uvals]
  wend <- as.numeric(z$weightsvector[uvals])
  xmat <- matrix(z$x, ncol=nparam(m,canonical=TRUE), byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- param_names(m,canonical=TRUE)
  rm(z,uvals)

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(wend)<d && mean(zy)<1/2){
    if(verbose) message("A sparse network with too many unique dyads encountered. Using case-control MPLE.")
    # Strip out the rows associated with ties.
    wend <- wend[zy==0]
    xmat <- xmat[zy==0,,drop=FALSE]
    zy <- zy[zy==0]

    ## Run a whitelist PL over all of the toggleable edges in the network.
    maxNumDyadTypes <- min(maxNumDyadTypes, e)

    z <- .Call("MPLE_wrapper",
               state,
               # MPLE settings
               as.double(to_ergm_Cdouble(elfd)),
               as.integer(.Machine$integer.max), # maxDyads
               as.integer(maxNumDyadTypes),
               PACKAGE="ergm")
    uvals <- z$weightsvector!=0
    zy.e <- z$y[uvals]
    wend.e <- as.numeric(z$weightsvector[uvals])
    xmat.e <- matrix(z$x, ncol=nparam(m,canonical=TRUE), byrow=TRUE)[uvals,,drop=FALSE]
    colnames(xmat.e) <- param_names(m,canonical=TRUE)
    rm(z,uvals)

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
