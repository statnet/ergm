#  File R/godfather.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
#=========================================================================
# This file contains the following 2 functions for computing changestat
# summaries of dynamic networks ??
#   <ergm.godfather>
#   <control.godfather>
#=========================================================================



#' A function to apply a given series of changes to a network.
#' 
#' Gives the network a series of proposals it can't refuse. Returns the
#' statistics of the network, and, optionally, the final network.
#' 
#' 
#' @param formula An \code{\link{ergm}}-style formula, with a
#'   \code{\link{network}} on its LHS.
#' @param changes Either a matrix with three columns: tail, head, and
#'   new value, describing the changes to be made; or a list of such
#'   matrices to apply these changes in a sequence. For binary network
#'   models, the third column may be omitted. In that case, the
#'   changes are treated as toggles. Note that if a list is passed, it
#'   must either be all of changes or all of toggles.
#' @template response
#' @param end.network Whether to return a network that
#'   results. Defaults to \code{FALSE}.
#' @param stats.start Whether to return the network statistics at
#'   \code{start} (before any changes are applied) as the first row of
#'   the statistics matrix.  Defaults to \code{FALSE}, to produce
#'   output similar to that of \code{\link[=simulate.ergm]{simulate}}
#'   for ERGMs when \code{output="stats"}, where initial network's
#'   statistics are not returned.
#' @param changes.only Whether to return network statistics or only
#'   their changes relative to the initial network.
#'
#' @templateVar mycontrol control.ergm.godfather
#' @template control
#' @template verbose
#'
#' @return If \code{end.network==FALSE} (the default), an
#'   \code{\link{mcmc}} object with the requested network statistics
#'   associed with the network series produced by applying the
#'   specified changes. Its \code{\link{mcmc}} attributes encode the
#'   timing information: so \code{\link{start}(out)} gives the time
#'   point associated with the first row returned, and
#'   \code{\link{end}(out)} out the last. The "thinning interval" is
#'   always 1.
#' 
#' If \code{end.network==TRUE}, return a \code{\link{network}} object,
#' representing the final network, with a matrix of statistics
#' described in the previous paragraph attached to it as an
#' \code{attr}-style attribute \code{"stats"}.
#' @seealso `tergm.godfather()` in \CRANpkg{tergm}, [simulate.ergm()],
#'   [simulate.formula()]
#' @examples
#' data(florentine)
#' ergm.godfather(flomarriage~edges+absdiff("wealth")+triangles,
#'                changes=list(cbind(1:2,2:3),
#'                             cbind(3,5),
#'                             cbind(3,5),
#'                             cbind(1:2,2:3)),
#'                stats.start=TRUE)
#' @export ergm.godfather
ergm.godfather <- function(formula, changes=NULL, response=NULL,
                           end.network=FALSE,
                           stats.start=FALSE,
                           changes.only=FALSE,
                           verbose=FALSE,
                           control=control.ergm.godfather()){
  on.exit(ergm_Cstate_clear())

  check.control.class("ergm.godfather", "ergm.godfather")

  if(!is.list(changes)) changes <- list(changes)

  nw <- ergm.getnetwork(formula)
  ergm_preprocess_response(nw,response)

  ncols <- sapply(changes, ncol)
  if(!all_identical(ncols) || ncols[1]<2 || ncols[1]>3 || (is.valued(nw)&&ncols[1]==2)) abort("Invalid format for list of changes. See help('ergm.godfather').")

  m <- ergm_model(formula, nw, term.options=control$term.options)
  state <- ergm_state(nw, model=m)
  state <- update(state, stats = if(changes.only) numeric(nparam(state,canonical=TRUE)) else summary(state))

  changem <- changes %>% map(~rbind(0L,.)) %>% do.call(rbind, .) # 0s are sentinels indicating next iteration.
  if(!stats.start) changem <- changem[-1,,drop=FALSE] # I.e., overwrite the initial statistic rather than advancing past it first thing.
  if(!is.directed(nw)) {
    tails <- changem[,1]
    heads <- changem[,2]
    changem[,1] <- pmin(tails, heads)
    changem[,2] <- pmax(tails, heads)
  }
  
  if(verbose) message_print("Applying changes...\n")
  z <-
    if(!is.valued(state))
      .Call("Godfather_wrapper",
            state,
            # Godfather settings
            as.integer(changem[,1]),
            as.integer(changem[,2]),
            if(ncol(changem)==3) as.integer(changem[,3]) else integer(0),
            as.logical(end.network),
            as.integer(verbose),
            PACKAGE="ergm")
    else
      .Call("WtGodfather_wrapper",
            state,
            # Godfather settings
            as.integer(changem[,1]),
            as.integer(changem[,2]),
            as.double(changem[,3]),
            as.logical(end.network),
            as.integer(verbose),
            PACKAGE="ergm")

  stats <- matrix(z$s, ncol=nparam(m,canonical=TRUE), byrow=TRUE)
  colnames(stats) <- param_names(m, canonical=TRUE)

  #' @importFrom coda mcmc
  stats <- mcmc(stats)
  
  if(end.network){ 
    if(verbose) cat("Creating new network...\n")
    newnetwork <- as.network(update(z$state))
    attr(newnetwork,"stats")<-stats
    newnetwork
  }else stats
}

#' Control parameters for [ergm.godfather()].
#'
#' Returns a list of its arguments.
#'
#' @template term_options
#' 
#' @export control.ergm.godfather
control.ergm.godfather<-function(term.options=NULL){
  control <- handle.controls("control.ergm.godfather")
  set.control.class("control.ergm.godfather")
}
