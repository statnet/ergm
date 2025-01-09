#  File R/godfather.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
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
#' @param object An [ergm()]-style formula, with a [`network`] on its
#'   LHS, an [ergm_model()] or the object appropriate to the method.
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
#' @param control Deprecated; arguments such as `term.options` can be
#'   passed directly.
#' @param formula Deprecated; replaced with `object` for consistency.
#' @template verbose
#'
#' @param ... additional arguments to [ergm_model()].
#'
#' @template basis
#'
#' @return If \code{end.network==FALSE} (the default), an
#'   [`mcmc`] object with the requested network statistics
#'   associed with the network series produced by applying the
#'   specified changes. Its [`mcmc`] attributes encode the
#'   timing information: so \code{\link{start}(out)} gives the time
#'   point associated with the first row returned, and
#'   \code{\link{end}(out)} out the last. The "thinning interval" is
#'   always 1.
#' 
#' If \code{end.network==TRUE}, return a [`network`] object,
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
ergm.godfather <- function(object, changes=NULL,
                           ...,
                           end.network=FALSE,
                           stats.start=FALSE,
                           changes.only=FALSE,
                           verbose=FALSE,
                           basis = NULL,
                           formula = NULL){
  ## TODO: Remove around ergm 4.9:
  if(!is.null(formula)){
    .Deprecate_once(msg = paste0("Argument ", sQuote("formula="), " to ", sQuote("ergm.godfather()"), "has been replaced with ", sQuote("formula="), "."))
    object <- formula
  }

  UseMethod("ergm.godfather", object)
}

#' @rdname ergm.godfather
#' @export
ergm.godfather.formula <- function(object, changes=NULL, response=NULL,
                                   ...,
                                   end.network=FALSE,
                                   stats.start=FALSE,
                                   changes.only=FALSE,
                                   verbose=FALSE,
                                   control=NULL,
                                   basis = ergm.getnetwork(object)){
  ergm_preprocess_response(basis,response)
  ## TODO: Remove this workaround around the 4.9 release.
  m <- if("term.options" %in% ...names()) ergm_model(object, basis, ...)
       else ergm_model(object, basis, term.options = control$term.options, ...)

  ergm.godfather(m, changes = changes, ...,
                 end.network = end.network,
                 stats.start = stats.start,
                 changes.only = changes.only,
                 verbose = verbose,
                 control = control,
                 basis = basis)
}

#' @rdname ergm.godfather
#' @note [ergm.godfather.ergm_model()] is a lower-level interface, providing
#'   an [ergm.godfather()] method for the [`ergm_model`] class. The `basis`
#'   argument is required.
#' @export
ergm.godfather.ergm_model <- function(object, changes=NULL,
                                   ...,
                                   end.network=FALSE,
                                   stats.start=FALSE,
                                   changes.only=FALSE,
                                   verbose=FALSE,
                                   control=NULL,
                                   basis = NULL){
  if(is.null(basis)) stop("This method requires the ", sQuote("basis="), " argument.")

  state <- ergm_state(basis, model=object)
  state <- update(state, stats = if(changes.only) numeric(nparam(state,canonical=TRUE)) else summary(state))

  s <- ergm.godfather(state, changes = changes, ...,
                      end.network = end.network, stats.start = stats.start,
                      changes.only = changes.only, verbose = verbose,
                      control = control)

  if(end.network) structure(as.network(s), stats = attr(s, "stats"))
  else s
}

#' @rdname ergm.godfather
#' @note [ergm.godfather.ergm_model()] is a lower-level interface, providing
#'   an [ergm.godfather()] method for the [`ergm_model`] class. The `basis`
#'   argument is required.
#' @export
ergm.godfather.ergm_state <- function(object, changes=NULL,
                                   ...,
                                   end.network=FALSE,
                                   stats.start=FALSE,
                                   verbose=FALSE,
                                   control=NULL){
  if(!is.null(control)) check.control.class("ergm.godfather", "ergm.godfather")

  if(!is.list(changes)) changes <- list(changes)

  ncols <- sapply(changes, ncol)
  if(!all_identical(ncols) || ncols[1]<2 || ncols[1]>3 || (is.valued(object)&&ncols[1]==2)) abort("Invalid format for list of changes. See help('ergm.godfather').")

  changem <- changes %>% map(~rbind(0L,.)) %>% do.call(rbind, .) # 0s are sentinels indicating next iteration.
  if(!stats.start) changem <- changem[-1,,drop=FALSE] # I.e., overwrite the initial statistic rather than advancing past it first thing.
  if(!is.directed(object$nw0)) {
    tails <- changem[,1]
    heads <- changem[,2]
    changem[,1] <- pmin(tails, heads)
    changem[,2] <- pmax(tails, heads)
  }
  
  if(verbose) message("Applying changes...")
  on.exit(ergm_Cstate_clear())
  z <-
    if(!is.valued(object))
      .Call("Godfather_wrapper",
            object,
            # Godfather settings
            as.integer(changem[,1]),
            as.integer(changem[,2]),
            if(ncol(changem)==3) as.integer(changem[,3]) else integer(0),
            as.logical(end.network),
            as.integer(verbose),
            PACKAGE="ergm")
    else
      .Call("WtGodfather_wrapper",
            object,
            # Godfather settings
            as.integer(changem[,1]),
            as.integer(changem[,2]),
            as.double(changem[,3]),
            as.logical(end.network),
            as.integer(verbose),
            PACKAGE="ergm")

  stats <- matrix(z$s, ncol=nparam(object, canonical=TRUE), byrow=TRUE)
  colnames(stats) <- param_names(object, canonical=TRUE)

  #' @importFrom coda mcmc
  stats <- mcmc(stats)
  
  if(end.network){ 
    if(verbose) message("Creating new network...")
    newnetwork <- update(z$state)
    attr(newnetwork,"stats") <- stats
    newnetwork
  }else stats
}
