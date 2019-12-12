#  File R/ergm.design.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
#################################################################################
# The <ergm.design> function functions as <ergm.Cprepare> would, but acts on the
# network of missing edges
#
# --PARAMETERS--
#   nw     : the network
#   model  : the model, as returned by <ergm_model>
#   verbose: whether the design matrix should be printed (T or F); default=FALSE
#
# --RETURNED--
#   fd
#      if 'nw' has missing edges, see the return list, 'Clist', from the
#                                 <ergm.Cprepare> function header
#      if 'nw' hasn't any missing edges, the list will only contain NULL
#                                 values for the 'tails' and 'heads' components,
#                                 a 0 for 'nedges' and 'dir' appropriately set
#
# The <ergm.Cprepare.miss> function functions constructs a static edgelist for
# the proposals that need it.
#
# --PARAMETERS--
#   nw     : the network
#
# --RETURNED--
#      a vector of length 1+Nmissing*2. The first element is the number of
#      missing edges, and the remainder a column-major edgelist
################################################################################

#' Obtain the set of informative dyads based on the network structure.
#'
#' Note that this function is not recommended for general use, since
#' it only supports only one way of specifying observational
#' structure---through `NA` edges. It is likely to be deprecated in
#' the future.
#'
#' @param nw a [`network`] object.
#' @param verbose verbosity; currently ignored.
#'
#' @return \code{ergm.design} returns a \code{\link{rlebdm}} of
#'   informative (non-missing, non fixed) dyads.
#' @export ergm.design
ergm.design <- function(nw, verbose=FALSE){
  if(!missing(verbose)) .Deprecate_once(msg="verbose= argument to ergm.design() has been deprecated and is ignored.")
  basecon <- ergm_conlist(~.attributes, nw)
  misscon <- if(!is.ergm_state(nw) && network.naedgecount(nw)) ergm_conlist(~.attributes+observed, nw)
  as.rlebdm(basecon, misscon, which="informative")
}
