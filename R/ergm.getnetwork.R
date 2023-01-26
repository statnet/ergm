#  File R/ergm.getnetwork.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
#################################################################################
# The <ergm.getnetwork> function ensures that the network in a given formula
# is valid; if so, the network is returned; if not, execution is halted with
# warnings
#
# --PARAMETERS--
#   formula     :  the formula as 'network ~ model.term(s)'
#   loopswarning:  whether warnings about loops should be printed (T or F);
#                  default=TRUE
#
# --RETURNED--
#   nw: the network from the formula IF (i) the formula was correctly structured
#       and (ii) the network is found within the formula's enviornment
#
###################################################################################

#' Acquire and verify the network from the LHS of an `ergm` formula
#' and verify that it is a valid network.
#' 
#' The function function ensures that the network in a given formula
#' is valid; if so, the network is returned; if not, execution is
#' halted with warnings.
#'
#' @param formula a two-sided formula whose LHS is a [`network`], an object that can be coerced to a [`network`], or an expression that evaluates to one.
#' @param loopswarning whether warnings about loops should be printed
#'   (`TRUE` or `FALSE`); defaults to `TRUE`.
#' 
#' @return A [`network`] object constructed by evaluating the LHS of
#'   the model formula in the formula's environment.
#' @export ergm.getnetwork
ergm.getnetwork <- function (formula, loopswarning=TRUE){
  nw <- eval_lhs.formula(formula)
  nw <- ensure_network(nw)

  if (loopswarning) {
    e <- as.edgelist(nw)
    if(any(e[,1]==e[,2])) {
      print("Warning:  This network contains loops")
    } else if (has.loops(as.network(nw,populate=FALSE))) {
      print("Warning:  This network is allowed to contain loops")
    }
  }
  nw
}

ensure_network <- function(nw){
  if(!is.network(nw) && !is.ergm_state(nw)){
    nw <- ERRVL(
      try(as.network(nw)),
      abort("A network object on the LHS of the formula or as a basis argument must be given")
    )
  }
  nw
}
