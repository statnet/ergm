#  File R/ergm.getnetwork.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
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

#' @rdname ergm.getmodel
#' @description The \code{ergm.getnetwork} function ensures that the
#'   network in a given formula is valid; if so, the network is
#'   returned; if not, execution is halted with warnings
#' @return \code{ergm.getnetwork} returns: \itemize{ \item the network
#'   from the formula IF (i) the formula was correctly structured and
#'   (ii) the network is found within the formula's enviornment }
#' @export ergm.getnetwork
ergm.getnetwork <- function (form, loopswarning=TRUE){
  nw <- ERRVL(
    try({
      tmp <- eval_LHS.formula(form)
      if(is.network(tmp)) tmp else as.network(tmp)
    }, silent = TRUE),
    stop("Invalid network on the LHS of the formula.")
  )

  if (loopswarning) {
    e <- as.edgelist(nw)
    if(any(e[,1]==e[,2])) {
      print("Warning:  This network contains loops")
    } else if (has.loops(nw)) {
      print("Warning:  This network is allowed to contain loops")
    }
  }
  nw
}
