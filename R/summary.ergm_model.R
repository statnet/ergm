#  File R/summary.ergm_model.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

#' Evaluate network summary statistics from an initialized ergm model
#' 
#' Returns a vector of the model's statistics for a given network or
#' an empty network. This is a low-level function that should not be
#' used by end-users, but may be useful to developers.
#'
#' @param object an [`ergm_model`] object.
#' @param nw a [`network`] whose statistics are to be evaluated. If
#'   `NULL`, returns empty network's statistics for that model.
#' @template response
#' @template dotdotdot
#' 
#' @seealso [summary_formula()]
#' @keywords internal
#' @export
summary.ergm_model <- function(object, nw=NULL, response=NULL,...){
  m <- object

  # Adjust to global values. This needs to happen before the C call,
  # so that an s_function, if exists, could override.
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgmTerm function
  # Read the comments at the top of InitErgm.R or InitErgmTerm.R for 
  # an explanation of the $emptynwstats mechanism
  gs <- numeric(nparam(m,canonical=TRUE))
  if(length(gs)==0) return(gs) # Escape if the model has 0 statistics.
  
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$terms[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }

  # If no actual network, we are done.
  if(is.null(nw)) return(gs)
  
  # Note that the empty network statistics are passed to the C
  # code. The reason is that if an s_??? function exists, it can
  # overwrite them, since it can compute the whole thing, while if
  # only the d_??? function exists, it needs to add on to empty
  # network statistics.
  
  Clist <- ergm.Cprepare(nw, m, response=response)
  
  # *** don't forget, tails are passes in first now, notheads  
  gs <- if(is.null(response))
         .C("network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(!is.null(Clist$lasttoggle)), as.integer(Clist$lasttoggle),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
         else
         .C("wt_network_stats_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads), as.double(Clist$weights), as.integer(!is.null(Clist$time)), as.integer(Clist$time), as.integer(!is.null(Clist$lasttoggle)), as.integer(Clist$lasttoggle),
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = as.double(gs),
            PACKAGE="ergm"
            )$gs
  names(gs) <- param_names(m,canonical=TRUE)

  gs
}


