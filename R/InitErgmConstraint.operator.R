#  File R/InitErgmConstraint.operator.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

#' @templateVar name Dyads
#' @title Constrain fixed or varying dyad-independent terms
#' @description This is an "operator" constraint that takes one or two [`ergmTerm`] dyad-independent formulas. For the terms in the `vary=` formula, only those that change at least one of the terms will be allowed to vary, and all others will be fixed. If both formulas are given, the dyads that vary either for one or for the other will be allowed to vary. Note that a formula passed to `Dyads` without an argument name will default to `fix=` .
#'
#' @usage
#' # Dyads(fix=NULL, vary=NULL)
#' @param fix,vary formula with only dyad-independent terms
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept operator
#' @concept directed
#' @concept undirected
InitErgmConstraint.Dyads<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("fix", "vary"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  fix <- a$fix; vary <- a$vary

  if(is.null(fix) & is.null(vary))
    ergm_Init_stop("Dyads constraint takes at least one argument, either ",sQuote("fix")," or ",sQuote("vary")," or both.")

  for(f in c(fix, vary)){
    f[[3]] <- f[[2]]
    f[[2]] <- nw
    if(!is.dyad.independent(f)) ergm_Init_stop("Terms passed to the Dyads constraint must be dyad-independent.")
  }

  list(
    free_dyads = function(){
      fd <- lapply(list(fix=fix,vary=vary),
                   function(f){
                     if(!is.null(f)){
                       f[[3]] <- f[[2]]
                       f[[2]] <- nw
                       m <- ergmMPLE(f, expand.bipartite=TRUE, output="array")$predictor
                       m <- (m != 0) %|% FALSE
                       if(!is.directed(nw)){
                         m <- m | aperm(m, c(2L,1L,3L))
                       }
                       lapply(seq_len(dim(m)[3]), function(i) as.rlebdm(m[,,i]))
                     }
                   })
      fd$fix <- if(length(fd$fix)) fd$fix %>% map(`!`) %>% reduce(`&`)
      fd$vary <- if(length(fd$vary)) fd$vary %>% reduce(`|`)
      fd <- Reduce(`|`, fd)

      compress(fd)
    },
    dependence = FALSE
  )
}
