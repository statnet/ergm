#  File R/InitErgmConstraint.operator.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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
InitErgmConstraint.Dyads<-function(nw, arglist, ..., verify_dind = TRUE){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("fix", "vary"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  fix <- a$fix; vary <- a$vary

  if(is.null(fix) & is.null(vary))
    ergm_Init_stop(sQuote("Dyads()"), " constraint takes at least one argument, either ",sQuote("fix"), " or ", sQuote("vary"), " or both.")

  if (verify_dind)
    for(f in c(fix, vary)){
      f[[3]] <- f[[2]]
      f[[2]] <- nw
      if(!is.dyad.independent(f)) ergm_Init_stop("Terms passed to the ", sQuote("Dyads"), " constraint must be dyad-independent.")
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


#' @templateVar name I
#' @title Substitute a formula into the constraints formula
#' @description This is a convenience operator that can be used to
#'   paste terms constructed elsewhere into a formula.
#'
#' @usage
#' # I(formula)
#' @param formula a constraints formula
#'
#' @note `formula` can also be a [term_list].
#'
#' @seealso [base::I()] (a.k.a. `AsIs`)
#'
#' @template ergmConstraint-general
#'
#' @concept operator
InitErgmConstraint.I <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  ergm_conlist(a$formula, nw, ...)
}

#' @templateVar name ChangeStats
#' @title Specified statistics must remain constant
#' @description This is an constraint operator that takes a [`ergmTerm`] formula, and prevents any changes to the network that modify its value. Unlike \ergmConstraint{ergm}{Dyads}{(fix = ...)}, the terms can be dyad-dependent and are calculated at the same time as the proposal rather than used to select proposable-dyads in the first place.
#'
#' @usage
#' # ChangeStats(fix, check_dind = TRUE)
#' @param fix an [`ergmTerm`] formula.
#' @param check_dind logical; if `fix` turns out to be dyad-independent, fall back to \ergmConstraint{ergm}{Dyads}{(fix)}.
#'
#' @section Dyadic dependence:
#'
#' If the constraint is dyad-independent, \ergmConstraint{ergm}{Dyads}{(fix)} will usually be faster (at least for estimation) and will make MPLE far more accurate, so if `ChangeStats` detects a dyad-independent constraint, it will fall back to `Dyads`. This can be overridden by setting `check_dind = FALSE`.
#'
#' @section Sampleability is not guaranteed:
#'
#' It is perfectly possible for this constraint to make sampling impossible. For example, `ChangeStats(~edges)` will prevent any proposals that change the number of edges in the network, but \pkg{ergm} has no way of knowing that two toggles (edge and non-edge) are now required, so it will keep trying making proposals and failing.
#'
#' More insidiously, it is in principle possible for the constraint to split the sample space into two parts such that it is not possible to go between them without passing through a state that breaks the constraint.
#'
#' Thus, it is recommended to use this constraint to, e.g., prevent certain motifs from forming.
#'
#'
#' @template ergmConstraint-general
#'
#' @concept operator
#' @concept directed
#' @concept undirected
InitErgmConstraint.ChangeStats <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("fix", "check_dind"),
                      vartypes = c("formula", "logical"),
                      defaultvalues = list(NULL, TRUE),
                      required = c(TRUE, FALSE))

  if(a$check_dind) {
    f <- a$fix
    f[[3]] <- f[[2]]
    f[[2]] <- nw
    if(is.dyad.independent(f)) {
      message(sQuote("ChangeStats()"), " constraint formula is dyad-independent; falling back to ", sQuote("Dyads()"), ".")
      return(InitErgmConstraint.Dyads(nw, arglist, ..., verify_dind = FALSE))
    }
  }


  list(
    fix = a$fix,
    dependence = TRUE,
    constrain = "changestats"
  )
}

EMPTY_TERM_LIST <- term_list(~.)[-1]

ergm_constrain_changestats <- function(arguments,
                                       aux_before = EMPTY_TERM_LIST,
                                       aux_after = EMPTY_TERM_LIST) {
  conlist <- arguments$constraints
  terms <- if (hasName(conlist, "changestats"))
             conlist[names(conlist) == "changestats"] |>
               map("fix") |>
               map(list_rhs.formula) |>
               unname() |>
               do.call(c, args = _)
           else EMPTY_TERM_LIST

  if (is(aux_before, "formula")) aux_before <- list_rhs.formula(aux_before)
  if (is(aux_after, "formula")) aux_after <- list_rhs.formula(aux_after)

  list(
    auxiliaries =
      c(aux_before,
        list_rhs.formula(~.submodel(terms)),
        aux_after),
    ChangeStat_pos = length(aux_before)
  )
}
