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
#' @param fix,vary binary ERGM formulas with only dyad-independent terms
#'
#' @note If used in a valued ERGM, this constraint will treat the
#'   `fix` and `vary` formulas as binary formulas. For example, terms
#'   such as \ergmTerm{ergm}{nodematch}{()} will work as expected, but
#'   \ergmTerm{ergm}{sum}{} or \ergmTerm{ergm}{atleast}{()} will
#'   produce an error.
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

  if (is.valued(nw)) {
    ergm_Init_message("Note that constraint ", sQuote("Dyads()"), " will treat the ERGM terms as binary even if the model is valued.")
    nw %ergmlhs% "response" <- NULL
  }

  if (verify_dind)
    for(f in c(fix, vary)){
      if (!is.dyad.independent(f, basis = nw))
        ergm_Init_stop("Terms passed to the ", sQuote("Dyads"), " constraint must be dyad-independent.")
    }

  list(
    free_dyads = function(){
      fd <- lapply(list(fix=fix,vary=vary),
                   function(f){
                     if(!is.null(f)){
                       m <- ergmMPLE(f, basis = nw, expand.bipartite=TRUE, output="array")$predictor
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
#' @note `formula` can also be a [`term_list`] or [`character`]. In
#'   the latter case, if there are multiple strings, they will be
#'   concatenated with `+`, and if they do not start with "~", one
#'   will be prepended. Its environment will be inherited from the
#'   top-level formula.
#'
#' @seealso [base::I()] (a.k.a. `AsIs`)
#'
#' @template ergmConstraint-general
#'
#' @concept operator
InitErgmConstraint.I <- function(nw, arglist, ..., env) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula,character,term_list"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  f <- a$formula

  if (is.character(f)) {
    if (length(f) > 1) f <- paste(f, collapse = " + ")

    if (!startsWith(trimws(f, "left"), "~")) f <- paste0("~", f)

    f <- as.formula(f, env)
  }

  ergm_conlist(a$formula, nw, ...)
}

#' @templateVar name ChangeStats
#' @title Specified statistics must remain constant
#' @description Given an [`ergmTerm`] formula, this constraint
#'   prevents any changes to the network that change its value. Unlike
#'   \ergmConstraint{ergm}{Dyads}{(fix = ...)}, the terms can be
#'   dyad-dependent and are calculated at the same time as the
#'   proposal rather than used to select proposable dyads in the first
#'   place.
#'
#' @usage
#' # ChangeStats(fix, check_dind = TRUE)
#' @param fix an [`ergmTerm`] formula.
#' @param check_dind logical; if `fix` turns out to be dyad-independent, fall back to \ergmConstraint{ergm}{Dyads}{(fix)}.
#'
#' @section Dyadic dependence:
#'
#' If the statistic is dyad-independent, \ergmConstraint{ergm}{Dyads}{(fix)} will usually be faster (at least for estimation) and will make MPLE far more accurate, so if `ChangeStats` detects a dyad-independent constraint, it will fall back to `Dyads`. This can be overridden by setting `check_dind = FALSE`. This fallback is never used for valued ERGMs.
#'
#' @section Sampleability is not guaranteed:
#'
#' It is perfectly possible for this constraint to make sampling
#' impossible. For example, `ChangeStats(~edges)` will prevent any
#' proposals that change the number of edges in the network, but
#' \pkg{ergm} has no way of knowing that two toggles (edge and
#' non-edge) are now required, so it will keep trying making proposals
#' and failing.
#'
#' More insidiously, it is in principle possible for the constraint to
#' split the sample space into two parts such that it is not possible
#' to go between them without passing through a state that breaks the
#' constraint.
#'
#' As of this writing, this constraint's implementation checks only
#' the difference between the current and the proposed network, not
#' any in between. Thus, for example, if the proposal *can* preserve
#' the number of edges, the (redundant) `ChangeStats(~edges)`
#' constraint will allow sampling to occur. An option to check after
#' every change may be added in the future.
#'
#' Thus, it is recommended to use this constraint to, e.g., prevent
#' certain motifs from forming.
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

  # For valued networks, dyad independence doesn't guarantee Dyads()
  # will work.
  if (!is.valued(nw) && a$check_dind) {
    if (is.dyad.independent(a$fix, basis = nw)) {
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

#' Helper function for proposals implementing the \ergmConstraint{ergm}{ChangeStats}{()} constraint
#'
#' This function inserts a `.submodel(terms)` auxiliary among the
#' proposal's other auxiliaries and records its position, where
#' `terms` are the statistics passed to the \ergmConstraint{ergm}{ChangeStats}{()}
#' constraint.
#'
#' @param arguments `arguments` argument of a call to
#'   `Init*ErgmProposal.*()` function.
#' @param aux_before,aux_after either [`formula`]s or a [`term_list`]s
#'   for other auxiliaries of the proposal.
#'
#' @return a list containing two elements: `auxiliaries` with the
#'   concatenated [`term_list`] and `ChangeStat_pos`, giving the
#'   position of the `.submodel(terms)` auxiliary or `NULL` if none
#'   was specified.
#'
#' @keywords internal
#' @export
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
        if (length(terms)) list_rhs.formula(~.submodel(terms)),
        aux_after),
    ChangeStat_pos = if (length(terms)) as.integer(length(aux_before))
  )
}
