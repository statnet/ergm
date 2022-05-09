#  File R/InitWtErgmTerm.operator.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################
InitWtErgmTerm.Passthrough <- function(nw, arglist, ...){
  out <- InitErgmTerm.Passthrough(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

#' @templateVar name B
#' @title Wrap binary terms for use in valued models
#' @description Wraps binary `ergm` terms for use in valued models, with `formula` specifying which terms
#'   are to be wrapped and `form` specifying how they are to be
#'   used and how the binary network they are evaluated on is to be constructed. 
#'   
#' @details For example, `B(~nodecov("a"), form="sum")` is equivalent to
#'   `nodecov("a", form="sum")` and similarly with
#'   `form="nonzero"` .
#'   
#'   When a valued implementation is available, it should be
#'   preferred, as it is likely to be faster.
#'
#' @usage
#' # valued: B(formula, form)
#' @param formula a one-sided [ergm()]-style formula whose RHS contains the
#'   binary ergm terms to be evaluated. Which terms may be used
#'   depends on the argument `form`
#' @param form One of three values:
#'   - `"sum"`: see section "Generalizations of
#'   binary terms" in [`ergmTerm`] help; all terms in `formula` must
#'   be dyad-independent.
#'   - `"nonzero"`: section "Generalizations of
#'   binary terms" in [`ergmTerm`] help; any binary `ergm` terms
#'   may be used in `formula` .
#'   - a one-sided formula value-dependent
#'   network. `form` must contain one "valued" `ergm` term, with
#'   the following properties:
#'     - dyadic independence;
#'     - dyadwise contribution of either 0 or 1; and
#'     - dyadwise contribution of 0 for a 0-valued dyad.
#'   
#'     Formally, this means that it is expressable as
#'     \deqn{g(y) = \sum_{i,j} f_{i,j}(y_{i,j}),}{sum[i,j] f[i,j] (y[i,j]),}
#'     where for all \eqn{i}, \eqn{j}, and \eqn{y},
#'     \eqn{f_{i,j}(y_{i,j})} is either 0 or 1 and, in particular,
#'     \eqn{f_{i,j}(0)=0}{f[i,j] (0)=0}.
#'   
#'     Examples of such terms include `nonzero` ,
#'     `ininterval()` , `atleast()` , `atmost()` ,
#'     `greaterthan()` , `lessthen()` , and `equalto()` .
#'   
#'     Then, the value of the statistic will be the value of the
#'     statistics in `formula` evaluated on a binary network that is
#'     defined to have an edge if and only if the corresponding
#'     dyad of the valued network adds 1 to the valued term in
#'     `form` .
#'
#' @template ergmTerm-general
#'
#' @concept operator
InitWtErgmTerm.B <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "form"),
                      vartypes = c("formula", "character,formula"),
                      defaultvalues = list(NULL, "sum"),
                      required = c(TRUE, FALSE))
  form <- if(is.character(a$form)) match.arg(a$form,c("sum","nonzero"))
          else a$form

  nwb <- nw
  nwb %ergmlhs% "response" <- NULL
  m <- ergm_model(a$formula, nwb, ..., offset.decorate=FALSE)
  ergm_no_ext.encode(m)

  if(!is.dyad.independent(m) && form=="sum") stop("Only dyad-independent binary terms can be imported with form 'sum'.")
  
  if(is(form, "formula")){
    form.name <- despace(deparse(ult(form)))
    name <- "import_binary_term_form"
    auxiliaries <- trim_env(~.binary.formula.net(form),"form")
  }else{
    form.name <- form
    name <- paste("import_binary_term",form,sep="_")
    auxiliaries <- if(form=="nonzero") trim_env(~.binary.nonzero.net)
  }

  mw <- wrap.ergm_model(m, nwb, ergm_mk_std_op_namewrap('B', form.name))
  if(form=="sum") mw$emptynwstats <- NULL

  c(list(name=name,
         submodel = m,
         auxiliaries=auxiliaries),
    mw)
}

InitWtErgmTerm..binary.nonzero.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_binary_nonzero_net", depenence=FALSE)
}

InitWtErgmTerm..binary.formula.net <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  if(!is.dyad.independent(m) || nparam(m)!=1) stop("The binary test formula must be dyad-independent and have exactly one statistc.")

  nw[,] <- FALSE
  gs <- summary(m, nw)
  if(gs!=0) stop("At this time, the binary test term must have the property that its dyadwise components are 0 for 0-valued relations. This limitation may be removed in the future.")
  
  c(list(name="_binary_formula_net", submodel=m, depenence=FALSE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, NULL))
}

# Arguments and outputs are identical to the binary version, except for the C routine names.
#' @rdname Sum-operator-ergmTerm
#' @usage
#' # valued: Sum(formulas, label)
InitWtErgmTerm.Sum <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.Sum
  term <- f(...)
  term$name <- "wtSum"
  term
}

#' @rdname Label-ergmTerm
#' @usage
#' # valued: Label(formula, label, pos)
InitWtErgmTerm.Label <- function(nw, arglist, ...){
  out <- InitErgmTerm.Label(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

#' @rdname Curve-ergmTerm
#' @usage
#' # valued: Curve(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf, cov=NULL)
InitWtErgmTerm.Curve <- function(nw, arglist, ...){
  out <- InitErgmTerm.Curve(nw, arglist, ...)
  out$name <- "wtpassthrough_term"
  out
}

InitWtErgmTerm..submodel_and_summary <- function(nw, arglist, ...){
  out <- InitErgmTerm..submodel_and_summary(nw, arglist, ...)
  out$name <- "_wtsubmodel_and_summary_term"
  out
}

#' @rdname Curve-ergmTerm
#' @usage
#' # valued: Parametrise(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf,
#' #           cov=NULL)
InitWtErgmTerm.Parametrise <- InitWtErgmTerm.Curve

#' @rdname Curve-ergmTerm
#' @usage
#' # valued: Parametrize(formula, params, map, gradient=NULL, minpar=-Inf, maxpar=+Inf,
#' #           cov=NULL)
InitWtErgmTerm.Parametrize <- InitWtErgmTerm.Curve

#' @rdname Exp-ergmTerm
#' @usage
#' # valued: Exp(formula)
InitWtErgmTerm.Exp <- function(nw, arglist, ...){
  out <- InitErgmTerm.Exp(nw, arglist, ...)
  out$name <- "wtExp"
  out
}

#' @rdname Log-ergmTerm
#' @usage
#' # valued: Log(formula, log0=-1/sqrt(.Machine$double.eps))
InitWtErgmTerm.Log <- function(nw, arglist, ...){
  out <- InitErgmTerm.Log(nw, arglist, ...)
  out$name <- "wtLog"
  out
}

#' @rdname Prod-ergmTerm
#' @usage
#' # valued: Prod(formulas, label)
InitWtErgmTerm.Prod <- InitErgmTerm.Prod
