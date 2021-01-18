#' @param response Either a character string, a formula, or `NULL` (the default), to specify the response attributes and whether the ERGM is binary or valued. Interpreted as follows: \describe{
#'
#' \item{`NULL`}{Model simple presence or absence, via a binary ERGM.}
#' 
#' \item{character string}{The name of the edge attribute whose value is to be modeled. Type of ERGM will be determined by whether the attribute is [`logical`] (`TRUE`/`FALSE`) for binary or [`numeric`] for valued.}
#' 
#' \item{a formula}{must be of the form `NAME~EXPR|TYPE` (with `|` being literal). `EXPR` is evaluated in the formula's environment with the network's edge attributes accessible as variables. The optional `NAME` specifies the name of the edge attribute into which the results should be stored, with the default being a concise version of `EXPR`. Normally, the type of ERGM is determined by whether the result of evaluating `EXPR` is logical or numeric, but the optional `TYPE` can be used to override by specifying a scalar of the type involved (e.g., `TRUE` for binary and `1` for valued).}
#' 
#' }
