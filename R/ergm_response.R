#  File R/ergm_response.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' Update the network and the response argument.
#'
#' @param nw a [`network`] object.
#' @template response
#'
#' @details
#' 1. If `response` is `NULL` or `logical`, drop all edge attributes except for `na` and return the network and the response as they are.
#' 2. If `response` is a character vector of length 1, drop all edge attributes in `nw` except for the one corresponding to `response`.
#' 3. If `response` is a formula, construct a name for it and assign to that name (as an edge attribute) the result of evaluating the formula environment; drop all other edge attributes. Return as response the name (possibly with the attribute for the formula attached). If the formula's RHS is of the form a|b use the logicalness of b in Step 4.
#' 4. Test if the resulting edge attribute is of mode [`logical`]. If so set `attr(response,'valued')` to `FALSE`, otherwise to `TRUE`.
#'
#' If both `nw` and `response` are ordinary variables (i.e., not expressions) in the parent frame, `nw` (whatever it was named) is overwritten with the updated network and `response` (again, whatever it was named) is deleted. This is for convenience and for making outdated code that relies on `response` fail immediately rather than introduce subtle bugs. Otherwise, the updated network is returned.
#'
#' @examples
#'
#' preproc_check_print <- function(nw, response){
#'   ergm_preprocess_response(nw, response) 
#'   str(list(
#'        valued = is.valued(nw),
#'        el = head(as.edgelist(nw, attrname=nw%ergmlhs%"response", output="tibble"),3)
#'   ))
#' }
#'
#' data(florentine)
#' preproc_check_print(flomarriage, NULL)
#'
#' flomarriage %e% "w" <- runif(network.edgecount(flomarriage))
#' flomarriage %e% "s" <- rep(c(-1,1), length.out=network.edgecount(flomarriage))
#'
#' # Edge attribute expression
#' preproc_check_print(flomarriage, ~w*s)
#'
#' # Named
#' preproc_check_print(flomarriage, wsprod~w*s)
#'
#' # Binary from valued
#' preproc_check_print(flomarriage, ~s>0)
#' 
#' # Default edge attribute mode is valued
#' flomarriage[,] <- 0 # Empty network
#' preproc_check_print(flomarriage, ~w*s)
#' 
#' # Force default edge attribute mode to binary
#' preproc_check_print(flomarriage, ~w|TRUE)
#' 
#' @keywords internal
#' @export
ergm_preprocess_response <- function(nw, response){
  nwname <- deparse(substitute(nw))
  respname <- deparse(substitute(response))

  inplace <- exists(nwname, parent.frame()) && exists(respname, parent.frame())

  clear_eattrs <- function(nw, except="na"){
    for(a in setdiff(list.edge.attributes(nw), except))
      nw <- delete.edge.attribute(nw, a)
    nw
  }

  if(!is.null(response)) nw %ergmlhs% "response" <- response
  else response <- nw %ergmlhs% "response"

  if(is.null(response) || is.logical(response)){
    nw <- clear_eattrs(nw)
  }else{
    binary_basis <- NULL
    if(is.character(response) && length(response)==1){
      el <- as_tibble(nw, attrnames=response, na.rm=FALSE)
      name <- response
    }else if(is(response,"formula")){
      if(length(ult(response))>1 && ult(response)[[1]]=="|"){ # I.e., RHS is a|b,
        binary_basis <- ult(response)[[3]] # use b as a basis, and
        ult(response) <- ult(response)[[2]] # use a as the expression.
      }

      el <- as_tibble(nw, attrnames=list.edge.attributes(nw), na.rm=FALSE)
      name <- if(length(response)==3) as.character(response[[2]]) else despace(deparse(ult(response), 500L))
      a <- if(nrow(el)) eval(ult(response), envir=c(list(.nw=nw), as.list(el)), enclos=environment(response))
           else c()
      el <- el[1:2]
      el[[name]] <- a
    }

    nw[,] <- 0
    binary <- is.logical(NVL(binary_basis, el[[name]]))

    if(nrow(el)!=0L){
      el <- el[is.na(el[[name]]) | el[[name]],]
      el$eal <- lapply(el[[name]], if(binary) function(a) list(na=is.na(a)) else function(a) list(is.na(a),a))
      nw <- add.edges(nw, tail=el$.tail, head=el$.head, names.eval=rep(list(c("na", if(!binary) name)), nrow(el)), vals.eval=el$eal)
    }

    nw %ergmlhs% "response" <- response <- if(!binary) structure(name, valued=!binary)
  }

  if(inplace){
    assign(nwname, nw, pos=parent.frame())
    rm(list=respname, pos=parent.frame())
    NULL
  }else nw
}
