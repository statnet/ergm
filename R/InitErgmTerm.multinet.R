#' A multinetwork network representation.
#'
#' A function for specifying the LHS of a multi-network (a.k.a. multilevel) ERGM.
#'
#' @param ... network specification, in one of two formats:
#' 
#'   1. An (optionally named) list of networks with same directedness and bipartedness (but possibly different sizes).
#'
#'   1. Several networks as (optionally named) arguments.
#'
#' @return A network object with multinetwork metadata.
#'
#' @seealso [Help on model specification][ergm-terms] for specific terms.
#' 
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- Networks(list(samplk1, samplk2))
#' ergm(monks ~ N(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- Networks(list(samplk1, samplk2))
#' ergm(monks ~ N(~edges))
#' 
#' @export
Networks <- function(...){
  args <- list(...)
  if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multinetwork specification. See help for information.")

  constraintsl <- lapply(nwl, get.network.attribute, "constraints")
  if(!all_identical(lapply(constraintsl, .unenv))) stop("Networks have differing constraint structures. This is not supported at this time.")
  obs.constraintsl <- lapply(nwl, get.network.attribute, "obs.constraints")
  if(!all_identical(lapply(obs.constraintsl, .unenv))) stop("Networks have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, blockID.vattr=".NetworkID", blockName.vattr=".NetworkName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints"), subnet.cache=TRUE)
  nw %n% "constraints" <-
      if(NVL(nwl[[1]]%n%"constraints",~.)==~.)
        ~blockdiag(".NetworkID")
      else
        append.rhs.formula(nwl[[1]]%n%"constraints", list(call("blockdiag",".NetworkID")), TRUE)
  if("obs.constraints" %in% list.network.attributes(nwl[[1]])) nw %n% "obs.constraints" <- nwl[[1]]%n%"obs.constraints"

  nw
}

## InitErgmTerm..layer.nets <- function(nw, arglist, response=NULL, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c(),
##                       vartypes = c(),
##                       defaultvalues = list(),
##                       required = c())
##   list(name="_layer_nets", coef.names=c(), inputs=unlist(.layer_vertexmap(nw)), dependence=FALSE)
## }

InitErgmTerm..subnets <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  list(name="_subnets", coef.names=c(), inputs=c(unlist(.block_vertexmap(nw, a$attrname))), dependence=FALSE)
}

InitErgmTerm.N <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula","weight","wname"),
                      vartypes = c("formula","character,numeric,function","character"),
                      defaultvalues = list(NULL,NULL,NULL),
                      required = c(TRUE,FALSE,FALSE))

  f <- a$formula

  nwl <- .split_constr_network(nw, ".NetworkID", ".NetworkName")
  nwnames <- names(nwl)
  nn <- length(nwl) 

  w <- switch(mode(a$weight),
              character = sapply(nwl, get.network.attribute, a$weight),
              numeric = rep(w, length.out=nn),
              `function` = sapply(nwl, a$weight),
              `NULL` = rep(1,nn))
  
  auxiliaries <- ~.subnets(".NetworkID")

  ms <- lapply(nwl, function(nw1){
    f <- nonsimp.update.formula(f, nw1~.)
    m <- ergm.getmodel(f, nw1, response=response,...)
    Clist <- ergm.Cprepare(nw1, m, response=response)
    list(model = m,
         inputs = pack.Clist_as_num(Clist),
         gs = ergm.emptynwstats.model(m))
  })

  inputs <- c(w, unlist(lapply(ms, "[[", "inputs")))

  ## FIXME: The following assumes that the formula evaluated on all the subnetworks produces the same set of statistics (or, at least, one with the same set of names).
  
  gs <- c(sapply(ms, "[[", "gs")%*%w)
  
  c(list(name="MultiNet", coef.names = paste0(NVL(a$wname,"sum"),"*",ms[[1]]$model$coef.names), inputs=inputs, dependence=!is.dyad.independent(ms[[1]]$model), emptynwstats = gs, auxiliaries = auxiliaries),
    passthrough.curved.ergm_model(ms[[1]]$m, function(x) paste0(NVL(a$wname,"sum"),"*",x)))
}
