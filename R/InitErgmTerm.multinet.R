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
        append_rhs.formula(nwl[[1]]%n%"constraints", list(call("blockdiag",".NetworkID")), TRUE)
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
                      vartypes = c("formula","character,numeric,function,formula","character"),
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
    f <- nonsimp_update.formula(f, nw1~.)
    m <- ergm.getmodel(f, nw1, response=response,...)
    Clist <- ergm.Cprepare(nw1, m, response=response)
    list(model = m,
         inputs = pack.Clist_as_num(Clist),
         gs = ergm.emptynwstats.model(m))
  })

  nparams <- sapply(lapply(ms, `[[`, "model"), nparam, canonical=FALSE)
  nstats <- sapply(lapply(ms, `[[`, "model"), nparam, canonical=TRUE)
  
  if(is(a$weight,"formula")){
    if(!all_identical(nparams)) stop("N() operator with linear model weights only supports models with the same numbers of parameters for every network. This may change in the future.")
    nparam <- nparams[1]
    
    nattrs <- Reduce(union, lapply(nwl, list.network.attributes))
    nattrs <- as.data.frame(lapply(nattrs, function(nattr) sapply(lapply(nwl, get.network.attribute, nattr), function(x) if(is.null(x) || length(x)!=1) NA else x)), col.names=nattrs)
    xm <- do.call(model.matrix, list(a$weight, nattrs), envir=environment(a$weight)) # Each network is a row.
    xl <- lapply(split(xm, row(xm)), rbind) # Rows of xl as a list.
    # What the following does: creates a nn-list of singleton lists
    # containing that network's row of xl, replicates it nparam times,
    # then binds them into nn-list of block-diagonal matrices. These
    # can now be used as covariates to vectorized MANOVA-style
    # parameters.
    Xl <- lapply(lapply(mapply(rep, lapply(xl,list), nparam, SIMPLIFY=FALSE), Matrix::.bdiag), as.matrix)

    inputs <- c(c(0,cumsum(nstats)[-length(nstats)]), unlist(lapply(ms, `[[`, "inputs")))

    map <- function(x, n, ...){
      unlist(mapply(ergm.eta, lapply(lapply(Xl, `%*%`, c(x)),c), lapply(lapply(ms, `[[`, "model"), `[[`, "etamap"), SIMPLIFY=FALSE))
    }
    gradient <- function(x, n, ...){
      do.call(cbind,mapply(crossprod, Xl, mapply(ergm.etagrad, lapply(lapply(Xl, `%*%`, c(x)),c), lapply(lapply(ms, `[[`, "model"), `[[`, "etamap"), SIMPLIFY=FALSE), SIMPLIFY=FALSE))
    }
    params <- rep(list(NULL), nparam*ncol(xm))
    names(params) <- paste0('N(',rep(colnames(xm), nparam),')*',rep(param_names(ms[[1]]$model, canonical=FALSE), each=ncol(xm)))
    coef.names <- paste0('N#',rep(seq_len(nn), nstats),'*',unlist(lapply(lapply(ms, `[[`, "model"), param_names, canonical=TRUE)))
    gs <- unlist(lapply(ms, `[[`, "gs"))

    list(name="MultiNets", coef.names = coef.names, inputs=inputs, dependence=!all(sapply(lapply(ms, `[[`, "model"), is.dyad.independent)), emptynwstats = gs, auxiliaries = auxiliaries, map = map, gradient = gradient, params = params)
  }else{
    if(!all_identical(nparams) || !all_identical(nstats)) stop("N() operator only supports models with the same numbers of parameters and coefficients for every network. This may change in the future.")
    
    inputs <- c(w, unlist(lapply(ms, `[[`, "inputs")))
    gs <- c(sapply(ms, `[[`, "gs")%*%w)
    
    c(list(name="MultiNet", coef.names = paste0("N(",NVL(a$wname,"sum"),")*",ms[[1]]$model$coef.names), inputs=inputs, dependence=!is.dyad.independent(ms[[1]]$model), emptynwstats = gs, auxiliaries = auxiliaries),
      passthrough.curved.ergm_model(ms[[1]]$m, function(x) paste0('N(',NVL(a$wname,"sum"),")*",x)))
  }
}
