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
                      varnames = c("formula","lm","subset","weights","contrasts","offset","label"),
                      vartypes = c("formula","formula","formula,logical,numeric,expression,call","formula,logical,numeric,expression,call","list","formula,logical,numeric,expression,call","character"),
                      defaultvalues = list(NULL,~1,TRUE,1,NULL,NULL,NULL),
                      required = c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE))

  f <- a$formula

  nwl <- .split_constr_network(nw, ".NetworkID", ".NetworkName")
  nwnames <- names(nwl)
  nn <- length(nwl)
  
  auxiliaries <- ~.subnets(".NetworkID")

  nattrs <- Reduce(union, lapply(nwl, list.network.attributes))
  nattrs <- as.data.frame(lapply(nattrs, function(nattr) sapply(lapply(nwl, get.network.attribute, nattr), function(x) if(is.null(x) || length(x)!=1) NA else x)), col.names=nattrs)

  subset <-
    if(mode(a$subset) %in% c("expression", "call")) eval(if(is(a$subset, "formula")) a$subset[[2]] else a$subset, envir = nattrs, enclos = environment(a$lm))
    else a$subset
  subset <- unwhich(switch(mode(subset),
                           logical = which(rep(subset, length.out = nn)),
                           numeric = subset),
                    nn)

  weights <- if(mode(a$weights) %in% c("expression", "call")) eval(if(is(a$weights, "formula")) a$weights[[2]] else a$weights, envir = nattrs, enclos = environment(a$lm))
             else a$weights
  weights <- rep(weights, length.out=nn)

  offset <- if(mode(a$offset) %in% c("expression", "call")) eval(if(is(a$offset, "formula")) a$offset[[2]] else a$offset, envir = nattrs, enclos = environment(a$lm))
             else a$offset

  subset[weights==0] <- FALSE
  
  ms <- lapply(nwl[subset], function(nw1){
    f <- nonsimp_update.formula(f, nw1~.)
    m <- ergm_model(f, nw1, response=response,...)
    list(model = m,
         inputs = to_ergm_Cdouble(m),
         gs = summary(m))
  })

  nm <- sum(subset)

  #' @importFrom purrr map_int
  nparams <- ms %>% map("model") %>% map_int(nparam, canonical=FALSE)
  nstats <-  ms %>% map("model") %>% map_int(nparam, canonical=TRUE)

  ### Linear model for parameters
  if(!all_identical(nparams)) stop("N() operator with linear model weights only supports models with the same numbers of parameters for every network. This may change in the future.")
  nparam <- nparams[1]

  # model.frame.
  xf <- do.call(stats::lm, list(a$lm, data=nattrs,contrast.arg=a$contrasts, offset = offset, subset=subset, weights=weights, na.action=na.fail, method="model.frame"), envir=environment(a$lm))
  offset <- model.offset(xf) %>% NVL(0) %>% matrix(nrow=nm, ncol=nparam) # offset is actually an nm*q matrix.
  weights <- model.weights(xf)
  if(!all(weights==1)) ergm_Init_abort("Network-level weights different from 1 are not supported at this time.")
  xm <- model.matrix(attr(xf, "terms"), xf, contrasts=a$contrasts)

  Xl <- xm %>%
    split(., row(.)) %>% # Rows of xl as a list.
    map(rbind) %>% # List of row vectors.
    map(list) %>% # List of singleton lists of row vectors.
    map(rep, nparam) %>% # List of lists with appropriate parameter numbers.
    map(Matrix::.bdiag) %>% # nm-list of block-diagonal matrices. 
    map(as.matrix) # FIXME: Maintain representation as sparse matrices?
    
  # Xl can now be used as covariates to vectorized MANOVA-style
  # parameters.

  ol <- offset %>% split(., row(.)) # Rows of offset as a list.

  nstats.all <- integer(nn)
  nstats.all[subset] <- nstats # So networks not in subset get 0 stats.
  inputs <- c(c(0,cumsum(nstats.all)), ms %>% map("inputs") %>% unlist())

  etamap <- function(x, n, ...){
    #' @importFrom purrr map map2 map2_dbl "%>%"
    Xl %>%
      map(`%*%`, c(x)) %>% # Evaluate X%*%c(x), where X is the predictor matrix and x is "theta".
      map(c) %>% # Output of the previous step is a matrix; convert to vector.
      map2(ol, `+`) %>% # Add on offset.
      map2(ms %>% map(c("model","etamap")), # Obtain the etamap.
           ergm.eta) %>% # Evaluate eta.
#      map2(weights, `*`) %>% # Weight.
      unlist()
  }
  etagradient <- function(x, n, ...){
    Xl %>%
      map(`%*%`, c(x)) %>% # Evaluate X%*%c(x), where X is the predictor matrix and x is "theta".
      map(c) %>% # Output of the previous step is a matrix; convert to vector.
      map2(ol, `+`) %>% # Add on offset.
      map2(ms %>% map(c("model","etamap")), # Obtain the etamap.
           ergm.etagrad) %>% # Evaluate eta gradient.
      map2(Xl, ., crossprod) %>%
#      map2(weights, `*`) %>% # Weight.
      do.call(cbind,.)
  }
  if(with(ms[[1]]$model$etamap,
          any(mintheta[!offsettheta]!=-Inf) || any(maxtheta[!offsettheta]!=+Inf))){
    warning("Submodel specified to N() operator with a linear model formula has parameter constraints. They will be ignored.")
  }
  params <- rep(list(NULL), nparam*ncol(xm))
  parnames <- colnames(xm)
  parnames <- ifelse(parnames=="(Intercept)", "1", parnames)
  names(params) <- paste0('N(',NVL3(a$label,paste0(.,","),""),rep(parnames, nparam),'):',rep(param_names(ms[[1]]$model, canonical=FALSE), each=ncol(xm)))
  coef.names <- paste0('N#',rep(seq_len(nm), nstats),':',unlist(lapply(lapply(ms, `[[`, "model"), param_names, canonical=TRUE)))
  gs <- unlist(lapply(ms, `[[`, "gs"))

  ## FIXME: The term needs to divide the each network's statistics by
  ## its weight above, to compensate for the fact that eta and
  ## gradient are multiplied.
  list(name="MultiNets", coef.names = coef.names, inputs=inputs, dependence=!all(sapply(lapply(ms, `[[`, "model"), is.dyad.independent)), emptynwstats = gs, auxiliaries = auxiliaries, map = etamap, gradient = etagradient, params = params)

  ## TODO: Re-add the optimised special case for when weights are fixed and constant, and all active models have the same numbers of coefficients.
  ## ### Test if weights are constant within the network.
  ## ### Fixed weights for networks
  ##   w <- switch(mode(a$weight),
  ##               character = sapply(nwl[subset], get.network.attribute, a$weight),
  ##               numeric = rep(w, length.out=sum(subset)),
  ##               `function` = sapply(nwl[subset], a$weight),
  ##               `NULL` = rep(1,nn))
  ##   if(!all_identical(nparams) || !all_identical(nstats)) stop("N() operator only supports models with the same numbers of parameters and coefficients for every network. This may change in the future.")
    
  ##   inputs <- c(w, unlist(lapply(ms, `[[`, "inputs")))
  ##   gs <- c(sapply(ms, `[[`, "gs")%*%w)
    
  ##   c(list(name="MultiNet", coef.names = paste0("N(",NVL(a$wname,"sum"),")*",ms[[1]]$model$coef.names), inputs=inputs, dependence=!is.dyad.independent(ms[[1]]$model), emptynwstats = gs, auxiliaries = auxiliaries),
  ##     passthrough.curved.ergm_model(ms[[1]]$m, function(x) paste0('N(',NVL(a$wname,"sum"),")*",x)))
  ## }
}
