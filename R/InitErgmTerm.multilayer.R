## TODO: LL-Constrained proposals.
## TODO: Check that noncommutative LL operators work as intended.
.unenv <- function(f){
  environment(f) <- NULL
  f
}

.despace <- function(s) gsub("[[:space:]]", "", s)

.lspec_coef.names <- function(Llist, collapse=TRUE){
  reprs <- sapply(seq_along(Llist), function(l){
    name <- names(Llist)[l]
    L <- Llist[[l]]
    s <- NVL3(L, toString(ergm_LayerLogic(.)), "")
    if(NVL(name,"")!="") s <- paste0(name,"=",s)
    s
  })
  if(collapse) paste0("L(",paste0(reprs, collapse=","),")") else reprs
}

#' Construct a "view" of a network.
#'
#' Returns a network with edges optionally filtered according to a
#' specified criterion and with edge attributes optionally computed
#' from other edge attributes.
#'
#' @param x a [`network`] object.
#' @param ... a list of attribute or filtering specifications. See
#'   Details.
#' @param .clear whether the edge attributes not set by this call
#'   should be deleted.
#' @param .sep when specifying via a character vector, use this as the
#'   separator for concatenating edge values.
#'
#' @details Attribute specification arguments have the form
#'   `<newattrname> = <expr>`, where `<newattrname>` specifies the
#'   name of the new edge attribute (or attribute to be overwritten)
#'   and `<expr>` can be one of the following:
#' \describe{
#'
#' \item{a function}{The function will be passed two arguments, the
#' edgelist [`tibble`] and the network, and must return a vector of
#' edge attribute values to be set on the edges in the order
#' specified.}
#'
#' \item{a formula}{The expression on the RHS of the formula will be
#' evaluated with names in it referencing the edge attributes. The
#' input network may be referenced as `.nw`. The expression's result
#' is expected to be a vector of edge attribute values to be set on
#' the edges in the order specified.}
#' 
#' \item{a character vector}{If of length one, the edge attribute with
#' that name will simply be copied; if greater than one, the attribute
#' values will be concatenated wtih the `.sep` argument as the
#' separator.}
#'
#' \item{an object enclosed in [I()]}{The object will be used directly
#' to set the edge attribute.}
#' }
#'
#' Filtering arguments are specified the same way as attribute
#' arguments, but they must be named arguments (i.e., must be passed
#' without the `=`) and must return a logical or numeric vector
#' suitable for indexing the edge list. Multiple filtering arguments
#' will be specified, and the edge will be kept if it satisfise
#' *all*. If the conjunction of the edge's original states and the
#' filtering results is ambiguous (i.e., `NA`), it will be set as
#' missing.
#'
#' @return A [`network`] object with modified edges and edge attributes.
#'
#' @examples
#' data(florentine)
#' flo <- flomarriage
#' flo[,,add.edges=TRUE] <- as.matrix(flomarriage) | as.matrix(flobusiness)
#' flo[,, names.eval="m"] <- as.matrix(flomarriage)==1
#' flobusiness[3,5] <- NA
#' flo[,, names.eval="b"] <- as.matrix(flobusiness)==1
#' flo
#' (flob <- network_view(flo, "b"))
#' (flobusiness) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness),as.matrix(flob))
#' }
#' 
#' (flob <- network_view(flo, ~b&m))
#' (flobusiness & flomarriage) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness & flomarriage),as.matrix(flob))
#' }
#'
#' as.matrix(flob <- network_view(flo, bm=~b+m), attrname="bm")
#' (as.matrix(flobusiness) + as.matrix(flomarriage)) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness)+as.matrix(flomarriage),as.matrix(flob, attrname="bm"))
#' }
#'
#' as.matrix(flob <- network_view(flo, ~b, bm=~b+m), attrname="bm")
#' as.matrix(flobusiness)*(1+as.matrix(flomarriage)) # for comparison
#' \dontshow{
#' if(require(testthat, quietly=TRUE))
#' testthat::expect_equivalent(as.matrix(flobusiness)*(1+as.matrix(flomarriage)),as.matrix(flob, attrname="bm"))
#' }
#' 
#' 
#' @export
network_view <- function(x, ..., .clear=FALSE, .sep="."){
  # Handle empty network
  if(network.edgecount(x,na.omit=FALSE)==0) return(x)
  
  exprs <- list(...)
  fes <- exprs[NVL(names(exprs),rep("",length(exprs)))==""]
  oes <- exprs[NVL(names(exprs),rep("",length(exprs)))!=""]

  #' @importFrom rlang abort
  evl <- function(e, el, x){
    switch(class(e),
           `function` = e(.el=el, .nw=x),
           formula = eval(e[[length(e)]], envir=c(list(.nw=x), as.list(el)), enclos=environment(e)),
           character = if(length(e)==1) el[[e]] else do.call(paste, c(as.list(el[e]), sep=.sep)),
           AsIs = e,
           abort("Unsupported specification for network_view."))
  }

  if(length(fes)){
    el <- as_tibble(x, attrnames=list.edge.attributes(x), na.rm=FALSE)
    keep <- !el$na
    for(e in fes){
      keep <- keep & evl(e, el, x)
    }
    del <- na.omit(el$.eid[!keep])
    nael <- el$.eid[is.na(keep)]
    if(length(del)) delete.edges(x, del)
    if(length(nael)) set.edge.attribute(x, "na", TRUE, nael)
    # This one applies to a rare situation where the edge is missing
    # but the filter says that it should be present.
    add <- el$.eid[el$na & !is.na(keep) & keep]
    if(length(add)) set.edge.attribute(x, "na", FALSE, add)
  }

  for(i in seq_along(oes)){
    el <- as_tibble(x, attrname=list.edge.attributes(x), na.rm=FALSE)
    
    e <- oes[[i]]
    nm <- names(oes)[[i]]
    
    newval <- evl(e, el, x)
    set.edge.attribute(x, nm, newval, el$.eid)
  }
  
  if(.clear) for(a in setdiff(list.edge.attributes(x), c(nm,"na"))) delete.edge.attribute(x, a)
  x
}


#' Returns a directed version of an undirected binary network
#'
#' @param x a [`network`] object.
#' @param rule a string specifying how the network is to be
#'   constructed.
direct.network <- function(x, rule=c("both", "upper", "lower")){
  rule <- match.arg(rule)

  el <- as.edgelist(x)
  el <- switch(rule,
               both = rbind(el, el[,2:1,drop=FALSE]),
               upper = cbind(pmin(el[,1],el[,2]),pmax(el[,1],el[,2])),
               lower = cbind(pmax(el[,1],el[,2]),pmin(el[,1],el[,2])))
  
  o <- network.initialize(network.size(x), directed=TRUE, bipartite=x%n%"bipartite", loops=has.loops(x), hyper=is.hyper(x), multiple=is.multiplex(x))
  o <- network.edgelist(el, o)
  nvattr.copy.network(o, x)
}

#' A multilayer network representation.
#'
#' A function for specifying the LHS of a multilayer (a.k.a. multiplex
#' a.k.a. multirelational a.k.a. multivariate) ERGM.
#'
#' @param ... layer specification, in one of three formats:
#' 
#'   1. An (optionally named) list of identically-dimensioned
#'      networks.
#'
#'   1. Several networks as (optionally named) arguments.
#'
#'   1. A single network, a chararacter vector, and an optional
#'      logical vector. Then, the layers are values of the named edge
#'      attributes. If the network is directed, the logical vector
#'      specifies which of the layers should be treated as undirected.
#'
#' @return A network object with layer metadata.
#'
#' @seealso [Help on model specification][ergm-terms] for specific terms.
#' 
#' @examples
#'
#' data(florentine)
#'
#' # Method 1: list of networks
#' flo <- Layer(list(m = flomarriage, b = flobusiness))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 2: networks as arguments
#' flo <- Layer(m = flomarriage, b = flobusiness)
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' # Method 3: edge attributes:
#' flo <- flomarriage | flobusiness
#' flo[,, names.eval="m"] <- as.matrix(flomarriage)
#' flo[,, names.eval="b"] <- as.matrix(flobusiness)
#' flo <- Layer(flo, c("m","b"))
#' ergm(flo ~ L(~edges, ~m)+L(~edges, ~b))
#'
#' @export
Layer <- function(...){
  args <- list(...)
  if(length(args)%in%2:3 && is(args[[1]], "network") && is(args[[2]], "vector") && (length(args)!=3 || is(args[[3]], "vector"))){
    nwl <-
      lapply(args[[2]], function(eattr){
        network_view(args[[1]], eattr)
      })
    if(length(args)==3){
      symm <- as.logical(args[[3]])
      for(i in which(symm)){
        # There is probably a more efficient way to do this, but we need
        # to compute nw1 anyway.
        nw1 <- symmetrize(nwl[[i]], rule="weak")
        nw2 <- symmetrize(nwl[[i]], rule="strong")
        if(!identical(as.vector(as.edgelist(nw1)),as.vector(as.edgelist(nw2))))
          stop("Layer specified to be treated as undirected is not symmetric.")
        nwl[[i]] <- nw1
      }
    }
    names(nwl) <- args[[2]]
  }else if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multilayer specification. See help for information.")

  # nwl may now be a list with networks of heterogeneous directedness.
  
  dir <- sapply(nwl, is.directed)
  symm <- if(all_identical(dir)) rep(FALSE, length(nwl)) else !dir

  nwl <- mapply(function(nw, symm) {
    nw %v% ".undirected" <- symm
    if(symm) direct.network(nw,rule="upper") else nw
  }, nwl, symm, SIMPLIFY=FALSE)

  # nwl is now a list of networks with homogeneous directedness, some
  # networks tagged with vertex attribute .undirected.
  
  if(any(!is.na(suppressWarnings(sapply(names(nwl), as.numeric))))) warning("Using numeric layer names is ambiguous.")

  constraintsl <- lapply(nwl, get.network.attribute, "constraints")
  if(!all_identical(lapply(constraintsl, .unenv))) stop("Layers have differing constraint structures. This is not supported at this time.")
  obs.constraintsl <- lapply(nwl, get.network.attribute, "obs.constraints")
  if(!all_identical(lapply(obs.constraintsl, .unenv))) stop("Layers have differing observation processes. This is not supported at this time.")
  
  nw <- combine_networks(nwl, blockID.vattr=".LayerID", blockName.vattr=".LayerName", ignore.nattr = c(eval(formals(combine_networks)$ignore.nattr), "constraints", "obs.constraints"), subnet.cache=TRUE)
  nw %n% "constraints" <-
      if(NVL(nwl[[1]]%n%"constraints",~.)==~.)
        ~blockdiag(".LayerID")
      else
        append_rhs.formula(nwl[[1]]%n%"constraints", list(call("blockdiag",".LayerID")), TRUE)
  
  if(any(symm)) nw %n% "constraints" <- append_rhs.formula(nw%n%"constraints", list(call("upper_tri",".undirected")), TRUE)

  if("obs.constraints" %in% list.network.attributes(nwl[[1]])) nw %n% "obs.constraints" <- nwl[[1]]%n%"obs.constraints"

  nw
}


## InitErgmTerm..layer.nets <- function(nw, arglist, response=NULL, ...){
##   a <- check.ErgmTerm(nw, arglist,
##                       varnames = c(),
##                       vartypes = c(),
##                       defaultvalues = list(),
##                       required = c())
##   list(name="_layer_nets", coef.names=c(), inputs=unlist(.block_vertexmap(nw, ".LayerID", TRUE)), dependence=FALSE)
## }

.depends_on_layers <- function(commands){
  coms <- commands[-1]
  if(any(coms==0)) coms <- coms[-(which(coms==0)+1)] # Drop all numeric literals (i.e., numbers preceded by 0).
  coms <- coms[coms>=1] # Drop all commands.
  unique(coms)
}

InitErgmTerm..layer.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("L"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  
  nwl <- .split_constr_network(nw,".LayerID",".LayerName")
  nwnames <- names(nwl)

  # Process layer specification
  namemap <- seq_along(nwl)
  names(namemap) <- nwnames
  
  ll <- to_ergm_Cdouble(ergm_LayerLogic(a$L, namemap))
  # Terms on this logical layer will induce dyadic independence if its
  # value depends on more than one other layer value.
  dependence <- length(.depends_on_layers(ll))>1
  
  if(test_eval.LayerLogic(ll, FALSE)) stop("Layer specifications that produce edges on the output layer for empty input layers are not supported at this time.", call.=FALSE)
  
  list(name="_layer_net", coef.names=c(), inputs=c(unlist(.block_vertexmap(nw, ".LayerID", TRUE)),if(is.directed(nw)) sapply(nwl, function(nw) (nw%v% ".undirected")[1]), ll), dependence=dependence)
}

LL_PREOPMAP <- list(
    # Unary operators
    c(`t` = -23)
  )
LL_POSTOPMAP <- list(
    # Unary operators
    c(`(` = NA,
      `!` = -1,
      `+` = NA,
      `-` = -16,
      `abs` = -17,
      `round` = -20,
      `sign` = -22),
    # Binary operators
    c(`&` = -2,
      `&&` = -2,
      `|` = -3,
      `||` = -3,
      `xor` = -4,
      `==` = -5,
      `!=` = -6,
      `<` = -7,
      `>` = -8,
      `<=` = -9,
      `>=` = -10,
      `+` = -11,
      `-` = -12,
      `*` = -13,
      `/` = -14,
      `%%` = -15,
      `^` = -18,
      `%/%` = -19,
      `round` = -21)
    )


#' Internal representation of Layer Logic
#'
#' @param formula A Layer Logic formula.
#' @param namemap A character vector giving the names, of the layers
#'   referenced, or `NULL`.
#'
#' @return A structure with nonce class
#'   `c("ergm_LayerLogic",class(formula))`, comprising the input
#'   `formula` and an attribute `namemap` containing the `namemap`.
#' @keywords internal
#' @export
ergm_LayerLogic <- function(formula, namemap=NULL){
  ## TODO: Check whether we should verify that this is a formula.
  structure(formula, namemap=namemap, class=c("ergm_LayerLogic", class(formula)))
}

#' @describeIn ergm_LayerLogic A method to generate coefficient names
#'   associated with the Layer Logic.
#' @param x An `ergm_LayerLogic` object.
#' @param ... Additional arguments, currently unused.
#' @export
toString.ergm_LayerLogic <- function(x, ...){
  class(x) <- keep(class(x), `!=`, "ergm_LayerLogic")
  fmt <- function(x)
    switch(class(x),
           formula = .despace(deparse(if(length(x)==2) x[[2]] else x)),
           character = x,
           list = paste0('(',paste(sapply(x,fmt),collapse=","),')'),
           as.character(x))
  fmt(x)
}

#' @describeIn ergm_LayerLogic A method to encode and serialize the
#'   Layer Logic into a postfix program understood by the C code.
#' @export
to_ergm_Cdouble.ergm_LayerLogic <- function(x, ...){
  formula <- x
  namemap <- attr(x, "namemap")
    
  lidMap <- function(l){
    switch(class(l),
           numeric = c(0,l),
           character =,
           name = if(regexpr('^[0-9]+$',l)!=-1) as.integer(as.character(l))
                  else namemap[as.character(l)])
  }

  preops <- 0
  postfix <- function(call, coml=c()){
    if(is.call(call)){
      op <- call[[1]]
      if(as.character(op) %in% unlist(lapply(LL_PREOPMAP, names))){
        preops <<- preops+1
        coml <- c(coml, LL_PREOPMAP[[length(call)-1]][[as.character(op)]])
        postop <- FALSE
      }else postop <- TRUE
      for(i in seq_along(call[-1])+1){
        coml <- c(coml, postfix(call[[i]]))
      }
      if(postop) coml <- c(coml, LL_POSTOPMAP[[length(call)-1]][as.character(op)])
    }else{
      coml <- c(coml, lidMap(call))
    }
    coml[!is.na(coml)]
  }
  
  com <- postfix(formula[[length(formula)]])
  c(sum(com!=0 & !com%in%unlist(LL_PREOPMAP)), com)
}

test_eval.LayerLogic <- function(commands, lv, lvr = lv){
  coms <- commands[-1]
  lv <- rep(lv, length.out=max(coms))
  stack <- c()
  if(sum(coms!=0 & !coms%in%unlist(LL_PREOPMAP))!=commands[1]) stop("Layer specification command vector specifies incorrect number of commands.", call.=FALSE)
  for(i in 1:commands[1]){
    com <- coms[1]
    if(com==0){
      coms <- coms[-1]
      com <- coms[1]
      stack <- c(com, stack)
    }else if(com==-1){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(!x0, stack)
    }else if(com==-2){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 && y0, stack)
    }else if(com==-3){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 || y0, stack)
    }else if(com==-4){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(xor(x0, y0), stack)
    }else if(com==-5){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 == y0, stack)
    }else if(com==-6){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 != y0, stack)
    }else if(com==-7){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 < y0, stack)
    }else if(com==-8){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 > y0, stack)
    }else if(com==-9){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 <= y0, stack)
    }else if(com==-10){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 >= y0, stack)
    }else if(com==-11){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 + y0, stack)
    }else if(com==-12){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 - y0, stack)
    }else if(com==-13){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 * y0, stack)
    }else if(com==-14){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 / y0, stack)
    }else if(com==-15){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %% y0, stack)
    }else if(com==-16){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(-x0, stack)
    }else if(com==-17){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(abs(x0), stack)
    }else if(com==-18){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 ^ y0, stack)
    }else if(com==-19){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 %/% y0, stack)
    }else if(com==-20){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(round(x0), stack)
    }else if(com==-21){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(round(x0, y0), stack)
    }else if(com==-22){
      x0 <- stack[1]; stack <- stack[-1]
      stack <- c(sign(x0), stack)
    }else if(com==-23){ 
      coms <- coms[-1]
      x0 <- coms[1]
      stack <- c(lvr[x0], stack)
    }else{
      stack <- c(lv[com], stack)
    }
    coms <- coms[-1]
  }
  if(length(stack)!=1) stop("Invalid layer specification command sequence.", call.=FALSE)
  stack
}

.all_layers_terms <- function(n, LHS=NULL){
  if(is.null(LHS))
    lapply(seq_len(n), function(i) as.formula(substitute(~i,list(i=as.name(i)))))
  else
    lapply(seq_len(n), function(i) as.formula(substitute(lhs~i,list(lhs=LHS, i=as.name(i)))))
}

.mk_.layer.net_auxform <- function(ll, nl){
  trmcalls <- .layers_expand_dot(ll, nl)
  # Get the formula as a list of term calls.
  trmcalls <- lapply(trmcalls, function(ltrm) call(".layer.net", ltrm))
  append_rhs.formula(~.,trmcalls)[-2]
}

.layers_expand_dot <- function(ll, nl){
  if(is(ll, "formula")) ll <- list(ll)
  # Replace . with all layers.
  do.call(c, lapply(ll, function(f) if(f[[length(f)]]=='.') .all_layers_terms(nl, LHS = if(length(f)==3) f[[2]]) else list(as.formula(f))))
}

InitErgmTerm.L <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "Ls"),
                      vartypes = c("formula", "formula,list"),
                      defaultvalues = list(NULL, ~.),
                      required = c(TRUE, FALSE))
  f <- a$formula

  nwl <- .split_constr_network(nw,".LayerID",".LayerID")
  nwnames <- names(nwl)

  Ls <- a$Ls
  if(is(Ls, "formula")) Ls <- list(Ls)
  Ls.dotexp <- .layers_expand_dot(Ls, length(nwl))
  

  auxiliaries <- .mk_.layer.net_auxform(Ls, length(nwl))  
  nltrms <- length(list_rhs.formula(auxiliaries))

  w <- rep(1,nltrms)
  have.LHS <- sapply(Ls.dotexp, length)==3
  w[have.LHS] <- as.numeric(sapply(lapply(Ls.dotexp[have.LHS], "[[", 2), eval,environment(Ls[[1]])))
  
  nw1 <- nwl[[1]]
  
  if(length(f)==2) f <- nonsimp_update.formula(f, nw1~.)
  else nw1 <- ergm.getnetwork(f)
  
  m <- ergm_model(f, nw1, response=response,...)

  dependence <- !is.dyad.independent(m) || !is.dyad.independent(nonsimp_update.formula(auxiliaries, nw~., from.new="nw"))

  
  inputs <- to_ergm_Cdouble(m)
  
  inputs <- c(nltrms, w, inputs)

  gs <- summary(m) * nltrms
  
  c(list(name="OnLayer", coef.names = paste0(.lspec_coef.names(list(a$Ls)),":",m$coef.names), inputs=inputs, dependence=dependence, emptynwstats = gs, auxiliaries = auxiliaries),
    passthrough.curved.ergm_model(m, function(x) paste0(.lspec_coef.names(list(a$Ls)),":",x)))
}

InitErgmTerm.lCMB <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("Ls"),
                      vartypes = c("formula,list"),
                      defaultvalues = list(~.),
                      required = c(FALSE))

  nwl <- .split_constr_network(nw,".LayerID",".LayerID")
  Ls <- a$Ls
  auxiliaries <- .mk_.layer.net_auxform(Ls, length(nwl))
  nltrms <- length(list_rhs.formula(auxiliaries))

  inputs <- c(nltrms)

  list(name="layerCMB", coef.names = paste0('lCMB(',.despace(deparse(Ls)),')'), inputs=inputs, dependence=TRUE, auxiliaries = auxiliaries)
}

################################################################################
InitErgmTerm.ldegree<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("d", "by", "Ls", "dir"),
                      vartypes = c("numeric", "character", "formula,list", "character"),
                      defaultvalues = list(NULL, NULL, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  d<-a$d; byarg <- a$by; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree() has only one value", call.=FALSE)
  }
  if(!is.null(byarg) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  if(is.null(byarg)) {
    if(length(d)==0){return(NULL)}
    coef.names <- paste("degree",d,sep="")
    name <- "degree"
    inputs <- c(d)
  } else if (homophily) {
    if(length(d)==0){return(NULL)}
    # See comment in d_degree_w_homophily function
    coef.names <- paste("deg", d, ".homophily.",byarg, sep="")
    name <- "degree_w_homophily"
    inputs <- c(d, nodecov)
  } else {
    if(ncol(du)==0) {return(NULL)}
    #  No covariates here, so "ParamsBeforeCov" unnecessary
    # See comment in d_degree_by_attr function
    coef.names <- paste("deg", du[1,], ".", byarg,u[du[2,]], sep="")
    name <- "degree_by_attr"
    inputs <- c(as.vector(du), nodecov)
  }

  if(is.null(a$Ls) || is.null(a$dir)) stop("A layer and direction specification is required for the ldegree term.")
  
  c(.process_layers_degree(nw, a, name=name,coef.names=coef.names, inputs=inputs, emptynwstats=emptynwstats), minval = 0)
}

################################################################################
InitErgmTerm.gwldegree<-function(nw, arglist,  ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("decay", "fixed", "attrname","cutoff", "levels", "Ls", "dir"),
                      vartypes = c("numeric", "logical", "character","numeric", "character,numeric,logical", "formula,list", "character"),
                      defaultvalues = list(NULL, FALSE, NULL, 30, NULL, NULL),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  cutoff<-a$cutoff

  if(is.null(a$Ls) || is.null(a$dir)) stop("A layer and direction specification is required for the gwldegree term.")

  # d <- 1:(network.size(nw)-1)
   maxesp <- min(cutoff,network.size(nw)-1)
  d <- 1:maxesp
  if (!is.null(attrname) && !fixed ) {
    stop("The gwldegree term cannot yet handle a nonfixed decay ",
            "term with an attribute. Use fixed=TRUE.", call.=FALSE)
    
  }
  if(!fixed){ # This is a curved exponential family model
    if(!is.null(a$decay)) warning("In term 'gwldegree': decay parameter 'decay' passed with 'fixed=FALSE'. 'decay' will be ignored. To specify an initial value for 'decay', use the 'init' control parameter.", call.=FALSE)
    ld<-length(d)
    if(ld==0){return(NULL)}
    c(.process_layers_degree(nw, a, name="odegree", coef.names=paste("gwldegree#",d,sep=""), inputs=c(d)),
      list(params=list(gwldegree=NULL,gwldegree.decay=decay), conflicts.constraints="degreedist"), GWDECAY)
  } else {
    if(is.null(a$decay)) stop("Term 'gwldegree' with 'fixed=TRUE' requires a decay parameter 'decay'.", call.=FALSE)

    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwldegree")
      u<-NVL(a$levels, sort(unique(nodecov)))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u,nomatch=length(u)+1) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwldegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(NULL)}
      #  No covariates here, so "ParamsBeforeCov" unnecessary
      name <- "gwldegree_by_attr"
      coef.names <- paste("gwodeg", decay, ".", attrname, u, sep="")
      inputs <- c(decay, nodecov)
    }else{
      name <- "gwldegree"
      coef.names <- paste("gwodeg.fixed.",decay,sep="")
      inputs <- c(decay)
    }
    c(.process_layers_degree(nw, a, name=name, coef.names=coef.names, inputs=inputs), conflicts.constraints="degreedist")
  }
}


################################################################################
InitErgmTerm.twostarL<-function(nw, arglist,  ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = c("Ls", "type", "distinct"),
                      vartypes = c("formula,list", "character", "logical"),
                      defaultvalues = list(NULL, NULL, TRUE),
                      required = c(TRUE, TRUE, FALSE))
  TYPES <- c("out", "in", "path")
  type <- match.arg(tolower(a$type), TYPES)
  typeID <- match(type, TYPES)

  Ls <- a$Ls
  if(is(Ls, "formula")) Ls <- list(Ls)
  Ls <- rep(Ls, length.out=2)

  nlayers <- length(unique(.peek_vattrv(nw, ".LayerID")))
  auxiliaries <- .mk_.layer.net_auxform(Ls, nlayers)
  reprs <- .lspec_coef.names(Ls, collapse=FALSE)
  coef.names <- paste0("twostarL(",
                        switch(type,
                               out = paste0(reprs, collapse="<>"),
                               `in` = paste0(reprs, collapse="><"),
                               path = paste0(reprs, collapse=">>")),
                       if(a$distinct) ",distinct",
                        ")")
  
  inputs <- c(typeID, a$distinct)
  list(name="twostarL", coef.names=coef.names, inputs=inputs, auxiliaries=auxiliaries, minval=0, dependence=TRUE)
}
