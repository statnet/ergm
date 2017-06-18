Layer <- function(...){
  args <- list(...)
  if(length(args)==2 && is(args[[1]], "network") && is(args[[2]], "vector")){
    nwl <-
      lapply(args[[2]], function(eattr){
        nw <- args[[1]]
        el <- as.edgelist(nw, eattr)
        nw[,]<-0
        nw[el[,-3]] <- el[,3]
        nw
      })
    names(nwl) <- args[[2]]
  }else if(all(sapply(args, is, "network"))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multilayer specification. See help for information.")
  
  combine_networks(nwl, blockID.vattr=".LayerID", blockName.vattr=".LayerName")
}

#' Calculate a vector that maps the global LHS network Vertex indices within-layer Vertex and a Vertex to layer lookup table.
#' @noRd
.layer_vertexmap <- function(nw){
  a <- .peek_vattrv(nw, ".LayerID")
  n <- length(a)
  bip <- nw %n% "bipartite"
  if(NVL(bip,0)){
    ea <- a[seq_len(bip)]
    aa <- a[bip+seq_len(n-bip)]
    el <- rle(ea)$lengths
    al <- rle(aa)$lengths
    if(!all.same(el) || !all.same(al)) stop("Layers must be networks of the same dimensions.", call.=FALSE)

    list(nl = length(el), lids = a, lmap = seq_len(n) - c((ea-1)*el[1], bip - el[1] + (aa-1)*al[1]))
    
  }else{
    l <- rle(a)$lengths
    if(!all.same(l)) stop("Layers must be networks of the same size.", call.=FALSE)
    c(nl = length(l), lids = a, lmap = seq_len(n) - (a-1)*l[1])
  }
}

InitErgmTerm..layer.nets <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  list(name="_layer_nets", coef.names=c(), inputs=unlist(.layer_vertexmap(nw)), dependence=FALSE)
}

InitErgmTerm..layer.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("layers"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  
  nwl <- uncombine_network(nw, split.vattr=".LayerID", names.vattr=".LayerName")
  nwnames <- names(nwl)

  # Process layer specification
  namemap <- seq_along(nwl)
  names(namemap) <- nwnames

  if(length(term.list.formula(a$layers[[2]]))!=1) stop("Currently, the .layer.net() auxiliary formula must have exactly one term.", call.=FALSE)
  
  ll <- pack.LayerLogic_formula_as_double(a$layers, namemap)

  if(any(sapply(ll, test_eval.LayerLogic, FALSE))) stop("Layer specifications that produce edges on the output layer for empty input layers are not supported at this time.", call.=FALSE)
  
  list(name="_layer_net", coef.names=c(), inputs=c(unlist(.layer_vertexmap(nw)),unlist(ll)), dependence=FALSE)
}

pack.LayerLogic_formula_as_double <- function(formula, namemap){
  OPMAP <- c(`(` = 0,
             `!` = -1,
             `&` = -2,
             `&&` = -2,
             `|` = -3,
             `||` = -3,
             `==` = -4,
             `!=` = -5,
             `xor` = -5)
  
  lterms <- term.list.formula(formula[[2]])

  lidMap <- function(l){
    switch(class(l),
           numeric = l,
           namemap[as.character(l)])
  }
  
  postfix <- function(call, coml=c()){
    if(is.call(call)){
      op <- call[[1]]
      for(i in seq_along(call[-1])+1){
        coml <- c(coml, postfix(call[[i]]))
      }
      coml <- c(coml, OPMAP[as.character(op)])
    }else{
      coml <- c(coml, lidMap(call))
    }
    coml[coml!=0]
  }
  
  o <- lapply(lterms, postfix)
  lapply(o, function(com) c(length(com), com))
}

test_eval.LayerLogic <- function(commands, lv){
  coms <- commands[-1]
  lv <- rep(lv, length.out=max(coms))
  stack <- c()
  if(length(coms)!=commands[1]) stop("Layer specification command vector specifies incorrect number of commands.", call.=FALSE)
  for(i in 1:commands[1]){
    com <- coms[i]
    if(com==-1){
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
      stack <- c(x0 == y0, stack)
    }else if(com==-5){
      x0 <- stack[1]; stack <- stack[-1]
      y0 <- stack[1]; stack <- stack[-1]
      stack <- c(x0 != y0, stack)
    }else{
      stack <- c(lv[com], stack)
    }
  }
  if(length(stack)!=1) stop("Invalid layer specification command sequence.", call.=FALSE)
  stack
}

.all_layers_terms <- function(n){
  if(n<1) return(NULL)
  if(n==1) return(1)
  o <- call("+", 1, 2)
  for(i in seq_len(n-2)+2) o <- call("+", o, i)
  o
}

.mk_.layer.net_auxform <- function(trms, nl){
  # Replace . with all layers.
  trms <- do.call(substitute, list(trms, list(`.`=.all_layers_terms(nl))))
  # Get the formula as a list of term calls.
  trms <- term.list.formula(trms[[2]])
  trmcalls <- lapply(trms, function(ltrm) as.formula(call("~", ltrm)))
  trmcalls <- lapply(trmcalls, function(ltrm) call(".layer.net", ltrm))
  append.rhs.formula(~.,trmcalls)[-2]
}

InitErgmTerm.OnLayer <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "layers"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, ~.),
                      required = c(TRUE, FALSE))

  f <- a$formula

  nwl <- uncombine_network(nw, split.vattr=".LayerID", names.vattr=".LayerName")
  nwnames <- names(nwl)

  layers <- a$layers
  auxiliaries <- .mk_.layer.net_auxform(layers, length(nwl))
  nltrms <- length(term.list.formula(auxiliaries[[2]]))

  nw1 <- nwl[[1]]
  
  if(length(f)==2) f <- nonsimp.update.formula(f, nw1~.)
  else nw1 <- ergm.getnetwork(f)
  
  m <- ergm.getmodel(f, nw1, response=response,...)
  
  Clist <- ergm.Cprepare(nw1, m, response=response)
  inputs <- pack.Clist_as_num(Clist)
  
  inputs <- c(nltrms, inputs)

  gs <- ergm.emptynwstats.model(m) * nltrms
  
  c(list(name="OnLayer", coef.names = paste0('OnLayer(',m$coef.names,',',deparse(a$layers),')'), inputs=inputs, dependence=!is.dyad.independent(m), emptynwstats = gs, auxiliaries = auxiliaries),
    passthrough.curved.ergm.model(m, function(x) paste0('OnLayer(',x,',',deparse(a$layers),')')))
}

InitErgmTerm.layerCMB <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("layers"),
                      vartypes = c("formula"),
                      defaultvalues = list(~.),
                      required = c(FALSE))

  nwl <- uncombine_network(nw, split.vattr=".LayerID", names.vattr=".LayerName")
  layers <- a$layers
  auxiliaries <- .mk_.layer.net_auxform(layers, length(nwl))
  nltrms <- length(term.list.formula(auxiliaries[[2]]))

  inputs <- c(nltrms)

  list(name="layerCMB", coef.names = paste0('layerCMB(',deparse(a$layers),')'), inputs=inputs, dependence=FALSE, auxiliaries = auxiliaries, emptynwstats = network.dyadcount(nwl[[1]], FALSE)*lfactorial(nltrms))
}
