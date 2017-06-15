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
  a <- .pop_vattrv(nw, ".LayerID")$vattr
  n <- length(a)
  bip <- nw %n% "bipartite"
  if(NVL(bip,0)){
    ea <- a[seq_len(bip)]
    aa <- a[bip+seq_len(n-bip)]
    el <- rle(ea)$lengths
    al <- rle(aa)$lengths
    if(!all.same(el) || !all.same(al)) stop("Layers must be networks of the same dimensions.", call.=FALSE)

    c(length(el), a, seq_len(n) - c((ea-1)*el[1], (aa-1)*al[1]))
    
  }else{
    l <- rle(a)$lengths
    if(!all.same(l)) stop("Layers must be networks of the same size.", call.=FALSE)
    c(length(l), a, seq_len(n) - (a-1)*l[1])
  }
}

InitErgmTerm..layer.nets <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  list(name="_layer_nets", coef.names=c(), inputs=.layer_vertexmap(nw), dependence=FALSE)
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

  # Process layer specification
  namemap <- seq_along(nwl)
  names(namemap) <- nwnames

  lIDs <- term.list.formula(a$layers[[2]])
  lIDs <- sapply(lIDs, function(l)
    switch(class(l),
           numeric = l,
           namemap[as.character(l)])
    )

  nl <- length(lIDs)           

  nwl <- nwl[lIDs]
  
  ml <- lapply(nwl, function(nw, ...){
    if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
    else nw <- ergm.getnetwork(f)
    
    ergm.getmodel(f, nw, response=response,...)
  }, ...)

  inputsl <- mapply(function(nw, m){
    Clist <- ergm.Cprepare(nw, m, response=response)
    inputs <- pack.Clist_as_num(Clist)
  }, nwl, ml, SIMPLIFY=FALSE)
  inputs <- c(nl, lIDs, unlist(inputsl))

  gsl <- lapply(ml, ergm.emptynwstats.model)
  gs <- Reduce(`+`, gsl)
  
  c(list(name="OnLayer", coef.names = paste0('OnLayer(',ml[[1]]$coef.names,',',deparse(a$layers),')'), inputs=inputs, dependence=!is.dyad.independent(ml[[1]]), emptynwstats = gs), auxiliaries = ~.layer.nets,
    passthrough.curved.ergm.model(ml[[1]], function(x) paste0('OnLayer(',x,',',deparse(a$layers),')')))
}
