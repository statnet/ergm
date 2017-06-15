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
  }else if(all(sapply(args, is, network))){
    nwl <- args
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, network))){
    nwl <- args[[1]]
  }else stop("Unrecognized format for multilayer specification. See help for information.")
  
  combine.networks(nwl, blockID.vattr=".LayerID", blockName.vattr=".LayerName")
}

.pop_vattrv <- function(nw, vattr){
  av <- get.vertex.attribute(nw, vattr, unlist=FALSE)
  a <- sapply(av, "[", 1)
  rest <- lapply(av, "[", -1)
  nw <- set.vertex.attribute(nw, vattr, rest)
  
  list(nw = nw, vattr = a)
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

