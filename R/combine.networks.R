#  File R/combine.networks.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

#' A single block-diagonal network created by combining multiple networks
#'
#' Given a list of compatible networks, the [combine_networks()] returns a single
#' block-diagonal network, preserving attributes that can be
#' preserved.
#'
#' @param nwl a list of [`network::network`]s to be combined; they
#'   must have similar fundamental properties: directedness and
#'   bipartedness, though their sizes (and the size of each bipartite
#'   group) can vary.
#'
#' @param ignore.nattr,ignore.vattr,ignore.eattr network, vertex, and
#'   edge attributes not to be processed as described below.
#'
#' @param blockID.vattr name of the vertex attribute into which to store
#'   the index of the network to which that vertex originally belonged.
#'
#' @param blockName.vattr if not `NULL`, the name of the vertex
#'   attribute into which to store the name of the network to which
#'   that vertex originally belonged.
#' 
#' @param detect.edgecov if `TRUE`, combine network attributes that
#'   look like dyadic covariate ([`ergm::edgecov`]) matrices into a
#'   block-diagonal matrix.
#'
#' @param standardized whether the passed networks can be assumed to
#'   have been run throgh [ergm::standardize.network()]: their
#'   internal representation fulfills certain constraints.
#'
#' @param keep.unshared.attr whether to keep those network, vertex,
#'   and edge attributes not shared by all networks in the list; if
#'   \code{TRUE}, positions corresponding to networks lacking the
#'   attribute are replaced with \code{NA}, \code{NULL}, or some other
#'   placeholder; incompatible with \code{detect.edgecov==TRUE}.
#'
#' @return a [`network::network`] with a block-diagonal structure (or
#'   its bipartite equivalent) comprising the networks passed in
#'   `nwl`. In particular,
#'
#' * the returned network's size is the sum of the input networks';
#' 
#' * its basic properties (directedness and bipartednes) are the same;
#'
#' * the input networks' sociomatrices (both edge presence and edge
#'   attributes) are the blocks in the sociomatrix of the returned
#'   network;
#'
#' * vertex attributes are concatenated;
#'
#' * edge attributes are assigned to their respective edges in
#'   the returned network;
#'
#' * network attributes are stored in a list; but if
#'   `detect.edgecov==TRUE`, those network attributes that have the
#'   same dimension as the sociomatrices of the constituent networks,
#'   they are combined into a single block-diagonal matrix that is
#'   then stored as that attribute.
#'
#' In addition, a two new vertex attibutes, specified by
#' `blockID.vattr` and (optionally) `blockName.vattr` contain,
#' respectively, the index in `nwl` of the network from which that
#' vertex came and its name, determined as follows:
#' 
#' 1. If `nwl` is a named list, names from the list are used.
#'
#' 2. If not 1, but the network has an attribute `title`, it is used.
#' 
#' 3. Otherwise, a numerical index is used.
#' 
#' If `blockID.vattr` already exists on the constituent networks, the
#' index is *prepended* to the attribute.
#' 
#' @examples
#'
#' data(samplk)
#' 
#' o1 <- combine_networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#' head(get.vertex.attribute(o1, ".NetworkID"))
#' o2 <- combine_networks(list(o1, o1))
#' image(as.matrix(o2))
#' head(get.vertex.attribute(o2, ".NetworkID", unlist=FALSE))
#'
#' data(florentine)
#' f1 <- combine_networks(list(business=flobusiness, marriage=flomarriage), blockName.vattr=".NetworkName") 
#' image(as.matrix(f1))
#' head(get.vertex.attribute(f1, ".NetworkID"))
#' head(get.vertex.attribute(f1, ".NetworkName"))
#' @export
combine_networks <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
  if(any(sapply(nwl, is.bipartite))) .combine_networks.bipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, standardized=standardized, keep.unshared.attr=keep.unshared.attr)
  else .combine_networks.unipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, standardized=standardized, keep.unshared.attr=keep.unshared.attr)
}


.combine_networks.unipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
  if(any(diff(sapply(nwl, is.directed)))) stop("All networks must have the same directedness.")
  if(keep.unshared.attr && detect.edgecov) stop("Detection of edge covariates is not compatible with retaining unshared attributes.")
  attrset <- if(keep.unshared.attr) union else intersect
  
  if(!standardized) nwl <- lapply(nwl, standardize.network)
  ns <- sapply(nwl, network.size)
  blks <- c(0, cumsum(ns))

  out <- network.initialize(sum(ns), directed=is.directed(nwl[[1]]))


  # Concatenate network attributes. If you run into what looks like a covariate matrix, combine it correctly.
  
  for(a in setdiff(Reduce(attrset,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==ns)
       && all(sapply(vl, ncol)==ns)
       && all_identical(sapply(vl, mode))){

      # A logical vector that extracts off-diagonal element of the ns*ns matrix.


      offdiags <- unlist(lapply(ns, function(n) c(diag(1,n)==0)))
      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the off-diagonal elements
      # of each covariate matrix, concatenate them into one vector,
      # remove the NAs, and take 0 (if it's present) or the minimum
      # value. (0 as filler can be helpful for sparse matrix
      # representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c"))[offdiags])
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(ns), sum(ns))
      mode(m) <- mode(vl[[1]])
      
      for(b in seq_along(vl)){
        inds <- blks[b]+seq_len(ns[b])
        m[inds, inds] <- vl[[b]]
      }

      vl <- m
    }
    
    out <- set.network.attribute(out, a, vl)
  }

  # Concatenate vertex attributes.
  
  for(a in setdiff(Reduce(attrset,lapply(nwl, list.vertex.attributes)),
                          ignore.vattr)){ # I.e., iterate through common attributes.
    out <- set.vertex.attribute(out, a,
                                do.call(c, lapply(nwl, get.vertex.attribute, a, unlist=FALSE))
                                )
  }

  # Add ties and attributes

  for(b in seq_along(nwl)){
    el <- rbind(as.edgelist(nwl[[b]]),as.edgelist(is.na(nwl[[b]])))
    eids <- apply(el, 1, function(e) get.edgeIDs(nwl[[b]], e[1], e[2], na.omit=FALSE))

    vals <- lapply(nwl[[b]]$mel,"[[","atl")[eids]
    names <- lapply(vals, names)

    out <- add.edges(out, el[,1]+blks[b], el[,2]+blks[b], names.eval=names, vals.eval=vals)
  }

  # Finally, add a vertex attribute specifying the blocks

  b <- rep(seq_along(ns),ns)
  if(blockID.vattr %in% list.vertex.attributes(out)){ # blockID.vattr already exists
    b <- mapply(c, # Concatenate
                b, # each element of b
                get.vertex.attribute(out, blockID.vattr, unlist=FALSE), # with the corresponding element of out %v% blockID.vattr.
                SIMPLIFY=FALSE)
  }
  out <- set.vertex.attribute(out, blockID.vattr, b)

  if(!is.null(blockName.vattr)){
    bn <-
      if(!is.null(names(nwl))) names(nwl)
      else if("title" %in% list.network.attributes(out)) out %v% "title"
      else seq_along(ns)  
    b <- rep(bn,ns)
    if(blockName.vattr %in% list.vertex.attributes(out)){ # blockID.vattr already exists
      b <- mapply(c, # Concatenate
                  b, # each element of b
                  get.vertex.attribute(out, blockName.vattr, unlist=FALSE), # with the corresponding element of out %v% blockID.vattr.
                  SIMPLIFY=FALSE)
    }
    out <- set.vertex.attribute(out, blockName.vattr, b)
  }
  
  out
}


.combine_networks.bipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
  if(!all(sapply(nwl, is.bipartite))) stop("This function operates only on bipartite networks.")
  if(any(sapply(nwl, is.directed))) stop("Bipartite directed networks are not supported at this time.")
  if(keep.unshared.attr && detect.edgecov) stop("Detection of edge covariates is not compatible with retaining unshared attributes.")
  attrset <- if(keep.unshared.attr) union else intersect

  if(!standardized) nwl <- lapply(nwl, standardize.network)
  ns <- sapply(nwl, network.size)
  es <- sapply(nwl, "%n%", "bipartite")
  eblks <- c(0, cumsum(es))
  bip <- eblks[length(eblks)]
  ablks <- cumsum(c(bip, ns-es))

  out <- network.initialize(sum(ns), directed=is.directed(nwl[[1]]), bipartite=bip)

  # Concatenate network attributes. If you run into what looks like a covariate matrix, combine it correctly.
  
  for(a in setdiff(Reduce(attrset,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==es)
       && all(sapply(vl, ncol)==ns-es)
       && all_identical(sapply(vl, mode))){

      # It doesn't matter what the "filler" elements are, as long as
      # adding them doesn't add another category and it's not NA. So,
      # what this does is as follows: grab the elements of each
      # covariate matrix, concatenate them into one vector, remove the
      # NAs, and take 0 (if it's present) or the minimum value. (0 as
      # filler can be helpful for sparse matrix representations.)
      dummyvals <- na.omit(unlist(lapply(vl, "c")))
      dummyval <- if(0 %in% dummyvals) 0 else min(dummyvals)
      m <- matrix(dummyval, sum(es), sum(ns-es))
      mode(m) <- mode(vl[[1]])
      
      for(b in seq_along(vl)){
        einds <- eblks[b]+seq_len(es[b])
        ainds <- ablks[b]+seq_len(ns[b]-es[b])
        m[einds, ainds-sum(es)] <- vl[[b]]
      }

      vl <- m
    }
    
    out <- set.network.attribute(out, a, vl)
  }

  # Concatenate vertex attributes.
  
  for(a in setdiff(Reduce(attrset,lapply(nwl, list.vertex.attributes)),
                   ignore.vattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.vertex.attribute, a, unlist=FALSE)
    v <- vector(mode(vl[[1]]), sum(ns))

    for(b in seq_along(vl)){
        einds <- eblks[b]+seq_len(es[b])
        ainds <- ablks[b]+seq_len(ns[b]-es[b])
        v[einds] <- vl[[b]][seq_len(es[b])]
        v[ainds] <- vl[[b]][es[b]+seq_len(ns[b]-es[b])]
    }
    
    out <- set.vertex.attribute(out, a, v)
  }

  # Add ties and attributes

  for(b in seq_along(nwl)){
    el <- rbind(as.edgelist(nwl[[b]]),as.edgelist(is.na(nwl[[b]])))
    eids <- apply(el, 1, function(e) get.edgeIDs(nwl[[b]], e[1], e[2], na.omit=FALSE))

    vals <- lapply(nwl[[b]]$mel,"[[","atl")[eids]
    names <- lapply(vals, names)

    out <- add.edges(out, el[,1]+eblks[b], el[,2]-es[b]+ablks[b], names.eval=names, vals.eval=vals)
  }

  # Finally, add a vertex attribute specifying the blocks

  b <- rep(rep(seq_along(ns),2),c(es,ns-es))
  if(blockID.vattr %in% list.vertex.attributes(out)){ # blockID.vattr already exists
    b <- mapply(c, # Concatenate
                b, # each element of b
                get.vertex.attribute(out, blockID.vattr, unlist=FALSE), # with the corresponding element of out %v% blockID.vattr.
                SIMPLIFY=FALSE)
  }
  out <- set.vertex.attribute(out, blockID.vattr, b)

  if(!is.null(blockName.vattr)){
    bn <-
      if(!is.null(names(nwl))) names(nwl)
      else if("title" %in% list.network.attributes(out)) out %v% "title"
      else seq_along(ns)
    b <- rep(rep(bn,2),c(es,ns-es))
    if(blockName.vattr %in% list.vertex.attributes(out)){ # blockID.vattr already exists
      b <- mapply(c, # Concatenate
                  b, # each element of b
                  get.vertex.attribute(out, blockName.vattr, unlist=FALSE), # with the corresponding element of out %v% blockID.vattr.
                  SIMPLIFY=FALSE)
    }
    out <- set.vertex.attribute(out, blockName.vattr, b)
  }
  
  out
}

.peek_vattrv <- function(nw, vattr){
  av <- get.vertex.attribute(nw, vattr, unlist=FALSE)
  sapply(av, "[", 1)
}

.pop_vattrv <- function(nw, vattr){
  av <- get.vertex.attribute(nw, vattr, unlist=FALSE)
  a <- sapply(av, "[", 1)
  rest <- lapply(av, "[", -1)
  nw <- set.vertex.attribute(nw, vattr, rest)
  
  list(nw = nw, vattr = a)
}


#' A [split()] method for [`network::network`] objects.
#'
#' Split a network into subnetworks on a factor.
#'
#' @param x a [`network::network`] object.
#'
#' @param f a vector of length equal to .
#'
#' @param f,drop see [split()]; note that `f` must have length equal to `network.size(x)`.
#'
#' @return A [`network.list`] containing the networks. Note that these
#'   networks will inherit all vertex and edge attributes, as well as
#'   relevant network attributes.
#' 
#' @seealso [network::get.inducedSubgraph()]
#' @export
split.network <- function(x, f, drop = FALSE, sep = ".", lex.order = FALSE, ...) 
{
  ### NOTE: This is taken from the split.default() implementation, but is trivial.
  if(!missing(...)) 
    .NotYetUsed(deparse(...), error = FALSE)
  if(is.list(f)) 
    f <- interaction(f, drop = drop, sep = sep, lex.order = lex.order)
  else if (!is.factor(f)) 
    f <- as.factor(f)
  else if (drop) 
    f <- factor(f)
  ### END Taken from split.default().
  
  o <- lapply(levels(f), function(l) network::get.inducedSubgraph(x, which(f==l)))
  class(o) <- c("network.list", class(o))
  o
}

#' Split up a network into a list of subgraphs
#'
#' Given a network created by [combine_networks()], [uncombine_networks()] returns a list of networks,
#' preserving attributes that can be preserved.
#'
#' @param nw a [`network::network`] created by [combine_networks()].
#'
#' @param ignore.nattr,ignore.vattr,ignore.eattr network, vertex, and
#'   edge attributes not to be processed as described below.
#'
#' @param split.vattr name of the vertex attribute on which to split.
#'
#' @param detect.edgecov if `TRUE`, split up network attributes that
#'   look like dyadic covariate ([`ergm::edgecov`]) matrices.
#'
#' @param names.vattr optional name of the vertex attribute to use as network
#'   names in the output list.
#'
#' @return a list of [`network::network`]s containing subgraphs on `split.vattr`. In particular,
#'
#' * their basic properties (directedness and bipartednes) are the same as those of the input network;
#'
#' * vertex attributes are split;
#'
#' * edge attributes are assigned to their respective edges in
#'   the returned networks.
#'
#' If `split.vattr` is a vector, only the first element is used and it's "popped".
#'
#' @seealso [split.network()]
#' @examples
#'
#' data(samplk)
#' 
#' o1 <- combine_networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#'
#' ol <- uncombine_network(o1)
#'
#' @export
uncombine_network <- function(nw, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), split.vattr=".NetworkID", names.vattr=NULL){
  tmp <- .pop_vattrv(nw, split.vattr); nw <- tmp$nw; f <- tmp$vattr
  if(!is.null(names.vattr)) tmp <- .pop_vattrv(nw, names.vattr); nw <- tmp$nw; nwnames <- tmp$vattr
  nwl <- split(nw, f)

  for(a in setdiff(list.network.attributes(nw),
                   ignore.nattr)){ # I.e., iterate through common attributes.
    nwl <- mapply(set.network.attribute, x=nwl, value=get.network.attribute(nw, a, unlist=FALSE), MoreArgs=list(attrname=a), SIMPLIFY=FALSE)
  }

  if(!is.null(names.vattr)) names(nwl) <- unique(nwnames)
  
  class(nwl) <- c("network.list", class(nwl))
  nwl
}

.split_constr_network <- function(nw, split.vattr=".NetworkID", names.vattr=".NetworkName"){
  uncombine_network(nw, split.vattr=split.vattr, names.vattr=names.vattr, ignore.nattr = c(eval(formals(uncombine_network)$ignore.nattr), "constraints", "obs.constraints"))
}

  #' Calculate a vector that maps the global LHS network Vertex indices within-layer Vertex and a Vertex to layer lookup table.
#' @noRd
.block_vertexmap <- function(nw, by=".NetworkID", same_dim=FALSE){
  a <- .peek_vattrv(nw, by)
  n <- length(a)
  bip <- nw %n% "bipartite"
  if(NVL(bip,0)){
    ea <- a[seq_len(bip)]
    aa <- a[bip+seq_len(n-bip)]
    el <- rle(ea)$lengths
    al <- rle(aa)$lengths
    if(same_dim) if(!all_identical(el) || !all_identical(al)) stop("Layers must be networks of the same dimensions.", call.=FALSE)

    eoff <- rep(cumsum(c(0,el[-length(el)])), el)
    aoff <- rep(cumsum(c(0,el[-length(al)])), al) + sum(el)
    
    o <- list(nb = length(el), bids = a, bmap = seq_len(n) - c(eoff,aoff), ns = rbind(el,al))
  }else{
    l <- rle(a)$lengths
    if(same_dim) if(!all_identical(l)) stop("Layers must be networks of the same size.", call.=FALSE)
    off <- rep(cumsum(c(0,l[-length(l)])), l)
    o <- list(nb = length(l), bids = a, bmap = seq_len(n) - off, ns = l)
  }
  o
}
