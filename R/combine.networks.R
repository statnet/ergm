#  File R/combine.networks.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

#' Test if all items in a vector or a list are identical.
#'
#' @param x a vector or a list
#'
#' @return `TRUE` if all elements of `x` are identical to each other.
#'
#' @seealso [base::identical()]
#'
#' @examples
#'
#' stopifnot(!all.same(1:3))
#'
#' stopifnot(all.same(list("a", "a", "a")))

all.same <- function(x){
  if(length(x)==0) return(TRUE)
  v0 <- x[1]
  for(v in x[-1]) if(!identical(v0,v)) return(FALSE)
  return(TRUE)
}

#' A single block-diagonal network created by combining multiple networks
#'
#' Given a list of compatible networks, the function returns a single
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
#' * the returned network's size is the sum of the input networks;
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
#' o1 <- combine.networks(list(samplk1, samplk2, samplk3))
#' image(as.matrix(o1))
#' head(get.vertex.attribute(o1, ".NetworkID"))
#' o2 <- combine.networks(list(o1, o1))
#' image(as.matrix(o2))
#' head(get.vertex.attribute(o2, ".NetworkID", unlist=FALSE))
#'
#' data(florentine)
#' f1 <- combine.networks(list(business=flobusiness, marriage=flomarriage), blockName.vattr=".NetworkName") 
#' image(as.matrix(f1))
#' head(get.vertex.attribute(f1, ".NetworkID"))
#' head(get.vertex.attribute(f1, ".NetworkName"))
#' @export
combine.networks <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
  if(any(sapply(nwl, is.bipartite))) .combine.networks.bipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, standardized=standardized, keep.unshared.attr=keep.unshared.attr)
  else .combine.networks.unipartite(nwl=nwl, ignore.nattr=ignore.nattr, ignore.vattr=ignore.vattr, ignore.eattr=ignore.eattr, blockID.vattr=blockID.vattr, blockName.vattr=blockName.vattr, detect.edgecov=detect.edgecov, standardized=standardized, keep.unshared.attr=keep.unshared.attr)
}


.combine.networks.unipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
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
       && all.same(sapply(vl, mode))){

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


.combine.networks.bipartite <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockID.vattr=".NetworkID", blockName.vattr=NULL, detect.edgecov=FALSE, standardized=FALSE, keep.unshared.attr=FALSE){
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
       && all.same(sapply(vl, mode))){

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
