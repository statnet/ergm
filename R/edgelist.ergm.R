edgelist.ergm <- function(x, ...) {
  UseMethod("edgelist.ergm")
}

edgelist.ergm.network <- function(x, ...) {
  edgelist.ergm(as.matrix.network.edgelist(x), directed=is.directed(x), ...)
}

edgelist.ergm.matrix <- function(x, directed=TRUE, ...) {
  if (NCOL(x)==2) {
    if (!directed) { # Change if necessary so [,1] < [,2]
      tmp <- apply(x, 1, function(a) a[1]>a[2]) 
      x[tmp,] <- x[tmp,2:1]
    }
    x <- x[order(x[,1], x[,2]), , drop=FALSE]
    x <- x[0!=apply(x, 1, diff), , drop=FALSE]
    return(unique(x))
  }
  else if (NROW(x)==NCOL(x)){ # Square matrix; assume adjacency
    nz <- x!=0
    if (missing(directed))
      directed <-  !all(x==t(x))
    return(edgelist.ergm(cbind(row(x)[nz], col(x)[nz]), directed=directed, ...))
  }
  else { # Assume bipartite undirected
    nz <- x!=0
    return(edgelist.ergm(cbind(row(x)[nz], NROW(x)+col(x)[nz]), directed=FALSE, ...))    
  }
}

edgelist.ergm.default <- function(x, ...) {
  if (is.null(x))
    return(NULL)
  stop("edgelist.ergm cannot convert argument to edgelist")
}

