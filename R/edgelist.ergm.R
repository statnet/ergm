#  File ergm/R/edgelist.ergm.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##########################################################################
# Each of the <edgelist.ergm.x> functions converts an object x into an
# edgelist.
##########################################################################

edgelist.ergm <- function(x, ...) {
  UseMethod("edgelist.ergm")
}


edgelist.ergm.default <- function(x, ...) {
  if (is.null(x))
    return(NULL)
  stop("edgelist.ergm cannot convert argument to edgelist")
}



edgelist.ergm.network <- function(x, ...) {
  edgelist.ergm(as.edgelist(x), directed=is.directed(x), ...)
}



edgelist.ergm.matrix <- function(x, directed=TRUE, check.uniqueness=TRUE, 
                                 check.sorted=TRUE, ...) {
  if (NCOL(x)==2) { # Assume an edgelist
    if (!directed) { # Change if necessary so [,1] < [,2]
      tmp <- apply(x, 1, function(a) a[1]>a[2]) 
      x[tmp,] <- x[tmp,2:1]
    }
    if (check.sorted) {
      x <- x[order(x[,1], x[,2]), , drop=FALSE]
    }
    x <- x[x[,1]!=x[,2], , drop=FALSE] # delete self-loops no matter what
    if (check.uniqueness)
      x <- unique(x)
    return(x)
  }
  else if (NROW(x)==NCOL(x)){ # Square matrix; assume adjacency
    nz <- x!=0
    if (missing(directed))
      directed <-  !all(x==t(x))
    if (!directed) 
      nz <- nz & (row(x)<col(x)) # upper triangle
    return(cbind(col(x)[t(nz)], row(x)[t(nz)])) 
    # Note:  We did not use cbind(row(x)[nz], col(x)[nz])
    # so output would be assured to have unique rows in sorted order
  }
  else { # Assume bipartite undirected
    nz <- x!=0
    return(cbind(NROW(x)+col(x)[t(nz)], row(x)[t(nz)]))
    # Note:  We did not use cbind(NROW(x)+row(x)[nz], col(x)[nz])
    # so output would be assured to have unique rows in sorted order
  }
}



