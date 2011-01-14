#=========================================================================
# This file contains the following 4 functions for converting various
# objects into edgelists:
#        <edgelist.ergm>
#        <edgelist.ergm.default>
#        <edgelist.ergm.network>
#        <edgelist.ergm.matrix>
#=========================================================================






##########################################################################
# Each of the <edgelist.ergm.x> functions converts an object x into an
# edgelist.
#
# --PARAMETERS--
#   x   :  either a network, a matrix, or null
#   ... :  additional paramters which may include:
#      directed        : whether the edgelist shall be directed (T or F);
#                        default=TRUE;
#      check.uniqueness: whether edges must be unique (T or F);
#                        default=TRUE;
#      check.sorted    : whether to sort the edgelist (T or F);
#                        default=TRUE;
#
# --RETURNED--
#   the edge list corresponding to x under the following assumptions about
#   x as a matrix:
#    (i) if x is originally 2 columns, it was already an edgelist
#    (2) if x is originally square, it is a non-bipartite adjacency matrix
#    (3) if x is originally rectangular, it is a bipartite adjacency matrix
#
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
  edgelist.ergm(as.matrix.network.edgelist(x), directed=is.directed(x), ...)
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



