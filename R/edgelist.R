edgelist <- function(x, ...) {
  UseMethod("edgelist")
}

edgelist.network <- function(x, ...) {
  edgelist(as.matrix.network.edgelist(x), directed=is.directed(x), ...)
}

edgelist.matrix <- function(x, directed=TRUE, ...) {
  if (NCOL(x)==2) {
    if (!directed) { # Change if necessary so [,1] < [,2]
      tmp <- apply(x, 1, function(a) a[1]>a[2]) 
      x[tmp,] <- x[tmp,2:1]
    }
    x <- x[order(x[,1], x[,2]), , drop=FALSE]
    x <- x[0!=apply(x, 1, diff), , drop=FALSE]
    return(unique(x))
  }
  else if (NROW(x)==NCOL(x)){
    nz <- x!=0
    if (missing(directed))
      directed <-  !all(x==t(x))
    return(edgelist(cbind(row(x)[nz], col(x)[nz]), directed=directed, ...))
  }
  stop("edgelist cannot convert argument to edgelist")
}

edgelist.default <- function(x, ...) {
  if (is.null(x))
    return(NULL)
  stop("edgelist cannot convert argument to edgelist")
}

