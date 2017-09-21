#' RLE-Compressed Boolean Dyad Matrix
#'
#' A simple class representing boolean (logical) square matrix
#' run-length encoded in a column-major order.
#'
#' @param x for [rlebdm()], an [rle()] object or a vector that is converted to one; it will be coerced to [logical()] before processing; for [as.rlebdm.matrix()], a matrix.
#' @param n the dimensions of the square matrix represented.
#' 
#' @export
rlebdm <- function(x, n){
  if(is(x, "rlebdm")) return(x)
  o <- as.rle(x)
  o$values <- as.logical(o$values)
  l <- n*n
  if(length(o)!=l){
    if(length(o)!=1) stop("Populating the matrix can only be done with a constant value at this time.")
    o <- rep(o, l, scale="run")
  }
  attr(o, "n") <- n
  class(o) <- c("rlebdm", class(o))
  o
}

#' @rdname rlebdm
#' @export
as.rlebdm <- function(x, ...) UseMethod("as.rlebdm")

#' @noRd
#' @export
as.rlebdm.NULL <- function(x, ...) NULL

#' @rdname rlebdm
#' @export
as.rlebdm.matrix <- function(x, ...){
  if(nrow(x)!=ncol(x)) stop("Input matrix must be square at this time.")
  rlebdm(x, nrow(x))
}

#' @rdname rlebdm
#' @export
as.rlebdm.edgelist <- function(x, ...){
  n <- attr(x, "n")
  ils <- lapply(lapply(lapply(seq_len(n), function(j) x[x[,2]==j,1]), unique), sort)
  o <- lapply(ils, function(il){
    o <- rle(c(rep(c(FALSE,TRUE), length(il)),FALSE))
    
    # Construct repetition counts: gaps between the i's, as well as
    # the gap before the first i and after the last i for that j,
    # and interleave it with 1s.
    lens <- c(rbind(diff(c(0,il,n+1))-1,1))
    lens <- lens[-length(lens)]
    rep(o, lens, scale='run')
  })
  # Concatenate the RLEs and compact.
  rlebdm(compact.rle(do.call(c, o)), attr(x, "n"))
}

#' @rdname rlebdm
#' @export
as.matrix.rlebdm <- function(x, ...){
  matrix(inverse.rle(x), attr(x, "n"), attr(x, "n"))
}

#' @rdname rlebdm
#'
#' @param compact whether to print the matrix compactly (dots and stars) or to print it as a logical matrix.
#' 
#' @export
print.rlebdm <- function(x, compact=TRUE, ...){
  x <- as.matrix(x)
  if(compact){
    x <- ifelse(x, "*", ".")
    x <- paste0(apply(x, 1, paste0, collapse=""),collapse="\n")
    cat(x,"\n")
  }else{
    print(x, ...)
  }
}

#' @rdname rlebdm
#' @export
`!.rlebdm` <- function(x){
  o <- NextMethod()
  rlebdm(o, attr(x, "n"))
}

#' @rdname rlebdm
#' @param e1,e2 arguments to the binary operations.
#' @export
`|.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`&.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`<.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`>.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`<=.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`>=.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`==.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
`!=.rlebdm` <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @noRd
as.rlebdm.conlist <- function(x, constraints.obs = NULL, which = c("free", "missing", "active"), ...){
  # FIXME: Probably don't need all these recursive calls.
  which <- match.arg(which)
  switch(which,
         free={
           y <- NULL
           for(con in x){
             if(!is.null(con$free_dyads)){
               y <- if(is.null(y)) con$free_dyads else y & con$free_dyads
             }
           }
           if(!is.null(y)) rlebdm(compact.rle(y), sqrt(length(y)))
         },
         missing={
           # Returns an RLE dyad matrix indicating the missing dyads in the
           # network (respecting the constraints).
           free_dyads <- as.rlebdm(x)
           free_dyads.obs <- as.rlebdm(constraints.obs)
           
           if(is.null(free_dyads)){
             if(is.null(free_dyads.obs)) NULL
             else free_dyads.obs
           }else{
             if(is.null(free_dyads.obs)) !free_dyads
             else (!free_dyads) | free_dyads.obs
           }
         },
         active={
           y <- as.rlebdm(x)
           if(is.null(constraints.obs)) y else y & !as.rlebdm(x,constraints.obs, which="missing")
         }
         )
}

#' @rdname rlebdm
as.edgelist.rlebdm <- function(x, ...){
  n <- attr(x, "n")
  starts <- cumsum(c(1,as.numeric(x$lengths)))
  starts <- starts[-length(starts)]

  starts <- starts[x$values!=0]
  lengths <- x$lengths[x$values!=0]
  ends <- starts + lengths - 1
  values <- x$values[x$values!=0]
  
  d <- do.call(rbind,
               mapply(function(s,e,v){
                 cbind(s:e,v)
               }, starts, ends, values, SIMPLIFY=FALSE))
  
  el <- cbind((d[,1]-1) %% n + 1, (d[,1]-1) %/% n + 1)
  if(!is.logical(values)) el <- cbind(el, d[,2])
  attr(el, "n") <- n
  el
}

pack_rlebdm_as_numeric <- function(x){
  cumlen <- cumsum(as.numeric(x$lengths[x$values==TRUE]))
  nruns <- length(cumlen)
  ndyads <- cumlen[nruns]
  runstart <- cumsum(c(1,as.numeric(x$lengths)))
  runstart <- runstart[-length(runstart)]
  runstart <- runstart[x$values==TRUE]

  c(ndyads, nruns, runstart, c(0,cumlen))
}
