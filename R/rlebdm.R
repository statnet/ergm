#  File R/rlebdm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################
#' RLE-Compressed Boolean Dyad Matrix
#'
#' A simple class representing boolean (logical) square matrix
#' run-length encoded in a column-major order.
#'
#' @param x for [rlebdm()], an [rle()] object or a vector that is converted to one; it will be coerced to [logical()] before processing; for [as.rlebdm.matrix()], a matrix.
#' @param n the dimensions of the square matrix represented.
#'
#' @seealso [as.rlebdm.ergm_conlist()]
#' @import statnet.common
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
#' @param ... additional arguments, currently unused.
#' @export
as.rlebdm <- function(x, ...) UseMethod("as.rlebdm")

#' @noRd
#' @export
as.rlebdm.rlebdm <- function(x, ...) x


#' @noRd
#' @export
as.rlebdm.NULL <- function(x, ...) NULL

#' @describeIn rlebdm
#'
#' Convert a square matrix of mode coercible to [`logical`] to an
#' [`rlebdm`].
#' 
#' @export
as.rlebdm.matrix <- function(x, ...){
  if(nrow(x)!=ncol(x)) stop("Input matrix must be square at this time.")
  rlebdm(c(x), nrow(x))
}

#' @describeIn rlebdm
#'
#' Convert an object of class [`edgelist`] to an [`rlebdm`] object
#' whose cells in the edge list are set to `TRUE` and whose other
#' cells are set to `FALSE`.
#' 
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

#' @describeIn rlebdm
#'
#' Convert an object of class [`network`] to an [`rlebdm`] object
#' whose cells corresponding to extant edges are set to `TRUE` and
#' whose other cells are set to `FALSE`.
#' 
#' @export
as.rlebdm.network <- function(x, ...){
  as.rlebdm(as.edgelist(x,...))
}


#' @rdname rlebdm
#' @export
as.matrix.rlebdm <- function(x, ...){
  matrix(inverse.rle(x), attr(x, "n"), attr(x, "n"))
}

#' @rdname rlebdm
#' @export
dim.rlebdm <- function(x){
  rep(attr(x, "n"),2)
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

#' Extract dyad-level ERGM constraint information into an [`rlebdm`] object
#'
#' A function to combine the `free_dyads` attributes of the
#' constraints appropriately to generate an [`rlebdm`] of dyads
#' toggleable and/or missing and/or informative under that combination
#' of constraints.
#'
#' @param x an [`ergm_conlist`] object: a list of initialised
#'   constraints. `NULL` is treated as a placeholder for no constraint
#'   (i.e., a constant matrix of `TRUE`).
#' @param constraints.obs observation process constraints; defaults to
#'   `NULL` for all dyads observed (i.e., a constant matrix of
#'   `FALSE`).
#' @param which which aspect of the constraint to extract:
#' \describe{
#' 
#' \item{`free`}{for dyads that *may* be toggled under the constraints
#'   `x`; ignores `constraints.obs`;}
#' 
#' \item{`missing`}{for dyads that are free but considered unobserved
#'   under the constraints; and}
#' 
#' \item{`informative`}{for dyads that are both free and observed.}
#' 
#' }
#' @param ... additional arguments, currently unused.
#'
#' @note For `which=="free"` or `"informative"`, `NULL` return value
#'   is a placeholder for a matrix of `TRUE`, whereas for
#'   `which=="missing"` it is a placeholder for a matrix of `FALSE`.
#' @note Each element in the constraint list has a sign, which
#'   determins whether the constraint further restricts (for `+`) or
#'   potentially relaxes restriction (for `-`).
#'
#' @seealso [ergm-constraints]
#'
#' @export
as.rlebdm.ergm_conlist <- function(x, constraints.obs = NULL, which = c("free", "missing", "informative"), ...){
  # FIXME: Probably don't need all these recursive calls.
  which <- match.arg(which)
  switch(which,
         free={
           y <- NULL
           for(con in x){
             if(!is.null(con$free_dyads)){
               y <-
                 if(is.null(y)) con$free_dyads
                 else{
                   if(NVL(con$sign, +1)==+1) y & con$free_dyads
                   else y | con$free_dyads
                 }
             }
           }
           if(!is.null(y)) rlebdm(compact.rle(y), sqrt(length(y)))
         },
         missing={
           free_dyads <- as.rlebdm(x)
           free_dyads.obs <- as.rlebdm(constraints.obs)
           
           if(is.null(free_dyads)){
             free_dyads.obs
           }else{
             NVL3(free_dyads.obs, free_dyads & .,  NULL)
           }
         },
         informative={
           y <- as.rlebdm(x)
           NVL3(constraints.obs, y & !as.rlebdm(x, ., which="missing"), y)
         }
         )
}

#' @describeIn rlebdm
#'
#' Convert an [`rlebdm`] object to an [`edgelist`]: a two-column
#' integer matrix giving the cells with `TRUE` values.
#'
#' @param prototype an optional network with network attributes that
#'   are transferred to the edgelist and will filter it (e.g., if the
#'   prototype network is given and does not allow self-loops, the
#'   edgelist will not have self-loops either,e ven if the dyad matrix
#'   has non-`FALSE` diagonal).
#' @importFrom network as.edgelist
#' @export
as.edgelist.rlebdm <- function(x, prototype=NULL, ...){
  dir <- NVL3(prototype, is.directed(.), TRUE)
  loop <- NVL3(prototype, has.loops(.), TRUE)
  bip <- NVL3(prototype, . %n% "bipartite", FALSE)
  
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

  if(!dir) el <- el[el[,1]<=el[,2],]
  if(!loop) el <- el[el[,1]!=el[,2],]
  if(bip) el <- el[(el[,1]<=bip) != (el[,2]<=bip),]

  if(!is.null(prototype)){
    attr(el, "directed") <- dir
    attr(el, "bipartite") <- bip
    attr(el, "loops") <- loop
  }
    
  attr(el, "n") <- n
  class(el) <- c("edgelist","matrix")
  el
}

#' @describeIn to_ergm_Cdouble
#'
#' Method for [`rlebdm`] objects.
#'
#' @return The [`rlebdm`] method returns a vector with the following:
#' * number of nonzero dyads,
#' * number of runs of nonzeros,
#' * starting positions of the runs, and
#' * cumulative lenght of the runs, prepended by 0.
to_ergm_Cdouble.rlebdm <- function(x, ...){
  cumlen <- cumsum(as.numeric(x$lengths[x$values==TRUE]))
  nruns <- length(cumlen)
  ndyads <- cumlen[nruns]
  runstart <- cumsum(c(1,as.numeric(x$lengths)))
  runstart <- runstart[-length(runstart)]
  runstart <- runstart[x$values==TRUE]

  as.double(c(ndyads, nruns, runstart, c(0,cumlen)))
}
