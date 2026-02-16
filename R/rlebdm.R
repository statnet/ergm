#  File R/rlebdm.R in package ergm, part of the Statnet suite of packages for
#  network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' RLE-Compressed Boolean Dyad Matrix
#'
#' A simple class representing boolean (logical) square matrix
#' run-length encoded in a column-major order.
#'
#' @param x for [rlebdm()], an [rle()] object or a vector that is converted to one; it will be coerced to [logical()] before processing; for [as.rlebdm.matrix()], a matrix.
#' @param n the dimensions of the square matrix represented.
#'
#' @examples
#' # From a vector
#' rlebdm(rep(rep(c(0,1),each=3),14)[seq_len(81)], 9)
#'
#' # From a constant
#' rlebdm(1, 3)
#'
#' # Large matrix (overflowing .Machine$integer.max)
#' big <- rlebdm(1, 50000)
#' unclass(big) # Represented as two runs
#' big # Only summary is printed
#' stopifnot(length(big)==50000^2)
#'
#' @seealso [as.rlebdm.ergm_conlist()]
#' @import rle
#' @import statnet.common
#' @keywords internal
#' @export
rlebdm <- function(x, n){
  if(is(x, "rlebdm")) return(x)
  o <- as.rle(x)
  o$values <- as.logical(o$values)
  l <- n^2 # 2 is numeric, so it upcasts n to numeric even if n is an integer.
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
  n <- as.integer(attr(x, "n")) # network size

  x <- as.matrix(x) # matrix edgelist
  storage.mode(x) <- "integer" # be sure tails and heads are stored as integers

  o <- order(x[,2L],x[,1L]) # order edges by head, then tail
  t <- x[o,1L] # tails
  h <- x[o,2L] # heads
  
  ne <- length(t) # number of edges

  vals <- logical(n + 2L*ne) # RLE values, defaulting to FALSE
  lens <- rep(n, n + 2L*ne) # RLE lengths, defaulting to n

  ci <- 1L # current index in vals and lens
  
  lh <- 1L # last head we handled
  lt <- 0L # last tail we handled
  
  for(i in seq_len(ne)){ # for each edge
    if(h[i] == lh){ # if this head is the same as the last head
      # add the FALSE run between this tail and the last tail
      lens[ci] <- t[i] - lt - 1L
    }else{ # else this head is different than the last head
      # add the final FALSE run for the last head
      lens[ci] <- n - lt
      
      # move the index forward to correspond to the current head
      ci <- ci + h[i] - lh

      # add the initial FALSE run for the current head
      lens[ci] <- t[i] - 1L
      
      # update the last head
      lh <- h[i]
    }

    # update the index for the (last) FALSE run added in the above conditional
    ci <- ci + 1L
    
    # add the TRUE run for the current tail
    lens[ci] <- 1L
    vals[ci] <- TRUE
    ci <- ci + 1L
    
    # update the last tail
    lt <- t[i]    
  }
  
  # add the final FALSE run, if needed
  if(ne > 0){
    lens[ci] <- n - lt
  }
  
  rlebdm(compress(structure(list(values=vals, lengths=lens), class = "rle")), n)
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

#' @describeIn rlebdm
#'
#' Strip `rlebdm`-specific attributes and class, returning a plain [`rle`] object.
#'
#' @export
as.rle.rlebdm <- function(x) structure(x, class = "rle", n = NULL)

#' @rdname rlebdm
#'
#' @param compact whether to print the matrix compactly (dots and stars) or to print it as a logical matrix.
#' 
#' @export
print.rlebdm <- function(x, compact=TRUE, ...){
  if(length(x) > getOption("max.print")){
    cat(sprintf(
      "Large Run-Length-Encoded Binary Dyad Matrix:\n  dimension: %0.0f*%0.0f\n  number of 1s: %0.0f/%0.0f\n  density: %f\n",
      nrow(x), ncol(x), sum(x), length(x), sum(x)/length(x)
    ))
    return(invisible(x))
  }
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
#' @param e1,e2 arguments to the unary (`e1`) or the binary (`e1` and `e2`) operators.
#'
#' @note The arithmetic operators are mathematical functions are
#'   implemented for the [`Ops`] and the [`Math`] group generics and
#'   therefore work for almost all of them automatically. To preserve
#'   the integrity of the data structure, the results are cast to
#'   logical before return.
#'
#' @export
Ops.rlebdm <- function(e1, e2){
  o <- NextMethod()
  rlebdm(o, attr(e1, "n"))
}

#' @rdname rlebdm
#' @export
Math.rlebdm <- function(x, ...){
  o <- NextMethod()
  rlebdm(o, attr(x, "n"))
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
#' @param constraints.obs an [`ergm_conlist`] object specifying the
#'   observation process constraints; defaults to `NULL` for all dyads
#'   observed (i.e., a constant matrix of `FALSE`).
#' @param which which aspect of the constraint to extract: \describe{
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
#' @seealso [`ergmConstraint`]
#'
#' @keywords internal
#' @export
as.rlebdm.ergm_conlist <- function(x, constraints.obs = NULL, which = c("free", "missing", "informative"), ...){
  # FIXME: Probably don't need all these recursive calls.
  which <- match.arg(which)
  switch(which,
         free={
           y <- NULL
           for(con in x){
             if(!is.null(free_dyads(con))){
               y <-
                 if(is.null(y)) free_dyads(con)
                 else{
                   if(NVL(con$sign, +1)==+1) y & free_dyads(con)
                   else y | free_dyads(con)
                 }
             }
           }
           if(!is.null(y)) compress(y)
         },
         missing={
           free_dyads.obs <- as.rlebdm(constraints.obs)
           
           if(is.null(x)){
             free_dyads.obs # Already compacted.
           }else{
             NVL3(free_dyads.obs, compress(as.rlebdm(x) & .),  NULL)
           }
         },
         informative={
           y <- as.rlebdm(x)
           EVL3(constraints.obs, compress(y & !as.rlebdm(x, ., which="missing")), y)
         }
         )
}

#' @describeIn rlebdm Compress the `rle` data structure in the
#'   `rlebdm` by merging successive runs with identical values.
#' @export
compress.rlebdm <- function(x, ...){
  y <- NextMethod()
  structure(y, n=attr(x, "n"), class=class(x))
}

#' @describeIn rlebdm
#'
#' Convert an [`rlebdm`] object to an [`edgelist`]: a two-column
#' integer matrix or [`tibble`] giving the cells with `TRUE` values.
#'
#' @param prototype an optional network with network attributes that
#'   are transferred to the edgelist and will filter it (e.g., if the
#'   prototype network is given and does not allow self-loops, the
#'   edgelist will not have self-loops either,e ven if the dyad matrix
#'   has non-`FALSE` diagonal).
#' @param output a string specifying whether the result should be a
#'   matrix or a [`tibble`].
#' @importFrom network as.edgelist
#' @seealso [as.edgelist()]
#' @export
as.edgelist.rlebdm <- function(x, prototype=NULL, ..., output=c("matrix", "tibble")){
  output <- match.arg(output)

  dir <- NVL3(prototype, is.directed(.), TRUE)
  loop <- NVL3(prototype, has.loops(.), TRUE)
  bip <- NVL3(prototype, b1.size(.), FALSE)

  n <- attr(x, "n")
  starts <- cumsum(c(1,as.numeric(x$lengths)))
  starts <- starts[-length(starts)]

  starts <- starts[x$values!=0]
  lengths <- x$lengths[x$values!=0]
  ends <- starts + lengths - 1
  values <- x$values[x$values!=0]
  
  d <- do.call(rbind,
               Map(function(s,e,v){
                 cbind(s:e,v)
               }, starts, ends, values))

  el <- tibble(.tail = as.integer((d[,1L]-1L) %% n + 1L),
               .head = as.integer((d[,1L]-1L) %/% n + 1L))

  if(!is.logical(values)) el <- cbind(el, .values = d[, 2L])

  if(!dir) el <- el[el[[1L]]<=el[[2L]], , drop=FALSE]
  if(!loop) el <- el[el[[1L]]!=el[[2L]], , drop=FALSE]
  if(bip) el <- el[(el[[1L]]<=bip) != (el[[2L]]<=bip), , drop=FALSE]

  if(!is.null(prototype)){
    attr(el, "directed") <- dir
    attr(el, "bipartite") <- bip
    attr(el, "loops") <- loop
  }
    
  attr(el, "n") <- n

  if(output=="matrix") el <- as.matrix(el)

  class(el) <- c(paste0(output, "_edgelist"), "edgelist", class(el))
  el
}

#' @include to_ergm_Cdouble.R
#'
#' @describeIn to_ergm_Cdouble
#'
#' Method for [`rlebdm`] objects.
#'
#' @return The [`rlebdm`] method returns a vector with the following:
#' * number of nonzero dyads,
#' * number of runs of nonzeros,
#' * starting positions of the runs, and
#' * cumulative lengths of the runs, prepended with 0.
#' @export
to_ergm_Cdouble.rlebdm <- function(x, ...){
  x <- compress(x) # Just in case.
  cumlen <- cumsum(as.numeric(x$lengths[x$values==TRUE]))
  nruns <- length(cumlen)
  ndyads <- cumlen[nruns]
  runstart <- cumsum(c(1,as.numeric(x$lengths)))
  runstart <- runstart[-length(runstart)]
  runstart <- runstart[x$values==TRUE]

  as.double(c(attr(x,"n"), ndyads, nruns, runstart, c(0,cumlen)))
}
