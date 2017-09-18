
.rle_dyad_matrix_from_el <- function(n, el, el_free){
  ## NB: Free dyad RLE matrix is stored in a column-major order for
  ## consistency with R.
  
  ils <- lapply(seq_len(n), function(j) fel[fel[,2]==j,1])
  o <- lapply(ils, function(il){
    # If el_free, construct an rle of the form c(FALSE, TRUE, FALSE,
    # ..., FALSE, TRUE, FALSE); otherwise, construct its logical
    # negation.
    o <- rle(c(rep(c(!el_free,el_free), length(il)),!el_free))
      
    # Construct repetition counts: gaps between the is', as well as
    # the gap before the first i and after the last i for that j,
    # and interleave it with 1s.
    lens <- c(rbind(diff(c(0,il,n+1))-1,1))
    lens <- lens[-length(lens)]
    rep(o, lens, scale='run')
  })
  # Concatenate the RLEs and compact.
  compact_rle(do.call(c, o))  
}

get.free.dyads <- function(constraints){
  y <- NULL
  for(con in constraints){
    if(!is.null(con$free.dyads)){
      y <- if(is.null(y)) con$free.dyads() else y & con$free.dyads()
    }
  }
  if(!is.null(y)) compact.rle(y)
}

get.miss.dyads <- function(constraints, constraints.obs){
# Returns an RLE dyad matrix indicating the missing dyads in the
# network (respecting the constraints).
  free.dyads <- get.free.dyads(constraints)
  free.dyads.obs <- get.free.dyads(constraints.obs)
  
  if(is.null(free.dyads)){
    if(is.null(free.dyads.obs)) NULL
    else free.dyads.obs
  }else{
    if(is.null(free.dyads.obs)) !free.dyads
    else (!free.dyads) | free.dyads.obs
  }
}

get.active.dyads <- function(constraints, constraints.obs){
  get.free.dyads(constraints) & ! get.miss.dyads(constraints,constraints.obs)
}

as.edgelist.rle <- function(x, n){
  starts <- cumsum(1,as.numeric(x$lengths))
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

.dyadcount.dyadrle <- function(x){
  Tlens <- x$lengths[x$values==TRUE]
  sum(as.numeric(Tlens))
}

pack_free.dyads_as_numeric <- function(fdrle){
  cumlen <- cumsum(as.numeric(fdrle$lengths[fdrle$values==TRUE]))
  nruns <- length(cumlen)
  ndyads <- cumlen[nruns]
  runstart <- cumsum(c(1,as.numeric(fdrle$lengths)))
  runstart <- runstart[-length(runstart)]
  runstart <- runstart[fdrle$values==TRUE]

  c(ndyads, nruns, runstart, c(0,cumlen))
}
