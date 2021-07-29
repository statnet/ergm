#  File R/InitErgmConstraint.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################

# Baseline constraint incorporating network attributes such as
# directedness, bipartitedness, and self-loops.
InitErgmConstraint..attributes <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)

  n <- network.size(nw)
  storage.mode(n) <- "integer"
  dir <- is.directed(nw)
  loops <- has.loops(nw)
  bip <- EVL(as.integer(nw%n%"bipartite"), FALSE)
  rm(nw) # All needed information has now been extracted from nw.

  list(
    free_dyads = function(){
      ## NB: Free dyad RLE matrix is stored in a column-major order for
      ## consistency with R.
      d <-
        if(dir){
          if(loops){
            compress(structure(list(lengths=rep(n,n), values=rep(TRUE,n)), class="rle"))
          }else{
            structure(list(lengths=c(1L,rep(c(n,1L),n-1L)), values=c(rep(c(FALSE, TRUE),n-1L),FALSE)), class="rle")
          }
        }else if(bip){
          b1 <- as.integer(bip)
          b2 <- n - b1
          compress(structure(list(lengths=c(rep(n,b1), rep(c(b1,b2),b2)), values=c(rep(FALSE, b1), rep(c(TRUE,FALSE),b2))),class="rle"))
        }else{
          if(loops){
            vals <- c(rep(c(TRUE,FALSE),n-1L),TRUE)
            lens <- integer(2L*(n-1L)+1L)
            for(i in seq_len(n-1L)){
              lens[2L*i-1L] <- i
              lens[2L*i] <- n-i
            }
            lens[2L*n-1L] <- n
          }else{
            vals <- c(rep(c(FALSE,TRUE),n-1L),FALSE)
            lens <- integer(2L*(n-1L)+1L)
            for(i in seq_len(n-1L)){
              lens[2L*i-1L] <- n-i+1L
              lens[2L*i] <- i
            }
            lens[2L*n-1L] <- 1L
          }          
          structure(list(lengths=lens,values=vals), class="rle")
        }
      rlebdm(d, n)
    },
    implies = ".attributes",
    dependence = FALSE)
}

InitErgmConstraint.edges<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, implies = "edges")
}

InitErgmConstraint.degrees<-InitErgmConstraint.nodedegrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, constrain = "degrees", implies = c("degrees", "edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
}

InitErgmConstraint.odegrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("odegrees", "edges", "odegreedist"))
}

InitErgmConstraint.idegrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("idegrees", "edges", "idegreedist"))
}

InitErgmConstraint.b1degrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE)
  list(dependence = TRUE, implies = c("b1degrees", "edges"))
}

InitErgmConstraint.b2degrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE)
  list(dependence = TRUE, implies = c("b2degrees", "edges"))
}

InitErgmConstraint.degreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, implies = c("degreedist", "edges", "idegreedist", "odegreedist"))
}

InitErgmConstraint.idegreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("idegreedist", "edges"))
}

InitErgmConstraint.odegreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("odegreedist", "edges"))
}

InitErgmConstraint.bd<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attribs", "maxout", "maxin", "minout", "minin"),
                      vartypes = c("matrix", "numeric,matrix", "numeric,matrix", "numeric,matrix", "numeric,matrix"),
                      defaultvalues = list(NULL, NA_integer_, NA_integer_, NA_integer_, NA_integer_),
                      required = c(FALSE, FALSE, FALSE, FALSE, FALSE))

   if(is.na(a$minout) && is.na(a$minin)) {
     constrain <- c("bd","bdmax")
   } else {
     constrain <- "bd"
   }

   list(constrain=constrain, attribs=a$attribs, maxout=a$maxout, maxin=a$maxin, minout=a$minout, minin=a$minin)
}

InitErgmConstraint.blocks <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "b1levels", "b2levels", "levels", "levels2"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                      defaultvalues = list(NULL, NULL, NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))
  attr <- a$attr; b1levels <- a$b1levels; b2levels <- a$b2levels; levels <- a$levels; levels2 <- a$levels2

  if(is.bipartite(nw)) {
    b1nodecov <- ergm_get_vattr(attr, nw, bip = "b1")
    b2nodecov <- ergm_get_vattr(attr, nw, bip = "b2")
    
    b1namescov <- ergm_attr_levels(b1levels, b1nodecov, nw, sort(unique(b1nodecov)))
    b2namescov <- ergm_attr_levels(b2levels, b2nodecov, nw, sort(unique(b2nodecov)))
    
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    
    levels2.list <- transpose(expand.grid(row = b1namescov, col = b2namescov, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = nr + 1:nc)
   
    levels2.sel <- ergm_attr_levels(levels2, list(row = b1nodecov, col = b2nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
  
    b1nodecov <- match(b1nodecov, b1namescov, nomatch = length(b1namescov) + 1)
    b2nodecov <- match(b2nodecov, b2namescov, nomatch = length(b2namescov) + 1)
      
    nodecov <- c(b1nodecov, b2nodecov)
                                           
    u[,2L] <- u[,2L] - nr
    amat <- matrix(TRUE, nrow = nr + 1, ncol = nc + 1)
    amat[as.matrix(u)] <- FALSE
    
    row_nodecov <- b1nodecov
    col_nodecov <- b2nodecov
    
  } else {
    nodecov <- ergm_get_vattr(attr, nw)
  
    u <- ergm_attr_levels(levels, nodecov, nw, sort(unique(nodecov)))
    namescov <- u 
    
    nr <- length(u)
    nc <- length(u)

    levels2.list <- transpose(expand.grid(row = u, col = u, stringsAsFactors=FALSE))
    indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
    uun <- as.vector(outer(u,u,paste,sep="."))
    
    if(!is.directed(nw)) {
        rowleqcol <- indices2.grid$row <= indices2.grid$col
        levels2.list <- levels2.list[rowleqcol]
        indices2.grid <- indices2.grid[rowleqcol,]
        uun <- uun[rowleqcol]
    }    
   
    levels2.sel <- ergm_attr_levels(levels2, list(row = nodecov, col = nodecov), nw, levels2.list)
    
    rows2keep <- match(levels2.sel,levels2.list, NA)
    rows2keep <- rows2keep[!is.na(rows2keep)]
  
    u <- indices2.grid[rows2keep,]
    uun <- uun[rows2keep]

    nodecov <- match(nodecov, namescov, nomatch = length(namescov) + 1)
    
    amat <- matrix(TRUE, nrow = nr + 1, ncol = nc + 1)
    amat[as.matrix(u)] <- FALSE
    if(!is.directed(nw)) amat <- amat & t(amat)
    
    row_nodecov <- nodecov
    col_nodecov <- nodecov
    
  }  

  constrain <- "blocks"
  
  n <- as.integer(network.size(nw))

  if(is.bipartite(nw)) {
    b1 <- as.integer(nw %n% "bipartite")
    b2 <- n - b1
  } else {
    b1 <- 0L
    b2 <- 0L
  }
  
  rm(nw) # All needed information has now been extracted from nw.

  free_dyads <- function() {
    rle_list <- vector(mode = "list", length = nc + 1)
    for(i in seq_len(nc + 1)) {
      rle_list[[i]] <- rle(c(amat[row_nodecov,i], logical(b2)))
    }

    lens <- vector(mode = "list", length = n)
    vals <- vector(mode = "list", length = n)

    for(i in seq_len(n)) {
      if(i > b1) {
        lens[[i]] <- rle_list[[col_nodecov[i - b1]]]$lengths
        vals[[i]] <- rle_list[[col_nodecov[i - b1]]]$values
      } else {
        lens[[i]] <- n
        vals[[i]] <- FALSE
      }
    }

    rlebdm(compress(structure(list(lengths=unlist(lens), values=unlist(vals)), class="rle")), n)
  }

  list(constrain = constrain,
       dependence = FALSE, 
       free_dyads = free_dyads,
       nodecov = nodecov,
       amat = amat)
}


InitErgmConstraint.hamming<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE)
}

InitErgmConstraint.observed <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(free_dyads = as.rlebdm(as.edgelist(is.na(nw))),
       dependence = FALSE, implies = c("observed"))
}

InitErgmConstraint.fixedas<-function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("present", "absent"),
                      vartypes = c("network,matrix", "network,matrix"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  present <- a$present; absent <- a$absent
  if(is.null(present) && is.null(absent))
    ergm_Init_abort(paste("fixedas constraint takes at least one argument, either present or absent or both."))

  list(
    free_dyads = function(){
      if(is.network(present)) present <- as.edgelist(present)
      if(is.network(absent)) absent <- as.edgelist(absent)

      # FixedEdgeList
      fixed <- as.edgelist(rbind(present,absent),
                           n=nw%n%"n",
                           directed=nw%n%"directed",
                           bipartite=nw%n%"bipartite",
                           loops=nw%n%"loops")
      if(any(duplicated(fixed))){
        ergm_Init_abort("Dyads cannot be fixed at both present and absent")
      }

      !as.rlebdm(fixed)
    },
    dependence = FALSE)
}


InitErgmConstraint.fixallbut<-function(nw, arglist,...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("free.dyads"),
                      vartypes = c("network,matrix"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  free.dyads <- a$free.dyads

  list(
    free_dyads = function(){
      if(is.network(free.dyads)) free.dyads <- as.edgelist(free.dyads)
      fixed <- as.edgelist(free.dyads,
                           n=nw%n%"n",
                           directed=nw%n%"directed",
                           bipartite=nw%n%"bipartite",
                           loops=nw%n%"loops")
      as.rlebdm(free.dyads)
    },
    dependence = FALSE)
}

InitErgmConstraint.dyadnoise<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("p01", "p10"),
                      vartypes = c("numeric,matrix", "numeric,matrix"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))
  p01 <- a$p01; p10 <- a$p10

  if(((length(p01) != 1 || length(p10) != 1) &&
      any(dim(as.matrix(nw, matrix.type="adjacency")) != c(dim(p01),dim(p10))))) # FIXME: Don't create an adjacency matrix unnecessarily.
    stop("p01 and p10 must be either scalars or matrices of the same dimension as the adjacency matrices of the LHS network.")

  list(p01=p01, p10=p10)
}

InitErgmConstraint.egocentric <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "direction"),
                      vartypes = c(ERGM_VATTR_SPEC, "character"),
                      defaultvalues = list(NULL, "both"),
                      required = c(FALSE, FALSE))
  attr <- a$attr; direction <- a$direction

  direction <- match.arg(direction, c("both", "out", "in"))
  if(!is.directed(nw) && direction!="both")
    stop("Directed egocentric constraint cannot be used for an undirected network.")

  list(
    free_dyads = function(){
      n <- network.size(nw)
      a <- ( # Are that node's dyads toggleable?
        if(is.null(attr)) get.vertex.attribute(nw, "na")
        else !as.vector(ergm_get_vattr(attr, nw, accept="logical"))
      )

      # Remember: column-major order.

      rlea <- rle(a)
      
      fd <- rlebdm(switch(direction,
                          `out` = rep(rlea, n),
                          `in` = rep(rlea, rep(n, length(rlea$lengths)), scale="run"),
                          `both` = compress(rep(rlea, n) & rep(rlea, rep(n, length(rlea$lengths)), scale="run"))), # The others are already compressed by rep().
                   n)
    },
    dependence = FALSE
  )
}

InitErgmConstraint.Dyads<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("fix", "vary"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  fix <- a$fix; vary <- a$vary

  if(is.null(fix) & is.null(vary))
    ergm_Init_abort(paste("Dyads constraint takes at least one argument, either",sQuote("fix"),"or",sQuote("vary"),"or both."))

  list(
    free_dyads = function(){
      fd <- lapply(list(fix=fix,vary=vary),
                   function(f){
                     if(!is.null(f)){
                       f[[3]] <- f[[2]]
                       f[[2]] <- nw
                       m <- ergmMPLE(f, expand.bipartite=TRUE, output="array")$predictor
                       m <- m!=0
                       m[is.na(m)] <- FALSE
                       if(!is.directed(nw)){
                         m <- m | aperm(m, c(2L,1L,3L))
                       }
                       lapply(seq_len(dim(m)[3]), function(i) as.rlebdm(m[,,i]))
                     }
                   })
      fd$fix <- if(length(fd$fix)) fd$fix %>% map(`!`) %>% reduce(`&`)
      fd$vary <- if(length(fd$vary)) fd$vary %>% reduce(`|`)
      fd <- Reduce(`|`, fd)

      compress(fd)
    },
    dependence = FALSE
  )
}
