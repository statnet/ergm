#  File R/InitErgmConstraint.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################

# Meta-constraint for a dot placeholder
InitErgmConstraint.. <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = FALSE)
}

# Meta-constraint selecting a specific proposal.
InitErgmConstraint..select <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("proposal"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  list(dependence = TRUE, proposal = a$proposal)
}

# Baseline constraint incorporating network attributes such as
# directedness, bipartitedness, and self-loops.
InitErgmConstraint..attributes <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)

  n <- network.size(nw)
  storage.mode(n) <- "integer"
  dir <- is.directed(nw)
  loops <- has.loops(nw)
  bip <- EVL(as.integer(nw%n%"bipartite"), FALSE)
  rm(nw, arglist, "...") # All needed information has now been extracted.

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

#' @templateVar name .dyads
#' @title A meta-constraint indicating handling of arbitrary dyadic constraints
#' @description This is a flag in the proposal table indicating that the proposal can enforce arbitrary combinations of dyadic constraints. It cannot be invoked directly by the user.
#'
#' @template ergmConstraint-general
NULL

#' @templateVar name edges
#' @title Preserve the edge count of the given network
#' @description Only networks
#'   having the same number of edges as the network passed
#'   in the model formula have non-zero probability.
#'
#' @usage
#' # edges
#'
#' @template ergmConstraint-general
InitErgmConstraint.edges<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, implies = "edges")
}

#' @templateVar name degrees
#' @title Preserve the degree of each vertex of the given network
#' @description Only networks
#'   whose vertex degrees are the same as those in the network passed
#'   in the model formula have non-zero probability. If the network is
#'   directed, both indegree and outdegree are preserved.
#'
#' @usage
#' # degrees
#'
#' @template ergmConstraint-general
#'
#' @concept directed
#' @concept undirected
InitErgmConstraint.degrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, constrain = "degrees", implies = c("degrees", "edges", "idegrees", "odegrees", "idegreedist", "odegreedist", "degreedist", "bd"))
}

#' @templateVar name degrees
#' @template ergmConstraint-rdname
#' @aliases nodedegrees-ergmConstraint
#' @usage
#' # nodedegrees
InitErgmConstraint.nodedegrees<-InitErgmConstraint.degrees

#' @templateVar name odegrees
#' @title Preserve outdegree for directed networks
#' @description For directed networks, preserve the outdegree of each vertex of the given
#'   network, while allowing indegree to vary
#'
#' @usage
#' # odegrees
#'
#' @template ergmConstraint-general
#'
#' @concept directed
InitErgmConstraint.odegrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("odegrees", "edges", "odegreedist"))
}

#' @templateVar name idegrees
#' @title Preserve indegree for directed networks
#' @description For directed networks, preserve the indegree of each vertex of the given
#'   network, while allowing outdegree to vary
#'
#' @usage
#' # idegrees
#'
#' @template ergmConstraint-general
#'
#' @concept directed
InitErgmConstraint.idegrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("idegrees", "edges", "idegreedist"))
}

#' @templateVar name b1degrees
#' @title Preserve the actor degree for bipartite networks
#' @description For bipartite networks, preserve the degree for the first mode of each vertex of the given
#'   network, while allowing the degree for the second mode to vary.
#'
#' @usage
#' # b1degrees
#'
#' @template ergmConstraint-general
#'
#' @concept bipartite
InitErgmConstraint.b1degrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE)
  list(dependence = TRUE, implies = c("b1degrees", "edges"))
}

#' @templateVar name b2degrees
#' @title Preserve the receiver degree for bipartite networks
#' @description For bipartite networks, preserve the degree for the second mode of each vertex of the given
#'   network, while allowing the degree for the first mode to vary.
#'
#' @usage
#' # b2degrees
#'
#' @template ergmConstraint-general
#'
#' @concept bipartite
InitErgmConstraint.b2degrees<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, bipartite=TRUE)
  list(dependence = TRUE, implies = c("b2degrees", "edges"))
}

#' @templateVar name degreedist
#' @title Preserve the degree distribution of the given network
#' @description Only networks
#'   whose degree distributions are the same as those in the network passed
#'   in the model formula have non-zero probability.
#'
#' @usage
#' # degreedist
#'
#' @template ergmConstraint-general
#'
#' @concept directed
#' @concept undirected
InitErgmConstraint.degreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE, implies = c("degreedist", "edges", "idegreedist", "odegreedist"))
}

#' @templateVar name idegreedist
#' @title Preserve the indegree distribution
#' @description Preserve the indegree distribution of the given network.
#'
#' @usage
#' # idegreedist
#'
#' @template ergmConstraint-general
#'
#' @concept directed
InitErgmConstraint.idegreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("idegreedist", "edges"))
}

#' @templateVar name odegreedist
#' @title Preserve the outdegree distribution
#' @description Preserve the outdegree distribution of the given network.
#'
#' @usage
#' # odegreedist
#'
#' @template ergmConstraint-general
#'
#' @concept directed
InitErgmConstraint.odegreedist<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist, directed=TRUE)
  list(dependence = TRUE, implies = c("odegreedist", "edges"))
}

#' @templateVar name bd
#' @title Constrain maximum and minimum vertex degree
#' @description Condition on the number of inedge or outedges posessed by a node. 
#' See Placing Bounds on Degrees section for more information. ([`?ergmConstraint`][ergmConstraint])
#'
#' @usage
#' # bd(attribs, maxout, maxin, minout, minin)
#' @param attribs a matrix of logicals with dimension `(n_nodes, attrcount)` for the attributes on which we are
#'   conditioning, where `attrcount` is the number of distinct attributes values to condition on.
#' @param maxout,maxin,minout,minin matrices of alter attributes with the same dimension as `attribs` when used
#'   in conjunction with `attribs`. Otherwise, vectors of integers specifying the relevant limits.
#'   If the vector is of length 1, the limit is applied to all nodes. If an individual entry is `NA`,
#'   then there is no restriction of that kind is applied.
#'
#' @template ergmConstraint-general
#'
#' @concept directed
#' @concept undirected
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

#' @templateVar name blocks
#' @title Constrain blocks of dyads defined by mixing type on a vertex attribute.
#' @description Any dyad whose toggle would produce a nonzero change statistic
#'              for a `nodemix` term with the same arguments will be fixed. Note
#'              that the `levels2` argument has a different default value for
#'              `blocks` than it does for `nodemix`.
#'
#' @usage
#' # blocks(attr=NULL, levels=NULL, levels2=FALSE, b1levels=NULL, b2levels=NULL)
#'
#' @template ergmConstraint-general
#' @template ergmTerm-attr
#' @param b1levels,b2levels,levels,level2 control what mixing types are fixed.
#'        `levels2` applies to all networks; `levels` applies to unipartite networks;
#'        `b1levels` and `b2levels` apply to bipartite networks (see Specifying Vertex
#'        attributes and Levels (`?nodal_attributes`) for details)
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmConstraint.blocks <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "b1levels", "b2levels", "levels", "levels2"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC,
                                   ERGM_LEVELS_SPEC, ERGM_LEVELS_SPEC),
                      defaultvalues = list(NULL, NULL, NULL, NULL, FALSE),
                      required = c(TRUE, FALSE, FALSE, FALSE, FALSE))

  if(is.bipartite(nw)) {
    row_nodecov <- ergm_get_vattr(a$attr, nw, bip = "b1")
    col_nodecov <- ergm_get_vattr(a$attr, nw, bip = "b2")

    row_levels <- ergm_attr_levels(a$b1levels, row_nodecov, nw, sort(unique(row_nodecov)))
    col_levels <- ergm_attr_levels(a$b2levels, col_nodecov, nw, sort(unique(col_nodecov)))

    offset <- length(row_levels) + 1L
  } else {
    all_nodecov <- ergm_get_vattr(a$attr, nw)
    row_nodecov <- col_nodecov <- all_nodecov

    all_levels <- ergm_attr_levels(a$levels, all_nodecov, nw, sort(unique(all_nodecov)))
    row_levels <- col_levels <- all_levels

    offset <- 0L
  }

  levels2_list <- transpose(expand.grid(row = row_levels,
                                        col = col_levels, stringsAsFactors = FALSE))
  indices2_grid <- expand.grid(row = seq_along(row_levels),
                               col = offset + seq_along(col_levels))

  if(!is.directed(nw) && !is.bipartite(nw)) {
    rows_leq_cols <- indices2_grid$row <= indices2_grid$col
    levels2_list <- levels2_list[rows_leq_cols]
    indices2_grid <- indices2_grid[rows_leq_cols,]
  }

  levels2_selected <- ergm_attr_levels(a$levels2,
                                       list(row = row_nodecov, col = col_nodecov),
                                       nw,
                                       levels2_list)

  rows_to_keep <- match(levels2_selected, levels2_list, nomatch = NA)
  rows_to_keep <- rows_to_keep[!is.na(rows_to_keep)]

  pairs_to_fix <- indices2_grid[rows_to_keep,]

  if(is.bipartite(nw)) {
    row_nodecov <- match(row_nodecov, row_levels, nomatch = length(row_levels) + 1)
    col_nodecov <- match(col_nodecov, col_levels, nomatch = length(col_levels) + 1)

    nodecov <- c(row_nodecov, col_nodecov + offset)
  } else {
    nodecov <- match(all_nodecov, all_levels, nomatch = length(all_levels) + 1)
  }

  size <- length(col_levels) + 1 + offset
  amat <- matrix(TRUE, nrow = size, ncol = size)
  amat[as.matrix(pairs_to_fix)] <- FALSE

  if(is.bipartite(nw)) {
    amat[,seq_len(offset)] <- FALSE
    amat[-seq_len(offset),] <- FALSE
  } else if(!is.directed(nw)) {
    amat <- amat & t(amat)
  }

  n <- as.integer(network.size(nw))

  rm(nw, arglist, "...") # All needed information has now been extracted.

  free_dyads <- function() {
    rle_list <- lapply(seq_len(NCOL(amat)), function(i) rle(c(amat[nodecov,i])))
    lens <- lapply(seq_len(n), function(i) rle_list[[nodecov[i]]]$lengths)
    vals <- lapply(seq_len(n), function(i) rle_list[[nodecov[i]]]$values)
    rlebdm(compress(structure(list(lengths = unlist(lens),
                                   values = unlist(vals)), class = "rle")), n)
  }

  list(constrain = "blocks",
       dependence = FALSE,
       free_dyads = free_dyads,
       nodecov = nodecov,
       amat = amat)
}


#' @templateVar name hamming
#' @title Preserve the hamming distance to the given network (BROKEN: Do NOT Use)
#' @description This constraint is currently broken. Do not use.
#'
#' @usage
#' # hamming
#'
#' @template ergmConstraint-general
#' @concept directed
#' @concept undirected
InitErgmConstraint.hamming<-function(nw, arglist, ...){
  stop("Constraint ", sQuote("hamming"), " is currently broken and will hopefully be fixed in a future release.")
  a <- check.ErgmTerm(nw, arglist)
  list(dependence = TRUE)
}

#' @templateVar name observed
#' @title Preserve the observed dyads of the given network
#' @description Preserve the observed dyads of the given network.
#'
#' @usage
#' # observed
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
InitErgmConstraint.observed <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist)
  list(free_dyads = as.rlebdm(as.edgelist(is.na(nw))),
       dependence = FALSE, implies = c("observed"))
}

#' @templateVar name fixedas
#' @title Preserve and preclude edges
#' @description Preserve the edges in 'present' and preclude the edges in 'absent'.
#'
#' @usage
#' # fixedas(present, absent)
#' @param present,absent edgelist or network
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
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

#' @templateVar name fixallbut
#' @title Preserve the dyad status in all but the given edges
#' @description Preserve the dyad status in all but `free.dyads`.
#'
#' @usage
#' # fixallbut(free.dyads)
#' @param free.dyads edgelist or network. Networks will be converted to the corresponding edgelist.
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
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
      else free.dyads <- as.edgelist(free.dyads,
                                     n=nw%n%"n",
                                     directed=nw%n%"directed",
                                     bipartite=nw%n%"bipartite",
                                     loops=nw%n%"loops")
      as.rlebdm(free.dyads)
    },
    dependence = FALSE)
}

#' @templateVar name dyadnoise
#' @title A soft constraint to adjust the sampled distribution for
#'   dyad-level noise with known perturbation probabilities
#' @description It is assumed that the observed LHS network is a noisy observation of
#'   some unobserved true network, with `p01` giving the dyadwise
#'   probability of erroneously observing a tie where the true network
#'   had a non-tie and `p10` giving the dyadwise probability of
#'   erroneously observing a nontie where the true network had a tie.
#'   
#' @usage
#' # dyadnoise(p01, p10)
#' @param p01,p10 can both be scalars or both be adjacency matrices of the same dimension as that of the
#'    LHS network giving these probabilities.
#'
#' @template ergmConstraint-general
#'
#' @note See Karwa et al. (2016) for an application.
#' @concept soft
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
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

#' @templateVar name egocentric
#' @title Preserve values of dyads incident on vertices with given attribute
#' @description Preserve values of dyads incident on vertices with attribute `attr` being `TRUE` or if `attrname` is `NULL` , the vertex attribute `"na"` being `FALSE`.
#'
#' @usage
#' # egocentric(attr=NULL, direction="both")
#' @template ergmTerm-attr
#' @param direction one of `"both"`, `"out"` and `"in"`, only applies to directed networks. `"out"` only preserves the out-dyads of those actors and `"in"` preserves their in-dyads.
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
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

#' @templateVar name Dyads
#' @title Constrain fixed or varying dyad-independent terms
#' @description This is an "operator" constraint that takes one or two [`ergmTerm`] dyad-independent formulas. For the terms in the `vary=` formula, only those that change at least one of the terms will be allowed to vary, and all others will be fixed. If both formulas are given, the dyads that vary either for one or for the other will be allowed to vary. Note that a formula passed to `Dyads` without an argument name will default to `fix=` .
#'
#' @usage
#' # Dyads(fix=NULL, vary=NULL)
#' @param fix,vary formula with only dyad-independent terms
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept operator
#' @concept directed
#' @concept undirected
InitErgmConstraint.Dyads<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("fix", "vary"),
                      vartypes = c("formula", "formula"),
                      defaultvalues = list(NULL, NULL),
                      required = c(FALSE, FALSE))
  fix <- a$fix; vary <- a$vary

  if(is.null(fix) & is.null(vary))
    ergm_Init_abort(paste("Dyads constraint takes at least one argument, either",sQuote("fix"),"or",sQuote("vary"),"or both."))

  for(f in c(fix, vary)){
    f[[3]] <- f[[2]]
    f[[2]] <- nw
    if(!is.dyad.independent(f)) ergm_Init_abort(paste("Terms passed to the Dyads constraint must be dyad-independent."))
  }

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
