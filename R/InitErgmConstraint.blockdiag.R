#  File R/InitErgmConstraint.blockdiag.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################

## FIXME: There is almost certainly a better way to do this.
.consensus.order <- function(x1, x2){
  o <- intersect(x1, x2)
  if(!all(x1[x1 %in% o] == x2[x2 %in% o])) stop("Current implementation of block-diagonal sampling requires the common blocks of egos and blocks of alters to have the same order. See ", sQuote("ergmConstraint?blockdiag"), "for more information.")
  o1 <- c(0, which(x1 %in% o),length(x1)+1)
  o2 <- c(0, which(x2 %in% o),length(x2)+1)
  n <- length(o1) - 1
  v <- c()

  sr <- function(from,to){from + seq_len(to-from + 1) - 1}

  for(i in seq_len(n)){
    v <- c(v, x1[sr(o1[i]+1,o1[i+1]-1)])
    v <- c(v, x2[sr(o2[i]+1,o2[i+1]-1)])
    v <- na.omit(c(v, x1[o1[i+1]]))
  }
  as.vector(v)
}

.double.rle <- function(a1, a2){
  e1 <- rle(a1)
  e2 <- rle(a2)

  o <- .consensus.order(e1$values, e2$values)

  l1 <- e1$lengths[match(o, e1$values)]
  l1 %[f]% is.na <- 0
  l2 <- e2$lengths[match(o, e2$values)]
  l2 %[f]% is.na <- 0

  list(values=o, lengths1=l1, lengths2=l2)
}

#' @templateVar name blockdiag
#' @title Block-diagonal structure constraint
#' @description Force a block-diagonal structure (and its bipartite analogue) on
#'   the network. Only dyads \eqn{(i,j)} for which
#'   `attr(i)==attr(j)` can have edges.
#'
#'   Note that the current implementation requires that blocks be
#'   contiguous for unipartite graphs, and for bipartite
#'   graphs, they must be contiguous within a partition and must have
#'   the same ordering in both partitions. (They do not, however,
#'   require that all blocks be represented in both partitions, but
#'   those that overlap must have the same order.)
#'
#'   If multiple block-diagonal constraints are given, or if
#'   `attr` is a vector with multiple attribute names, blocks
#'   will be constructed on all attributes matching.
#'
#' @usage
#' # blockdiag(attr)
#' @template ergmTerm-attr
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @import rle
InitErgmConstraint.blockdiag<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  contigmsg <- paste0("Current implementation of the block-diagonal constraint requires that the blocks be contiguous. See ", sQuote("ergmConstraint?blockdiag"), " for more information.")

  list(attr=a$attr,
       free_dyads = {
         n <- network.size(nw)
         storage.mode(n) <- "integer"
         a <- c(ergm_get_vattr(a$attr, nw)) # Strip attributes, which confuse rle().
         if(NVL(nw%n%"bipartite",0)){
           bip <- nw %n% "bipartite"
           ea <- a[seq_len(bip)]
           aa <- a[bip+seq_len(n-bip)]
           if(anyDuplicated(rle(ea)$values) || anyDuplicated(rle(aa)$values)) stop(contigmsg)

           tmp <- .double.rle(ea, aa)
           el <- tmp$lengths1
           al <- tmp$lengths2

           o <- rlebdm(c(rep(rle(FALSE), bip*n, scale="run"),
                         do.call(c,rep(
                                     Map(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                         el, cumsum(el)),
                                     al)
                                 )), n)
           # Future-proofing: in case it's bipartite directed, add
           # both thte blocks and their transposes. (If undirected,
           # it'll get filtered out by the .attributes constraints.)
           ot <- rlebdm(c(do.call(c,rep(
                                      Map(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bip+bend-blen, blen, n-bip-bend), scale="run")},
                                          al, cumsum(al)),
                                      el)
                                  ),
                          rep(rle(FALSE), (n-bip)*n, scale="run")), n)
           compress(o | ot)
         }else{
           a <- rle(a)
           if(anyDuplicated(a$values)) stop(contigmsg)
           rlebdm(compress(do.call(c,rep(
                                       Map(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                           a$lengths, cumsum(a$lengths)),
                                       a$lengths)
                                   )), n)
         }
       },
       dependence = FALSE)
}
