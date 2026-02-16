#  File R/InitErgmConstraint.blockdiag.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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
#' @details For bipartite graphs, blocks must be have the same
#'   ordering in both partitions. (They do not, however, require that
#'   all blocks be represented in both partitions, but those that
#'   overlap must have the same order.)
#'
#'   If multiple block-diagonal constraints are given, or if
#'   `attr` is a vector with multiple attribute names, blocks
#'   will be constructed on all attributes matching.
#'
#' @usage
#' # blockdiag(attr, noncontig = "merge")
#' @template ergmTerm-attr
#' @param noncontig character: what to do if the blocks are not contiguous? \describe{
#'
#' \item{`"merge"`}{A placeholder option, to treat them as the same block. It may be implemented in the future, and an explicit `"stop"` option added later.}
#'
#' \item{`"split"`}{Treat them as separate blocks.}
#'
#' }
#'
#' @template ergmConstraint-general
#'
#' @concept dyad-independent
#' @concept directed
#' @concept undirected
#' @import rle
InitErgmConstraint.blockdiag<-function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "noncontig"),
                      vartypes = c(ERGM_VATTR_SPEC, "character"),
                      defaultvalues = list(NULL, "merge"),
                      required = c(TRUE, FALSE))

  a$noncontig <- match.arg(a$noncontig, c("merge", "split"))

  contigmsg <- paste0("Current implementation of the block-diagonal constraint requires that either the blocks be contiguous or be treated as separate; this may change in the future. See ", sQuote("ergmConstraint?blockdiag"), " for more information.")

  list(attr=a$attr, warn_noncontig = a$warn_noncontig,
       free_dyads = {
         n <- network.size(nw)
         storage.mode(n) <- "integer"
         check_noncontig <- a$noncontig == "merge"
         a <- c(ergm_get_vattr(a$attr, nw)) # Strip attributes, which confuse rle().
         if (is.bipartite(nw)) {
           bip <- b1.size(nw)
           ea <- a[seq_len(bip)]
           aa <- a[bip+seq_len(n-bip)]
           if (check_noncontig &&
               (anyDuplicated(rle(ea)$values) || anyDuplicated(rle(aa)$values)))
             ergm_Init_stop(contigmsg)

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
           if (check_noncontig && anyDuplicated(a$values))
             ergm_Init_stop(contigmsg)
           rlebdm(compress(do.call(c,rep(
                                       Map(function(blen,bend){rep(rle(c(FALSE,TRUE,FALSE)), c(bend-blen, blen, n-bend), scale="run")},
                                           a$lengths, cumsum(a$lengths)),
                                       a$lengths)
                                   )), n)
         }
       },
       dependence = FALSE)
}
