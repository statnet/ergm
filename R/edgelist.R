#  File R/edgelist.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2015 Statnet Commons
#######################################################################
as.edgelist <- function(x, ...) UseMethod("as.edgelist")

as.edgelist.network <- function(x, attrname = NULL, as.sna.edgelist = FALSE, inverted = NULL, ...){
  as.edgelist(as.matrix.network.edgelist(x, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...), n=network.size(x), directed=is.directed(x), bipartite=if(is.bipartite(x)) x%n%"bipartite" else FALSE, loops=has.loops(x), inverted=NVL(inverted, NVL(x%n%"inverted", FALSE)))
}

## .copy.el.attr <- function(to, from){
##   attr(to,"n") <- attr(from,"n")
##   attr(to,"directed") <- attr(from,"directed")
##   attr(to,"bipartite") <- attr(from,"bipartite")
##   attr(to,"loops") <- attr(from,"loops")
##   to
## }

## .check.el.attr <- function(x, y){
##   stopifnot(ncol(x) == ncol(y),
##             attr(x,"n") == attr(y,"n"),
##             attr(x,"directed") == attr(y,"directed"),
##             attr(x,"bipartite") == attr(y,"bipartite"),
##             attr(x,"loops") == attr(y,"loops"))
## }

## is.inverted <- function(x, ...) UseMethod("is.inverted")
## is.inverted.edgelist <- function(x) attr(x, "inverted")

## dyad.count <- function(x, ...) UseMethod("dyad.count")
## dyad.count.network <- function(x, na.omit=FALSE) network.dyadcount(x, na.omit)
## dyad.count.edgelist <- function(x, ...)
##   with(attributes(x),
##        if(bipartite) (n-m)*m*(if(directed) 2 else 1)
##        else n*(n-1)/(if(directed) 1 else 2) + (if(loops) n else 0)
##        )

## edge.count <- function(x, ...) UseMethod("edge.count")
## edge.count.network <- function(x, ...) network.edgecount(x)
## edge.count.edgelist <- function(x, ...)
##   with(attributes(x),
##        if(inverted) dyad.count(x) - nrow(x) else nrow(x)
##        )

as.edgelist.matrix <- function(x, n, directed=TRUE, bipartite=FALSE, loops=FALSE, inverted=FALSE, ...){
  if(!directed) x[,1:2] <- cbind(pmin(x[,1],x[,2]),pmax(x[,1],x[,2]))
  if(!loops) x <- x[x[,1]!=x[,2],,drop=FALSE]
  if(bipartite) x <- x[(x[,1]<=bipartite)!=(x[,2]<=bipartite),,drop=FALSE]
  x <- unique(x[order(x[,1],x[,2]),,drop=FALSE])
  attr(x,"n") <- n
  attr(x,"directed") <- directed
  attr(x,"bipartite") <- bipartite
  attr(x,"loops") <- loops
  attr(x,"inverted") <- inverted
  x
}


## union <- function(x, y, ...) UseMethod("union")
## union.defualt <- function(x, y) base::union(x, y)
## union.edgelist <- function(x, y, ...){
##   .check.el.attr(x,y)

##   if(edge.count(x) + edge.count(y) <= dyad.count(x)/2){
##     if(is.inverted(x)){
##       if(is.inverted(y)){
        
##       }else{
        
##       }
##     }else{
##       if(is.inverted(y)){
        
##       }else{
##         out <- unique(rbind(as.edgelist(x),as.edgelist(y)))
##       }
##     }
##     attr(out, "inverted") <- FALSE
##   }else{ # Use DeMorgan's Rule: (x
##   }
  
##   .copy.el.attr(out,x)
## }


## intersect <- function(x, y, ...) UseMethod("intersect")
## intersect.default <- function(x, y) base::intersect(x, y)
## intersect.edgelist <- function(x, y, ...){
##   .check.el.attr(x,y)

##   if(edge.count(x) + edge.count(y) <= dyad.count(x)/2){
##     if(is.inverted(x)){
##       if(is.inverted(y)){
        
##       }else{
        
##       }
##     }else{
##       if(is.inverted(y)){
        
##       }else{
##         out <- unique(rbind(as.edgelist(x),as.edgelist(y)))
##       }
##     }
##     attr(out, "inverted") <- FALSE
##   }else{ # Use DeMorgan's Rule: (x
##   }
  
##   .copy.el.attr(out,x)
## }
