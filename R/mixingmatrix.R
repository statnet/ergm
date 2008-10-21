#  File ergm/R/mixingmatrix.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
mixingmatrix <- function(nw, attrname) {
  if(!is.network(nw)){
    stop("mixingmatrix() requires a network object")
  }                                                                
  nodecov <- get.node.attr(nw, attrname, functionname="mixingmatrix")
  u<-sort(unique(nodecov))
  # nodecovnum <- match(nodecov, u)
  el <- as.matrix.network.edgelist(nw)
  type <- "directed"
  if (is.bipartite(nw)) { # must have heads < tails now
    if (is.directed(nw)) 
      cat("Warning:  Bipartite networks are currently\n",
          "automatically treated as undirected\n")
    type <- "bipartite"
    rowswitch <- apply(el, 1, function(x) x[1]>x[2])
    el[rowswitch, 1:2] <- el[rowswitch, 2:1]
    nb1 <- get.network.attribute(nw,"bipartite")
    u<-sort(unique(nodecov[1:nb1]))
    From <- c(u, nodecov[el[,1]])
    u<-sort(unique(nodecov[(nb1+1):network.size(nw)]))
    To <- c(u, nodecov[el[,2]])
  }else{
    From <- c(u, nodecov[el[,1]])
    To <- c(u, nodecov[el[,2]])
  }
  tabu <- table(From, To)  # Add u,u diagonal to ensure each 
  # value is represented, then subtract it later
  diag(tabu) <- diag(tabu) - 1
  if(!is.directed(nw) && !is.bipartite(nw)){
    type <- "undirected"
    tabu <- tabu + t(tabu)
    diag(tabu) <- diag(tabu)/2
  }
  ans <- list(type=type, matrix=tabu)
  class(ans) <- "mixingmatrix"
  ans
}

print.mixingmatrix <- function(x, ...) {
  m <- x$mat
  rn <- rownames(m)
  cn <- colnames(m)  
  if (x$type == "undirected") {
    dimnames(m) <- list(rn, cn)
    cat("Note:  Marginal totals can be misleading\n",
        "for undirected mixing matrices.\n")
  } else {
    total <- apply(m,1,sum)
    m <- cbind(m,total)
    total <- apply(m,2,sum)
    m <- rbind(m,total)
    rn <- c(rn, "Total")
    cn <- c(cn, "Total")
    if (x$type == "bipartite")
      dimnames(m) <- list(B1 = rn,B2 = cn)
    else
      dimnames(m) <- list(From = rn,To = cn)
  }
  print(m)
}



