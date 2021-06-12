#  File R/InitErgmTerm.coincidence.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
InitErgmTerm.coincidence<-function(nw, arglist, ..., version=packageVersion("ergm")) {
  if(version <= as.package_version("3.9.4")){
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("d","active"),
                        vartypes = c("matrix","numeric"),
                        defaultvalues = list(NULL,0),
                        required = c(FALSE,FALSE))
	levels <- if(!is.null(a$d)) I(transpose(data.frame(node1 = a$d[,1], node2 = a$d[,2]))) else NULL
  }else{
    a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                        varnames = c("levels","active"),
                        vartypes = c(ERGM_LEVELS_SPEC,"numeric"),
                        defaultvalues = list(NULL,0),
                        required = c(FALSE,FALSE))
	levels <- a$levels
  }
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- network.size(nw)-nb1

  # coincidence was originally implemented in such a way that the second column varies fastest
  # this behavior is preserved here, although it violates the convention used by other terms
  levels2.grid <- expand.grid(node2 = 1:nb2, node1 = 1:nb2)
  levels2.grid <- levels2.grid[,c(2,1)] # make the second column vary fastest
  levels2.grid <- levels2.grid[levels2.grid$node1 < levels2.grid$node2,]

  if(!is.null(levels)) {
	levels2.list <- transpose(levels2.grid)
    levels2.sel <- ergm_attr_levels(levels, list(node1 = 1:nb2, node2 = 1:nb2), nw, levels2.list)
	
	active <- !is.na(match(levels2.list, levels2.sel))
  } else if(a$active > 0){
   active=summary(nw ~ coincidence(active=0)) > a$active
  }
  if(!is.null(levels) || (a$active > 0)){
   coef.names <- paste("coincidence.",levels2.grid[active,1],".", levels2.grid[active,2], sep="")
   b <- cumsum(active)
   b[active==0] <- 0
  }else{
   b <- 1:nrow(levels2.grid)
   coef.names <- paste("coincidence.",levels2.grid[,1],".", levels2.grid[,2], sep="")
  }
  name <- "coincidence"
  inputs <- c(b)
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
}
