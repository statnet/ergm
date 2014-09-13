#  File R/InitErgmTerm.coincidence.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
InitErgmTerm.coincidence<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=TRUE,
                      varnames = c("d","active"),
                      vartypes = c("numeric","numeric"),
                      defaultvalues = list(NULL,0),
                      required = c(FALSE,FALSE))
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- network.size(nw)-nb1
  base <- 1:nb2 
  base <- cbind(rep(1:nb2, rep(nb2, nb2)),
                rep(1:nb2, nb2))
  base <- base[base[, 2] > base[, 1], ]
  if(!is.null(a$d)){
   active <- !is.na(match(base [,1] + nb2*base[,2], a$d[,1]+nb2*a$d[,2]))
  }else{if(a$active > 0){
   active=summary(nw ~ coincidence(active=0)) > a$active
  }}
  if(!is.null(a$d) | (a$active > 0)){
   coef.names <- paste("coincidence.",base[active,1],".", base[active,2], sep="")
   b <- cumsum(active)
   b[active==0] <- 0
  }else{
   b <- 1:nrow(base)
   coef.names <- paste("coincidence.",base[,1],".", base[,2], sep="")
  }
  name <- "coincidence"
  inputs <- c(b)
  list(name=name, coef.names=coef.names, inputs=inputs, dependence=TRUE)
}
