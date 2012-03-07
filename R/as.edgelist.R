#  File ergm/R/as.edgelist.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
as.edgelist <- function(nw, attrname = NULL, as.sna.edgelist = FALSE,...){
  el <- as.matrix.network.edgelist(nw, attrname=attrname, as.sna.edgelist=as.sna.edgelist,...)
  if(!is.directed(nw)) el[,1:2] <- cbind(pmin(el[,1],el[,2]),pmax(el[,1],el[,2]))
  el
}
