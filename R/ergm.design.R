#  File ergm/R/ergm.design.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#################################################################################
# The <ergm.design> function functions as <ergm.Cprepare> would, but acts on the
# network of missing edges
################################################################################

ergm.design <- function(nw, model, verbose=FALSE){
  if(network.naedgecount(nw)==0){
    Clist.miss <- list(tails=NULL, heads=NULL, nedges=0, dir=is.directed(nw))
  }else{
    Clist.miss <- ergm.Cprepare(is.na(nw), model)
    if(verbose){
      cat("Design matrix:\n")
      print(summary(is.na(nw)))
    }
  }
  Clist.miss
}

ergm.Cprepare.miss <- function(nw)
  ergm.Cprepare.el(is.na(nw))

