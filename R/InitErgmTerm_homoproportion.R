#  File ergm/R/InitErgmTerm_homoproportion.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#########################################################
InitErgmTerm.homoproportion<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "keep"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
  u <- sort(unique(nodecov))
  if (!is.null(a$keep)) {
    u <- u[a$keep]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov==(length(u)+1)
  nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
  ui <- seq(along=u)
  ## Extract
  ng <- tabulate(nodecov)
  nn <- sum(ng)
  multfactor <- (nn*(nn-1))/sum(ng*(ng-1))
  ### Construct the list to return
  coef.names <- paste("homoproportion", paste(a$attrname,collapse="."), sep=".")
  inputs <- c(nodecov,multfactor)
  list(name="homoproportion",                            #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs,
       emptynwstats=0,
       dependence = TRUE, # So we don't use MCMC if not necessary
       minval = 0
       )
}
