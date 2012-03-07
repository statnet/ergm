#  File ergm/R/ergm.godfather.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# <ergm.godfather>:  make the network a proposal it can't refuse. 
# Each toggle has a timestamp, and this function forces the network to make
# all of the changes at each unique timestamp value (in increasing order)
# keeping track of the change statistics that result. Thus, the final
# product is a matrix of change statistics in which the number of rows is
# determined by the # of unique timestamps and the number of columns is
# determined by the ERGM terms as usual.
############################################################################

ergm.godfather <- function(formula, timestamps=NULL, toggles=NULL, sim=NULL,
                           start=NULL, end=NULL,
                           accumulate=FALSE,
                           final.network=FALSE,
                           verbose=FALSE,
                           control=control.godfather()) {
  if(is.null(sim)){
    if(is.null(timestamps) | is.null(toggles)){
      stop("Both 'timestamps' and 'toggle' are required arguments if 'sim' ",
           "is not given.")
    }
    if(is.null(start)) start<-min(timestamps)
    if(is.null(end)) end<-max(timestamps)
  }else{
    if(nrow(sim$changed)==0){
      stop("There are no changes (or too many changes) to compute!")
    }else{
      timestamps <- sim$changed[,1]
      toggles <- sim$changed[,2:3]
      start <- attr(sim$changed,"start")
      end <- attr(sim$changed,"end")
    }
  }

  nw <- ergm.getnetwork(formula)
  m <- ergm.getmodel(formula, nw, initialfit=TRUE)
  Clist <- ergm.Cprepare(nw, m)
  m$obs <- summary(m$formula)
  ots <- order(timestamps)
  toggles <- matrix(toggles[ots,],ncol=2)
  timestamps <- timestamps[ots]
  mincol = apply(toggles,1,which.min)
  toggles[mincol==2,] <- toggles[mincol==2,2:1] # make sure col1 < col2
  maxedges <- final.network * max(control$GF.init.maxedges, 5*Clist$nedges)
  obsstat <- summary(formula)  
  z <- .C("godfather_wrapper",
          as.integer(Clist$tails), as.integer(Clist$heads), 
          as.integer(Clist$nedges),
          as.integer(Clist$n),
          as.integer(Clist$dir), as.integer(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring), 
          as.integer(length(timestamps)), as.integer(timestamps),
          as.integer(toggles[,1]), as.integer(toggles[,2]),
          as.integer(start), as.integer(end),
          as.double(Clist$inputs),
          s = double((2+end-start) * Clist$nstats),
          newnwtails = integer(maxedges+1),
          newnwheads = integer(maxedges+1),
          as.integer(accumulate),
          as.integer(verbose),
          as.integer(maxedges), 
          PACKAGE="ergm")  
  stats <- matrix(z$s + obsstat, ncol=Clist$nstats, byrow=TRUE)
  colnames(stats) <- m$coef.names

  out <- list(stats = stats, timestamps = c(NA, start:end))
  if (final.network) { 
    cat("Creating new network...\n")
    out$newnetwork <- newnw.extract(nw,z)
  }
  out
}




####################################################################
# The <control.godfather> function allows for tuning of the
# <ergm.godfather> function
####################################################################

control.godfather<-function(GF.init.maxedges=100000
              ){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))
    control
  }
