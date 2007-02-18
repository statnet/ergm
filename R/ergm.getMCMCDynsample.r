ergm.getMCMCDynsample <- function(g, model, model.dissolve, 
                                  MHproposal, eta0, MCMCparams, 
                                  verbose, BD) {
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of MHproposal, MCMCparams, verbose, and BD should certainly
# be part of Clist.  But this is a project for another day!
#
#   Check for truncation of the returned edge list
#
  Clist <- ergm.Cprepare(g, model)
  Clist.dissolve <- ergm.Cprepare(g, model.dissolve)
  maxedges <- max(5000, Clist$nedges)
  z <- list(newnwhead=maxedges+1)
  while(z$newnwhead[1] > maxedges){
    maxedges <- 5*maxedges
    z <- .C("MCMCDyn_wrapper",
          as.double(Clist$heads), as.double(Clist$tails), 
          as.double(Clist$nedges), as.double(Clist$n),
          as.integer(Clist$dir), as.double(Clist$bipartite),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring),
          as.character(Clist$snamestring),
          as.character(MHproposal$type), as.character(MHproposal$package),
          as.double(Clist$inputs), as.double(eta0),
          as.integer(Clist.dissolve$nterms),
          as.character(Clist.dissolve$fnamestring),
          as.character(Clist.dissolve$snamestring),
          as.double(Clist.dissolve$inputs),
          as.double(MCMCparams$samplesize),
          s = double(MCMCparams$samplesize * Clist$nparam),
          as.double(MCMCparams$burnin), as.double(MCMCparams$interval),
          newnwhead = integer(maxedges), newnwtail = integer(maxedges), 
          diffnwhead = integer(maxedges), diffnwtail = integer(maxedges), 
          numdissolve=integer(1),
          dissolvehead = integer(maxedges), dissolvetail = integer(maxedges), 
          as.integer(verbose), 
          as.double(MCMCparams$gamma), as.integer(MCMCparams$dyninterval), 
          as.integer(BD$attribs), 
          as.integer(BD$maxout), as.integer(BD$maxin),
          as.integer(BD$minout), as.integer(BD$minin),
          as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
          as.double(maxedges),
          as.double(0.0), as.double(0.0), 
          as.double(0.0), as.integer(0),
          PACKAGE="statnet") 
    }
    if(z$newnwhead[1]>1){
#    newedgelist <- matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE)
     newedgelist <- cbind(z$newnwhead[2:z$newnwhead[1]],z$newnwtail[2:z$newnwhead[1]])
    }else{
     newedgelist <- matrix(0, ncol=2, nrow=0)
    }
    if(z$numdissolve>0){
     dissolveedgelist <- cbind(z$dissolvehead[1:z$numdissolve],
                               z$dissolvetail[1:z$numdissolve])
    }else{
     dissolveedgelist <- matrix(0, ncol=2, nrow=0)
    }
#   Next create the network of differences from the origianl one
    if(z$diffnwhead[1]>1){
     diffedgelist <- cbind(z$diffnwhead[2:z$diffnwhead[1]],z$diffnwtail[2:z$diffnwhead[1]])
    }else{
     diffedgelist <- matrix(0, ncol=2, nrow=0)
    }
    statsmatrix <- matrix(z$s, nrow=MCMCparams$samplesize,
                          ncol=Clist$nparam,
                          byrow = TRUE)
    colnames(statsmatrix) <- model$coef.names
#
# recenter statsmatrix by mean statistics if necessary
#
  ms <- MCMCparams$meanstats
  if(!is.null(ms)) {
    if (is.null(names(ms)) && length(ms) == length(model$coef.names))
      names(ms) <- model$coef.names
    obs <- MCMCparams$orig.obs
    obs <- obs[match(colnames(statsmatrix), names(obs))]
    ms  <-  ms[match(names(obs), names(ms))]
    matchcols <- match(names(ms), names(obs))
    if (any(!is.na(matchcols))) {
      ms[!is.na(matchcols)] <- ms[!is.na(matchcols)] - obs[matchcols[!is.na(matchcols)]]
      statsmatrix[,!is.na(matchcols)] <- sweep(as.matrix(
         statsmatrix[,!is.na(matchcols)]), 2, ms[!is.na(matchcols)], "-")
    }
  }
  list(statsmatrix=statsmatrix, newedgelist=newedgelist, meanstats=ms,
       diffedgelist=diffedgelist, dissolveedgelist=dissolveedgelist)
}
