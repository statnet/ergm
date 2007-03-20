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
  maxchanges <- max(MCMCparams$maxchanges, Clist$nedges)/5
  MCMCparams$maxchanges <- MCMCparams$maxchanges/5
  z <- list(newnwhead=maxchanges+1)
  while(z$newnwhead[1]  >= maxchanges || 
        z$dissnwhead[1] >= maxchanges ||
        z$diffnwhead[1] >= maxchanges){
    maxchanges <- 5*maxchanges
    MCMCparams$maxchanges <- 5*MCMCparams$maxchanges
    if(verbose){cat(paste("MCMCDyn workspace is",maxchanges,"\n"))}
    z <- .C("MCMCDyn_wrapper",
          as.integer(Clist.dissolve$order.code),
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
          newnwhead = integer(maxchanges), newnwtail = integer(maxchanges), 
          diffnwtime = integer(maxchanges), diffnwhead = integer(maxchanges), diffnwtail = integer(maxchanges), 
          dissnwtime = integer(maxchanges), dissnwhead = integer(maxchanges), dissnwtail = integer(maxchanges),
          as.integer(verbose), 
          as.double(MCMCparams$gamma), as.integer(MCMCparams$dyninterval), 
          as.integer(BD$attribs), 
          as.integer(BD$maxout), as.integer(BD$maxin),
          as.integer(BD$minout), as.integer(BD$minin),
          as.integer(BD$condAllDegExact), as.integer(length(BD$attribs)), 
          as.double(maxchanges),
          as.double(0.0), as.double(0.0), 
          as.double(0.0), as.integer(0),
          PACKAGE="statnet") 
    }
#   cat(paste("z$diffnwhead = ",maxchanges,z$diffnwhead[1],"\n"))
#   cat(paste("z$dissnwhead = ",maxchanges,z$dissnwhead[1],"\n"))
#   cat(paste("z$newwhead = ",maxchanges,z$newnwhead[1],"\n"))
    if(z$newnwhead[1]>1){
#    newedgelist <- matrix(z$newnw[2:z$newnw[1]], ncol=2, byrow=TRUE)
     newedgelist <- cbind(z$newnwhead[2:z$newnwhead[1]],z$newnwtail[2:z$newnwhead[1]])
    }else{
     newedgelist <- matrix(0, ncol=2, nrow=0)
    }
    if(z$dissnwhead[1]>1){
     dissedgelist <- cbind(z$dissnwtime[2:z$dissnwtime[1]],z$dissnwhead[2:z$dissnwhead[1]],z$dissnwtail[2:z$dissnwhead[1]])
    }else{
     dissedgelist <- matrix(0, ncol=3, nrow=0)
    }
#   Next create the network of differences from the origianl one
    if(z$diffnwhead[1]>1){
     diffedgelist <- cbind(z$diffnwtime[2:z$diffnwtime[1]],z$diffnwhead[2:z$diffnwhead[1]],z$diffnwtail[2:z$diffnwhead[1]])
    }else{
     diffedgelist <- matrix(0, ncol=3, nrow=0)
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
       changed=diffedgelist, dissolved=dissedgelist,
       maxchanges=MCMCparams$maxchanges)
}
