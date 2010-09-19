ergm.getglobalstats <- function(nw, m, response=NULL) {
  Clist <- ergm.Cprepare(nw, m, response=response)
  #
  #    Calculate the global statistics
  #
  
  gs <- (
         if(is.null(response))
         .C("network_stats_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), 
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = double(Clist$nstats),
            PACKAGE="ergm"
            )$gs
         else
         .C("wt_network_stats_wrapper",
            as.integer(Clist$heads), as.integer(Clist$tails), as.double(Clist$weights), 
            as.integer(Clist$nedges), as.double(Clist$baseline_weight),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite), 
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring), 
            as.double(Clist$inputs),
            gs = double(Clist$nstats),
            PACKAGE="ergm"
            )$gs
         )
  names(gs) <- m$coef.names
  
  #
  # Adjust to global values
  #
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgm function
  # Read the comments at the top of InitErgm.R or InitErgmTerm.R for 
  # an explanation of the $emptynwstats mechanism
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$term[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }
  # Note:  "duration" is not in the CRAN version.  
  # It looks impossible to handle "duration" using the $emptynwstats, so 
  # this is done separately:
  tase <- grep("duration",names(gs))
  if(length(tase) >0){
    gs[tase] <- -gs[tase]
  }

## Please don't add ad hoc global-stat-changing code here.  Use
## emptynwstats in the InitErgm function instead.



# Old hammingmix code.  Can be deleted when no longer needed to
# show how this used to work.
#  tase <- grep("hammingmix\\.",names(gs))
#  if(length(tase) > 0){
#    for(i in seq(along=m$terms)){
#     if(m$terms[[i]]$name=="hammingmix"){
#       ng0 <-  m$terms[[i]]$inputs[4]
#       nu <-  m$terms[[i]]$inputs[1]
#       nw %v% "nodecov" <- m$terms[[i]]$inputs[-c(1:(2*ng0+4+2*nu))]
#       gs[tase] <- summary(nw ~ nodemix("nodecov"))+gs[tase]
#     }
#    }
#  }

  gs
}


