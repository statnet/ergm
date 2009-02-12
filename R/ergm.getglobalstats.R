ergm.getglobalstats <- function(nw, m) {
  Clist <- ergm.Cprepare(nw, m)
  #
  #    Calculate the global statistics
  #
  gs <- .C("network_stats_wrapper",
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
  names(gs) <- m$coef.names
  
  #
  # Adjust to global values
  #
                                                                
  # New method:  Use $emptynwstats added to m$terms by the InitErgm function
  # For example, check the InitErgm.degree function.
  i <- 1
  for (j in 1:length(m$terms)) {
    tmp <- m$term[[j]]
    k <- tmp$inputs[2] # Number of statistics for this model term
    if (!is.null(tmp$emptynwstats)) {
      gs[i:(i+k-1)] <- gs[i:(i+k-1)] + tmp$emptynwstats
    }
    i <- i + k
  }
  
  # Old method:  do adjustments on case-by-case basis 
  # This really needs to be eliminated eventually by inserting the
  # correct $emptynwstats values into the appropriate InitErgm functions

  nnodes <- network.size(nw)
  if(is.directed(nw)){
    dyads <- nnodes*(nnodes-1)
  }else{
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw,"bipartite")
      nb2 <- network.size(nw) - nb1
      dyads <- nb2*nb1
      #   temporary! add these back later
      #   dyads <- (dyads*(dyads-1))/2
    }else{
      dyads <- (nnodes*(nnodes-1))/2
    }
  }
  
  # It looks impossible to handle "duration" using the $emptynwstats, so 
  # this is done separately:
  tase <- grep("duration",names(gs))
  if(length(tase) >0){
    gs[tase] <- -gs[tase]
  }
  
  # SEARCH_ON_THIS_TO_TRACK_DOWN_TRIADCENSUS_CHANGE
  # to undo triadcensus change, uncomment next 8 lines:
#  tase <- match("triadcensus.003",names(gs))
#  if(!is.na(tase)){
#    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6-gs[tase]
#  }
#  tase <- match("triadcensus.0",names(gs))
#  if(!is.na(tase)){
#    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6-gs[tase]
#  }


## Please stop adding ad hoc global-stat-changing code here.  Use
## emptynwstats in the InitErgm function instead.

  tase <- grep("hamming\\.",names(gs))
  if(length(tase) > 0){
    for(i in seq(along=m$terms)){
     if(m$terms[[i]]$name=="hamming"){
       ng0 <-  m$terms[[i]]$inputs[4]
       gs[tase] <- ng0+gs[tase]
     }
    }
  }
  tase <- grep("hammingdyadcov\\.",names(gs))
  if(length(tase) > 0){
    for(i in seq(along=m$terms)){
     if(m$terms[[i]]$name=="hammingdyadcov"){
       ng0 <-  m$terms[[i]]$inputs[4]
       dimnw=dim(as.sociomatrix(nw))
       g0 <- matrix(m$terms[[i]]$inputs[4+(1:(2*ng0))],ncol=2)
       g0 <- cbind(g0[,2],g0[,1]-dimnw[1])
       covm <- array(m$terms[[i]]$inputs[-c(1:(2*ng0+4))], dim=dimnw)
#      ss <- 0
#      for(i in 1:nrow(g0)){
#       ss <- ss + covm[g0[i, 1], g0[i, 2]]
#      }
       g00 <- g0[,1]+dimnw[1]*(g0[,2]-1)
       gs[tase] <- sum(covm[g00])+gs[tase]
     }
    }
  }
  tase <- grep("hammingfixmix\\.",names(gs))
  if(length(tase) > 0){
    for(i in seq(along=m$terms)){
     if(m$terms[[i]]$name=="hammingfixmix"){
       ng0 <-  m$terms[[i]]$inputs[4]
       gs[tase] <- ng0+gs[tase]
     }
    }
  }

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


