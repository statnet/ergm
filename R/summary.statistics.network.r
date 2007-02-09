summary.formula <- function(object, ...){
  current.warn <- options()$warn
  options(warn=0)
  trms<-terms(object)
  if (trms[[1]]!="~")
    stop ("Formula must be of form 'y ~ model'.")
  if(length(trms)<3){return(object)}
  parent <- sys.parent()
  rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  while(inherits(rhs,"try-error") & parent > 1){
    parent <- parent - 1
    rhs <- try(eval(trms[[2]],parent), silent = TRUE)
  }
  options(warn=current.warn)
  UseMethod("summary.statistics",object=rhs)
}


summary.statistics.network <- function(object,...,basis=NULL)
{
  current.warn <- options()$warn
  options(warn=0)
  formula <- object
  trms <- terms(formula)
  if(is.network(basis)){
    nw <- basis
    formula <- as.formula(paste(c("basis",as.character(formula)),collapse=" "))
    trms <- terms(formula)
  }else{
    if(length(trms)>2){
      parent <- sys.parent()
      nw <- try(eval(trms[[2]],parent), silent = TRUE)
      while(inherits(nw,"try-error") & parent > 1){
        parent <- parent - 1
        nw <- try(eval(trms[[2]],parent), silent = TRUE)
      }
      if(class(nw) =="network.series")
        nw <- nw$networks[[1]]
      nw <- as.network(nw)
    }else{
      stop("Must specify a network object")
    }
  }
  nnodes <- network.size(nw)
  if(is.directed(nw)){
    dyads <- nnodes*(nnodes-1)
  }else{
    if(is.bipartite(nw)){
      nactors <- get.network.attribute(nw,"bipartite")
      nevents <- network.size(nw) - nactors
      dyads <- nevents*nactors
      #   temporary! add these back later
      #   dyads <- (dyads*(dyads-1))/2
    }else{
      dyads <- (nnodes*(nnodes-1))/2
    }
  }
  
  m <- ergm.getmodel(formula, nw, drop=FALSE)
  Clist <- ergm.Cprepare(nw, m)
  #
  #    Calculate the global statistics
  #
  gs <- -.C("MCMC_global",
           as.double(Clist$heads), as.double(Clist$tails), 
           as.double(Clist$nedges), as.double(Clist$n),
           as.integer(Clist$dir), as.double(Clist$bipartite), 
           as.integer(Clist$nterms), 
           as.character(Clist$fnamestring), as.character(Clist$snamestring), 
           as.double(Clist$inputs),
           gs = double(Clist$nparam),
           PACKAGE="statnet"
           )$gs
  names(gs) <- m$coef.names
  #
  # Adjust to global values
  #
  tdegree0  <- match( "degree0",names(gs)) 
  if(!is.na(tdegree0)){
    gs[tdegree0] <- gs[tdegree0] + Clist$n
  }
  tidegree0  <- grep( "idegree0",names(gs)) 
  if(any(tidegree0 > 0)){
    gs[tidegree0] <- gs[tidegree0] + Clist$n
  }
  todegree0  <- grep( "odegree0",names(gs)) 
  if(any(todegree0 > 0)){
    gs[todegree0] <- gs[todegree0] + Clist$n
  }
  tdegree0  <- match( "adegree0",names(gs)) 
  if(!is.na(tdegree0)){
    gs[tdegree0] <- gs[tdegree0] + nactors
  }
  tdegree0  <- match( "edegree0",names(gs)) 
  if(!is.na(tdegree0)){
    gs[tdegree0] <- gs[tdegree0] + nevents
  }
  tspartner0 <- match("spartner0",names(gs))
  if(!is.na(tspartner0)){
    gs[tspartner0] <- gs[tspartner0] + Clist$nedges 
  }
  tsesp0 <- match("esp0",names(gs))
  if(!is.na(tsesp0)){
    gs[tsesp0] <- gs[tsesp0] + Clist$nedges 
  }
  #
  tgeodeg <- grep("geodegree",names(gs))
  if(length(tgeodeg) >0){
    gs[tgeodeg] <- Clist$n + gs[tgeodeg]
  }
  #
  tisolates <- grep("isolates",names(gs))
  if(length(tisolates) >0){
    gs[tisolates] <- Clist$n + gs[tisolates]
  }
  ts <- grep("sender[1-9]",names(gs))
  if(length(ts) > 0){
    gs[ts] <- sum(as.sociomatrix(nw)[1,]) + gs[ts]
  }
  ts <- grep("receiver[1-9]",names(gs))
  if(length(ts) > 0){
    gs[ts] <- sum(as.sociomatrix(nw)[,1]) + gs[ts]
  }
  tdsp <- grep("dsp0",names(gs))
  if(length(tdsp) >0){
    if(is.bipartite(nw)){
      gs[tdsp] <- nactors*(nactors-1)/2 + nevents*(nevents-1)/2 + gs[tdsp]
    }else{
      gs[tdsp] <- dyads + gs[tdsp]
    }
  }
  tesa <- grep("esa0",names(gs))
  if(length(tesa) >0){
   if(is.bipartite(nw)){
    gs[tesa] <- nevents*(nevents-1)/2 + gs[tesa]
   }else{
    gs[tesa] <- dyads + gs[tesa]
   }
  }
  twdeg <- grep("ewdegree",names(gs))
  if(length(twdeg) >0){
    gs[twdeg] <- Clist$n + gs[twdeg]
  }
  tase <- grep("duration",names(gs))
  if(length(tase) >0){
    gs[tase] <- -gs[tase]
  }
  # tgeosdeg <- grep("geospartner",names(gs))
  # if(length(tgeosdeg) >0){
  #   gs[tgeosdeg] <- Clist$nedges + gs[tgeosdeg]
  # }
  # twdeg <- grep("gwdegree",names(gs))
  # if(length(twdeg) >0){
  #   gs[twdeg] <- Clist$n + gs[twdeg]
  # }
  # tgwesp <- grep("gwesp",names(gs))
  # if(length(tgwesp) >0){
  #   gs[tgwesp] <- Clist$nedges + gs[tgwesp]
  # }
  tase <- grep("heideriandynamic",names(gs))
  if(length(tase) >0){
    gs[tase] <- summary(nw~asymmetric)+gs[tase]
  }
  tase <- match("transitivity",names(gs))
  if(!is.na(tase)){
    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6+gs[tase]
  }
  tase <- match("transitive",names(gs))
  if(!is.na(tase)){
    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6+gs[tase]
  }
  tase <- grep("intransitivedynamic",names(gs))
  if(length(tase) >0){
    gs[tase] <- summary(nw~intransitive)+gs[tase]
#   gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6-gs[tase]
  }
  tase <- match("triadcensus.003",names(gs))
  if(!is.na(tase)){
    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6-gs[tase]
  }
  tase <- match("triadcensus.0",names(gs))
  if(!is.na(tase)){
    gs[tase] <- nnodes*(nnodes-1)*(nnodes-2)/6-gs[tase]
  }
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
  tase <- grep("hammingmix\\.",names(gs))
  if(length(tase) > 0){
    for(i in seq(along=m$terms)){
     if(m$terms[[i]]$name=="hammingmix"){
       ng0 <-  m$terms[[i]]$inputs[4]
       nu <-  m$terms[[i]]$inputs[1]
       nw %v% "nodecov" <- m$terms[[i]]$inputs[-c(1:(2*ng0+4+2*nu))]
       gs[tase] <- summary(nw ~ mix("nodecov"))+gs[tase]
     }
    }
  }
  options(warn=current.warn)
  gs
}

summary.statistics.ergm <- function(object, ..., basis=NULL)
{
  summary.statistics.network(object$formula, ..., basis=basis)
}
