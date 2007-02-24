simulatedyn <- function(object, nsim=1, seed=NULL, ...,theta0,
                        burnin=1, interval=1,
                        sequential=TRUE,
                        proposaltype="formationTNT",
                        dissolve=NULL, gamma=0.01,
                        algorithm.control=list(),
                        drop=FALSE,
                        verbose=FALSE) {
  out.list <- list()
  out.mat <- numeric(0)
  formula <- object

  ## Defaults :
  con <- list(boundDeg=NULL, drop=drop,
              proposalpackage="statnet",
              dyninterval=1000,
              maxchanges=20000,
              summarizestats=FALSE
             )

  con[(namc <- names(algorithm.control))] <- algorithm.control
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
  nw <- ergm.getnetwork(formula)
  if(class(nw) =="network.series"){
    nw <- nw$networks[[1]]
  }
  nw <- as.network(nw)
  if(!is.network(nw)){
    stop("A network object on the LHS of the formula must be given")
  }
  # Resolve conditioning
  # This is laborious to cover the various partial specifications
  #

######################################
# Almost all of the following is commented out for now and should be modified.
#BD <- con$boundDeg
#  if(is.null(BD)){
#   BD <- ergm.boundDeg(NULL)
    BD <- ergm.boundDeg(con$boundDeg, nnodes=network.size(nw))
#  }
#
#  if (BD$condAllDegExact==TRUE && proposaltype != "conddeg") {
##     cat("Warning:  If condAllDegExact is set to TRUE inside boundDeg,")
##     cat("then switching must be chosen.  Setting proposaltype == "conddeg" now.\n")
#    proposaltype <- "conddeg"
#  }
#
#  if(mixed||conditional){
#    proposalnumber <- c(2, 3, 7, 8)[match(proposaltype,
#                                          c("conddegdist", "conddeg",
#                                            "condoutdeg", "condindeg"),
#                                          nomatch=1)]
#  }else{
#    proposalnumber <- match(proposaltype,
#                            c("toggle","conddegdist", "conddegdistswitch",
#                              "conddeg", "nodeedges", "node",
#                              "condoutdeg", "condindeg", "constantedges",
#                              "tnt"),
#                            nomatch=1)
#  }
#  
################################################
  m <- ergm.getmodel(formula, nw, drop=con$drop)
  distanceMetric <- 0
  Clist <- ergm.Cprepare(nw, m)
  
  model.dissolve <- ergm.getmodel.dissolve(dissolve, nw)
  Clist.dissolve <- ergm.Cprepare(nw, model.dissolve)
#
  MCMCsamplesize <- 1
  verb <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(is.bipartite(nw)){
   if(proposaltype=="randomtoggle"){proposaltype <- "Bipartiterandomtoggle"}
   if(proposaltype=="ConstantEdges"){proposaltype <- "BipartiteConstantEdges"}
   if(proposaltype=="TNT"){proposaltype <- "BipartiteTNT"}
   if(proposaltype=="HammingConstantEdges"){proposaltype <- "BipartiteHammingConstantEdges"}
   if(proposaltype=="Hamming"){proposaltype <- "BipartiteHamming"}
   if(proposaltype=="formation"){proposaltype <- "BipartiteFormation"}
   if(proposaltype=="formationTNT"){proposaltype <- "BipartiteFormationTNT"}
  }
  if(missing(theta0)) {
    theta0 <- rep(0,Clist$nparam)
    warning("No parameter values given, using Bernouli network\n\t")
  }
  
  if(is.null(seed)){seed <- sample(10000000, size=1)}
  set.seed(as.integer(seed))
    
    maxchanges <- max(con$maxchanges, Clist$nedges)
    use.burnin <- nsim
    z <- list(newnwhead=maxchanges+1)
    while(z$newnwhead[1] > maxchanges){
     maxchanges <- 5*maxchanges
     z <- .C("MCMCDyn_wrapper",
             as.double(Clist$heads), as.double(Clist$tails), 
             as.double(Clist$nedges), as.double(Clist$n),
             as.integer(Clist$dir), as.double(Clist$bipartite),
             as.integer(Clist$nterms), 
             as.character(Clist$fnamestring),
             as.character(Clist$snamestring), 
             as.character(proposaltype),
             as.character(con$proposalpackage),
             as.double(Clist$inputs),
             as.double(theta0),
             as.integer(Clist.dissolve$nterms),
             as.character(Clist.dissolve$fnamestring),
             as.character(Clist.dissolve$snamestring),
             as.double(Clist.dissolve$inputs),
             as.double(MCMCsamplesize),
             s = double(MCMCsamplesize * Clist$nparam),
             as.double(use.burnin), as.double(interval), 
             newnwhead = integer(maxchanges), newnwtail = integer(maxchanges),
             diffnwtime = integer(maxchanges), diffnwhead = integer(maxchanges), diffnwtail = integer(maxchanges),
             dissnwtime = integer(maxchanges), dissnwhead = integer(maxchanges), dissnwtail = integer(maxchanges),
             as.integer(verb),
             as.double(gamma), as.integer(con$dyninterval),
             as.integer(BD$attribs), 
             as.integer(BD$maxout), as.integer(BD$maxin), as.integer(BD$minout),
             as.integer(BD$minin), as.integer(BD$condAllDegExact),
             as.integer(length(BD$attribs)), 
             as.double(maxchanges), 
             as.double(0.0), as.double(0.0), 
             as.double(0.0), as.integer(0),
             PACKAGE="statnet")
    }
#
#   Next update the network to be the final (possibly conditionally)
#   simulated one
#
    if(z$newnwhead[1]>1){
     newedgelist <- cbind(z$newnwhead[2:z$newnwhead[1]],z$newnwtail[2:z$newnwhead[1]])
    }else{
     newedgelist <- matrix(0, ncol=2, nrow=0)
    }
#   Next create the network of differences from the origianl one
    if(z$diffnwhead[1]>1){
     diffedgelist <- cbind(z$diffnwtime[2:z$diffnwtime[1]],z$diffnwhead[2:z$diffnwhead[1]],z$diffnwtail[2:z$diffnwhead[1]])
    }else{
     diffedgelist <- matrix(0, ncol=3, nrow=0)
    }
    if(z$dissnwhead[1]>1){
     dissedgelist <- cbind(z$dissnwtime[2:z$dissnwtime[1]],z$dissnwhead[2:z$dissnwhead[1]],z$dissnwtail[2:z$dissnwhead[1]])
    }else{
     dissedgelist <- matrix(0, ncol=3, nrow=0)
    }
    out.list <- diffedgelist
    out.mat <- rbind(out.mat,z$s[(Clist$nparam+1):(2*Clist$nparam)])

    out.list <- list(formula = formula, networks = NULL,
                       changed=diffedgelist, 
                       dissolved=dissedgelist, 
                       stats = out.mat, coef=theta0)
    class(out.list) <- "network.series"

  return(out.list)
}
