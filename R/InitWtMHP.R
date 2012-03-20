InitWtMHP.Poisson <- function(arguments, nw, response) {
  MHproposal <- list(name = "Poisson", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePoisson"
  }
  MHproposal
}
#ergm.MHP.table("c", "Poisson", "",  0, "random", "Poisson")


InitWtMHP.ZIPoisson <- function(arguments, nw, response) {
  if(! "p0" %in% names(arguments)){
    arguments$p0 <- max(sum(nw %e% response > 0)/network.dyadcount(nw) - exp(-sum(nw %e% response)/network.dyadcount(nw)),0)
    if(arguments$p0==0) warning("The data do not appear to zero-inflated and are likely to be zero-deflated.")
    cat("Using adaptive jump-to-0 probability of ",arguments$p0,".\n")
  }else if(!is.numeric(arguments$p0) || length(arguments$p0)!=1 || arguments$p0<0 || arguments$p0 >=1) {
    stop(paste("Invalid jump-to-0 argument to the Zero-Inflated Poisson proposal: must be either omited or a single number in the interval [0,1)."))
  }
  MHproposal <- list(name = "ZIPoisson", inputs=as.double(arguments$p0), package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteZIPoisson"
  }
  MHproposal
}
#ergm.MHP.table("c", "Poisson", "",  0, "0inflated", "ZIPoisson")

InitWtMHP.PoissonNonObserved <- function(arguments, nw, response) {
  MHproposal <- list(name = "PoissonNonObserved", inputs=ergm.Cprepare.miss(nw), package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePoissonNonObserved"
  }
  MHproposal
}
#ergm.MHP.table("c", "Poisson", "observed",  0, "random", "PoissonNonObserved")

InitWtMHP.DescRank <- function(arguments, nw, response) {
  MHproposal <- list(name = "CompleteOrdering", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "CompleteOrderingBipartite"
  }
  MHproposal
}
#ergm.MHP.table("c", "DescRank", "",  0, "random", "DescRank")

InitWtMHP.DescRankEquivalent <- function(arguments, nw, response) {
  MHproposal <- list(name = "CompleteOrderingEquivalent", package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "CompleteOrderingEquivalentBipartite"
  }
  
  # Construct the data structure for not-too-inefficient sampling
  n <- if(is.bipartite(nw)) nw %n% "bipartite" else network.size(nw)
  m <- as.matrix(nw, matrix.type = "adjacency", attrname = response)
  diag(m)<-NA

  m<-apply(m, 1,
        function(x){
          cs<-js<-c()
          for(c in sort(unique(x))){ # For each equivalence class
            j<-which(x==c)
            njs<-length(j)
            if(njs<2) next
            cs<-c(cs,njs)
            js<-c(js,j)
          }
          cs<-c(0,cumsum(cs))+1
          cs<-cs+length(cs)
          c(length(cs)-1,cs,js)
        }
        )

  is<-c(0,cumsum(sapply(m,length))[-length(m)])
  is<-is+length(is)

  MHproposal$inputs<-c(is,unlist(m))
  
  MHproposal
}
#ergm.MHP.table("c", "DescRank", "ranks",  0, "random", "DescRankEquivalent")

InitWtMHP.StdNormal <- function(arguments, nw, response) {
  MHproposal <- list(name = "StdNormal", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteStdNormal"
  }
  MHproposal
}
#ergm.MHP.table("c", "StdNormal", "",  0, "random", "StdNormal")

InitWtMHP.StdNormalRank <- function(arguments, nw, response) {
  if(!is.directed(nw) && !is.bipartite(nw)) stop("StdNormRank: The Standard Normal proposal with rank-constraint only works with directed or bipartite networks.")

  if(is.bipartite(nw)){
    stop("Bipartite rank constraint not yet implemented")
    MHproposal <- list(name = "BipartiteStdNormalRank", inputs=c(m), package="ergm")
  }else{
    m<-as.matrix(nw,attrname=response,matrix.type="adjacency")
    diag(m)<-NA
    # This produces a matrix of "classes" of alters for each ego.
    classes<-apply(m,1,function(x) as.numeric(factor(x)))
    diag(classes)<-0

    MHproposal <- list(name = "StdNormalRank", inputs=c(classes), package="ergm")
  }
  MHproposal
}
#ergm.MHP.table("c", "StdNormal", "ranks",  0, "random", "StdNormalRank")
