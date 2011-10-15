InitWtMHP.Poisson <- function(arguments, nw, model, response) {
  MHproposal <- list(name = "Poisson", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePoisson"
  }
  MHproposal
}

InitWtMHP.PoissonNonObserved <- function(arguments, nw, model, response) {
  MHproposal <- list(name = "PoissonNonObserved", inputs=ergm.Cprepare.miss(nw), package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartitePoissonNonObserved"
  }
  MHproposal
}

InitWtMHP.DescRank <- function(arguments, nw, model, response) {
  MHproposal <- list(name = "CompleteOrdering", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "CompleteOrderingBipartite"
  }
  MHproposal
}

InitWtMHP.DescRankEquivalent <- function(arguments, nw, model, response) {
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


InitWtMHP.StdNormal <- function(arguments, nw, model, response) {
  MHproposal <- list(name = "StdNormal", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteStdNormal"
  }
  MHproposal
}

InitWtMHP.StdNormalRank <- function(arguments, nw, model, response) {
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

