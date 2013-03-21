InitWtMHP.DescRank <- function(arguments, nw, response) {
  MHproposal <- list(name = "CompleteOrdering", inputs=NULL)
  MHproposal
}

InitWtMHP.DescRankEquivalent <- function(arguments, nw, response) {
  MHproposal <- list(name = "CompleteOrderingEquivalent")
  
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

InitWtMHP.StdNormal <- function(arguments, nw, response) {
  MHproposal <- list(name = "StdNormal", inputs=NULL)
  MHproposal
}

InitWtMHP.StdNormalRank <- function(arguments, nw, response) {
  if(!is.directed(nw) && !is.bipartite(nw)) stop("StdNormRank: The Standard Normal proposal with rank-constraint only works with directed or bipartite networks.")

  if(is.bipartite(nw)){
    stop("Bipartite rank constraint not yet implemented")
    MHproposal <- list(name = "BipartiteStdNormal", inputs=c(m))
  }else{
    m<-as.matrix(nw,attrname=response,matrix.type="adjacency")
    diag(m)<-NA
    # This produces a matrix of "classes" of alters for each ego.
    classes<-apply(m,1,function(x) as.numeric(factor(x)))
    diag(classes)<-0

    MHproposal <- list(name = "StdNormalRank", inputs=c(classes))
  }
  MHproposal
}

InitWtMHP.DiscUnif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnif", inputs=c(a,b))
  MHproposal
}

InitWtMHP.DiscUnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "DiscUnifNonObserved", inputs=c(a,b,ergm.Cprepare.miss(nw)))
  MHproposal
}

InitWtMHP.Unif <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "Unif", inputs=c(a,b))
  MHproposal
}

InitWtMHP.UnifNonObserved <- function(arguments, nw, response) {
  a <- NVL(arguments$reference$a, -Inf)
  b <- NVL(arguments$reference$b, Inf)
  if(!is.finite(a) || !is.finite(b)) stop('Uniform reference measures that are not bounded are not implemented at this time. Specifiy a and b to be finite.')
  MHproposal <- list(name = "UnifNonObserved", inputs=c(a,b,ergm.Cprepare.miss(nw)))
  MHproposal
}
