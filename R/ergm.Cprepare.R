##########################################################################
# The <ergm.Cprepare> function builds an object called Clist that contains
# all the necessary ingredients to be passed to the C functions
#
# --PARAMETERS--
#   nw:  a network object
#   m :  a model object, as returned by <ergm.getmodel>
#
# --RETURNED--
#   Clist:  a list of parameters used by several of the fitting routines
#           containing
#            n           :  the size of the network
#            dir         :  whether the network is directed (T or F)
#            bipartite   :  whether the network is bipartite (T or F)
#            ndyads      :  the number of dyads in the network
#            nedges      :  the number of edges in this network
#            tails       :  the vector of tail nodes; tail nodes are
#                               the 1st column of the implicit edgelist,
#                               so either the lower-numbered nodes in an
#                               undirected graph, or the out nodes of a
#                               directed graph, or the b1 nodes of a bi-
#                               partite graph
#            heads           :  the vector of head nodes; head nodes are
#                               the 2nd column of the implicit edgelist,
#                               so either the higher-numbered nodes in an
#                               undirected graph, or the in nodes of a
#                               directed graph, or the b2 nodes of a bi-
#                               partite graph
#            nterms      :  the number of model terms
#            nstats      :  the total number of change statistics
#                           for all model terms
#            inputs      :  the concatenated vector of 'input's from each
#                           model term as returned by <InitErgmTerm.X> or
#                           <InitErgm.X>
#            fnamestring :  the concatenated string of model term names
#            snamestring :  the concatenated string of package names that
#                           contain the C function 'd_fname'; default="ergm"
#                           for each fname in fnamestring
#            maxpossibleedges :  the maximum number of edges to allocate
#                                space for
#            stergm.order.code:  a numeric code indicating which dissolution
#                                and formation process is to be used, where
#                                    1 = DissThenForm
#                                    2 = DissAndForm or FormAndDiss
#                                    3 = FormThenDiss
#                                    4 = FormOnly
#                                    5 = DissOnly
#                                the default is 0, which will lead to an
#                                error in the C code
#
##########################################################################

ergm.Cprepare <- function(nw, m, response=NULL)
{
  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw$gal$bipartite
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- n * (n-1) / (2-dir)
  e<-as.matrix.network(nw,matrix.type="edgelist",attrname=response)
  Clist$maxpossibleedges <- min(max(1e+6, 2*nrow(e)), Clist$ndyads)
  if(length(e)==0){
    Clist$nedges<-0
    Clist$tails<-NULL
    Clist$heads<-NULL
    ## Make sure weights is not NULL if response!=NULL, even if it's
    ## empty, since it's used to decide whether MCMC or WtMCMC is
    ## called.
    if(!is.null(response)) Clist$weights<-numeric(0)
  }else{
    if(!is.matrix(e)){e <- matrix(e, ncol=2+!is.null(response))}

    ## Delete 0 edges.
    if(!is.null(response)) e<-e[e[,3]!=0,,drop=FALSE]
    
    Clist$nedges<-dim(e)[1]
    # *** Ensure that for undirected networks, tail<head.
    if(dir){
      Clist$tails<-e[,1]
      Clist$heads<-e[,2]
    }else{
      Clist$tails<-pmin(e[,1],e[,2])
      Clist$heads<-pmax(e[,1],e[,2])
    }
    if(!is.null(response)){
      Clist$weights<-e[,3]
    }
  }
  mo<-m$terms 
  
  Clist$nterms<-length(mo)
  Clist$nstats<-0
  Clist$fnamestring<-""
  Clist$snamestring<-""
  Clist$inputs<-numeric(0)
  if (Clist$nterms>0) {
    for(i in 1:Clist$nterms) {
      term_i <- mo[[i]]
      Clist$fnamestring <- paste(Clist$fnamestring, term_i$name)
      # This lets "pkgname" play the same role as "soname":
      Clist$snamestring <- paste(Clist$snamestring, 
                                 if (!is.null(term_i$soname)) {
                                   term_i$soname
                                 } else if (!is.null(term_i$pkgname)) {
                                   term_i$pkgname
                                 } else {
                                   "ergm"
                                 } )
      Clist$inputs <- c(Clist$inputs, term_i$inputs)
      Clist$nstats <- Clist$nstats + term_i$inputs[2]
    }
  }
  while (substring(Clist$fnamestring, 1, 1)==" ")
    Clist$fnamestring <- substring(Clist$fnamestring, 2)
  while (substring(Clist$snamestring, 1, 1)==" ")
    Clist$snamestring <- substring(Clist$snamestring, 2)

  if("stergm.order" %in% names(m)) Clist$stergm.order.code <- switch(m$stergm.order,
                                                              DissThenForm=1,
                                                              DissAndForm=2,
                                                              FormAndDiss=2,
                                                              FormThenDiss=3,
                                                              FormOnly=4,
                                                              DissOnly=5,
                                                              0)
  
  Clist
}


