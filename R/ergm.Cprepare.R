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
#            heads       :  the vector of head nodes
#            tails       :  the vector of tail nodes
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

ergm.Cprepare <- function(nw, m) 
{
  # Build an object called Clist that contains all the necessary
  # ingredients to be passed to the C function.
  n <- network.size(nw)
  dir <- is.directed(nw)
  Clist<-list(n=n, dir=dir)
  bip <- nw$gal$bipartite
  if (is.null(bip)) bip <- 0
  Clist$bipartite <- bip
  Clist$ndyads <- n * (n-1) / (2-dir)
  e<-as.matrix.network(nw,matrix.type="edgelist")
  Clist$maxpossibleedges <- min(max(1e+6, 2*nrow(e)), Clist$ndyads)
  if(length(e)==0){
    Clist$nedges<-0
    Clist$heads<-NULL
    Clist$tails<-NULL
  }else{
    if(!is.matrix(e)){e <- matrix(e, ncol=2)}
    Clist$nedges<-dim(e)[1]
    # Ensure that for undirected networks, head<tail.
    if(dir){
      Clist$heads<-e[,1]
      Clist$tails<-e[,2]
    }else{
      Clist$heads<-pmin(e[,1],e[,2])
      Clist$tails<-pmax(e[,1],e[,2])
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
      Clist$fnamestring <- paste(Clist$fnamestring, mo[[i]]$name)
      Clist$snamestring <- paste(Clist$snamestring, 
                                 ifelse(is.null(mo[[i]]$soname), "ergm",
                                        mo[[i]]$soname))
      Clist$inputs <- c(Clist$inputs, mo[[i]]$inputs)
      Clist$nstats <- Clist$nstats + mo[[i]]$inputs[2]
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


