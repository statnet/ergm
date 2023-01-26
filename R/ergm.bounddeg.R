#  File R/ergm.bounddeg.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################



#' Initializes the parameters to bound degree during sampling
#' 
#' Not normally called directly by user, \code{ergm.bounddeg} initializes the
#' list of parameters used to bound the degree during the Metropolis Hastings
#' sampling process, and issues warnings if the original network doesn't meet
#' the constraints specified by 'bounddeg'.
#' 
#' In some modeling situations, the degree of certain nodes are constrained to
#' lie in a certain range (rather than their theoretically possible range of 0
#' to n-1).  Such sample space constraints may be incorporated into the ergm
#' modeling process, and if so then the MCMC routine is prevented from visiting
#' network states that violate any of these bounds.
#' 
#' In case there are categories of nodes and degree bounds for each set of
#' categories, such constraints may be incorporated as well.  For instance, if
#' the nodes are girls and boys, and there is a maximum of 5 out-ties to boys
#' and a maximum of 5 out-ties to girls for each node, we would define p to be
#' 2, and the nxp matrix attribs would have TRUE in the first column (say) for
#' exactly those nodes that are boys and TRUE in the second column for only the
#' girls.  The maxout matrix would consist of all 5s in this case, and the
#' other arguments would be left as their default values.
#' 
#' Since the observed network is generally the beginning of the Markov chain,
#' it must satisfy all of the degree constraints itself; thus, this function
#' returns an error message if any bound is violated by the observed network.
#' 
#' @param arguments the `arguments` argument passed to the `InitErgmProposal.*()` function; the sub-sublist `arguments$constraints$bd` should be a list of parameters which may contain the following for a
#' network of size n nodes: \itemize{ \item attribs: an nxp matrix, where entry
#' ij is TRUE if node i has attribute j, and FALSE otherwise; default=an nx1
#' matrix of 1's \item maxout : an nxp matrix, where entry ij is the maximum
#' number of out degrees for node i to nodes with attribute j; default=an nxp
#' matrix of the value (n-1) \item maxin : defined similarly to maxout, but
#' ignored for undirected networks; default=an nxp matrix of the value (n-1)
#' \item minout : defined similarly to maxout; default=an nxp matrix of 0's
#' \item minin : defined similarly to maxout, but ignored for undirected
#' networks; default=an nxp matrix of 0's }
#' @param nw the orginal \code{network} specified to \code{ergm} in 'formula'
#' @return a list of parameters used to bound degree during sampling
#' \item{condAllDegExact}{ always `FALSE`}
#' \item{attribs}{ as defined above}
#' \item{maxout}{ as defined above}
#' \item{maxin}{ as defined above}
#' \item{minout}{ as defined above}
#' \item{minin}{ as defined above}
#' @seealso \code{\link{ergm-proposals}}
#' @keywords internal
#' @export
ergm_bd_init <- function(arguments,nw){
  bounddeg <- arguments$constraints$bd
  nnodes=network.size(nw)
  if(is.null(bounddeg) ||
     all(sapply(bounddeg,function(x){length(x)==1 && is.na(x)}))) {
    return(NULL)
  } else {
    attribs <- bounddeg$attribs
    maxout <- bounddeg$maxout
    maxin <- bounddeg$maxin
    minout <- bounddeg$minout
    minin <- bounddeg$minin
    if (is.null(attribs) || all(attribs==0)){ 
      if(any(!is.null(c(minin,minout,maxout,maxin)))){ 
        attribs <- matrix(1,ncol=1,nrow=nnodes)
        # Get degree for each node
        el <- as.edgelist(nw)
        if (!is.directed(nw)) {
          outdeg <- tabulate(as.vector(el), nbins=nnodes)
        } else {
          outdeg <- tabulate(el[,1], nbins=nnodes)
          indeg <- tabulate(el[,2], nbins=nnodes)
        }
      }else{
        attribs <- 0
      }
    }else{   #Tabulate degrees by attributes
      el <- as.edgelist(nw)
      if(!is.directed(nw)){
        el<-rbind(el[,1:2],el[,2:1])    #Need edges going both directions
        outdeg<-apply(attribs,2,function(z){tabulate(el[z[el[,2]],1], nbins=nnodes)})
      }else{
        outdeg<-apply(attribs,2,function(z){tabulate(el[z[el[,2]],1], nbins=nnodes)})
        indeg<-apply(attribs,2,function(z){tabulate(el[z[el[,1]],2], nbins=nnodes)})
      }
    }

    # Convert degree bounds into matrices if they aren't already, and
    # if they are NULL, make them NA. Also, ensure they are in the
    # right range (though bipartite networks' degrees could be further constrained).
    minin  <- pmax(pmin(matrix(NVL(minin , NA), ncol=ncol(attribs), nrow=nnodes), nnodes-1), 0)
    minout <- pmax(pmin(matrix(NVL(minout, NA), ncol=ncol(attribs), nrow=nnodes), nnodes-1), 0)
    maxin  <- pmax(pmin(matrix(NVL(maxin , NA), ncol=ncol(attribs), nrow=nnodes), nnodes-1), 0)
    maxout <- pmax(pmin(matrix(NVL(maxout, NA), ncol=ncol(attribs), nrow=nnodes), nnodes-1), 0)

    minin [is.na(minin )] <- 0
    minout[is.na(minout)] <- 0
    maxin [is.na(maxin )] <- nnodes-1
    maxout[is.na(maxout)] <- nnodes-1
    if (any(outdeg>maxout | outdeg<minout) || (is.directed(nw) && any(indeg>maxin | indeg<minin))) {
      message("Warning:  Initial network does not satisfy degree constraints. ",
          "Proceeding anyway, but final network may not satisfy constraints.")
    }
    attribs[is.na(attribs)] <- 0
    dependence <- TRUE
    constrains <- "bd"    
  }
  list(condAllDegExact=FALSE,
       attribs=as.integer(attribs),
       maxout=as.integer(maxout),
       maxin=as.integer(maxin),
       minout=as.integer(minout),
       minin=as.integer(minin),
       dependence=dependence,
       constrains=constrains)
}

