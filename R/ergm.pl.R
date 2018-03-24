#  File R/ergm.pl.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#' @rdname ergm.mple
#' @description \code{ergm.pl} is an even more internal workhorse
#'   function that prepares many of the components needed by
#'   \code{ergm.mple} for the regression rountines that are used to
#'   find the MPLE estimated ergm. It should not be called directly by
#'   the user.
#' @param theta.offset a numeric vector used to specify the offset
#'   (i.e., fixed) coefficients when `ignore.offset=FALSE`.
#' @param ignore.offset If \code{FALSE} (the default), columns
#'   corresponding to terms enclosed in \code{offset()} are not
#'   returned with others but are instead processed by multiplying
#'   them by their corresponding coefficients (which are fixed, by
#'   virtue of being offsets) and the results stored in a separate
#'   column.
#' @return \code{ergm.pl} returns a list containing: \itemize{ \item
#'   xmat : the compressed and possibly sampled matrix of change
#'   statistics \item zy : the corresponding vector of responses,
#'   i.e. tie values \item foffset : ?? \item wend : the vector of
#'   weights for 'xmat' and 'zy' \item numobs : the number of dyads
#'   \item xmat.full: the 'xmat' before sampling; if no sampling is
#'   needed, this is NULL \item zy.full : the 'zy' before sampling; if
#'   no sampling is needed, this is NULL \item foffset.full : ?? \item
#'   theta.offset : a numeric vector whose ith entry tells whether the
#'   the ith curved coefficient?? was offset/fixed; -Inf implies the
#'   coefficient was fixed, 0 otherwise; if the model hasn't any
#'   curved terms, the first entry of this vector is one of
#'   log(Clist$nedges/(Clist$ndyads-Clist$nedges))
#'   log(1/(Clist$ndyads-1)) depending on 'Clist$nedges' \item
#'   maxMPLEsamplesize: the 'maxMPLEsamplesize' inputted to
#'   \code{ergm.pl} }


ergm.pl<-function(Clist, fd, m, theta.offset=NULL,
                    maxMPLEsamplesize=1e+6,
                    control, ignore.offset=FALSE,
                    verbose=FALSE) {
  bip <- Clist$bipartite
  n <- Clist$n

  maxNumDyadTypes <- min(control$MPLE.max.dyad.types,
                         ifelse(bip>0, bip*(n-bip), 
                                ifelse(Clist$dir, n*(n-1), n*(n-1)/2)))
                        
  # May have to think harder about what maxNumDyadTypes should be if we 
  # implement a hash-table approach to compression.
  # *** don't forget, pass in tails first now, not heads
  z <- .C("MPLE_wrapper",
          as.integer(Clist$tails), as.integer(Clist$heads),
          as.integer(Clist$nedges),
          as.double(pack_rlebdm_as_numeric(fd)),
          as.integer(n), 
          as.integer(Clist$dir),     as.integer(bip),
          as.integer(Clist$nterms), 
          as.character(Clist$fnamestring), as.character(Clist$snamestring),
          as.double(Clist$inputs),
          y = integer(maxNumDyadTypes),
          x = double(maxNumDyadTypes*Clist$nstats),
          weightsvector = integer(maxNumDyadTypes),
          as.integer(.Machine$integer.max), # maxDyads
          as.integer(maxNumDyadTypes),
          PACKAGE="ergm")
  uvals <- z$weightsvector!=0
  if (verbose) {
    message(paste("MPLE covariate matrix has", sum(uvals), "rows."))
  }
  zy <- z$y[uvals]
  wend <- as.numeric(z$weightsvector[uvals])
  informative.ties <- sum(wend[zy==1])
  xmat <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
  colnames(xmat) <- m$coef.names
  rm(z,uvals)

  # If we ran out of space, AND we have a sparse network, then, use
  # case-control MPLE.
  if(sum(wend)<sum(fd) && mean(zy)<1/2){
    if(verbose) message("A sparse network with too many unique dyads encountered. Using case-control MPLE.")
    # Strip out the rows associated with ties.
    wend <- wend[zy==0]
    xmat <- xmat[zy==0,,drop=FALSE]
    zy <- zy[zy==0]

    el <- as.edgelist(cbind(Clist$tails, Clist$heads), n, directed=TRUE, bipartite=FALSE, loops=TRUE) # This will be filtered by fd anyway.
    ## Run a whitelist PL over all of the toggleable edges in the network.
    presentrle <- as.rlebdm(el) & fd
    z <- .C("MPLE_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads),
            as.integer(Clist$nedges),
            as.numeric(pack_rlebdm_as_numeric(presentrle)),
            as.integer(n), 
            as.integer(Clist$dir),     as.integer(bip),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring), as.character(Clist$snamestring),
            as.double(Clist$inputs),
            y = integer(maxNumDyadTypes),
            x = double(maxNumDyadTypes*Clist$nstats),
            weightsvector = integer(maxNumDyadTypes),
            as.integer(.Machine$integer.max), # maxDyads
            as.integer(maxNumDyadTypes),
            PACKAGE="ergm")
    uvals <- z$weightsvector!=0
    zy.e <- z$y[uvals]
    wend.e <- as.numeric(z$weightsvector[uvals])
    xmat.e <- matrix(z$x, ncol=Clist$nstats, byrow=TRUE)[uvals,,drop=FALSE]
    colnames(xmat.e) <- m$coef.names
    rm(z,uvals)

    # Divvy up the sampling weight of the ties:
    informative.ties <- wend.e.total <- sum(presentrle)
    wend.e <- wend.e / sum(wend.e) * wend.e.total

    # Divvy up the sampling weight of the nonties:
    wend <- wend / sum(wend) * (sum(fd)-wend.e.total)

    zy <- c(zy,zy.e)
    wend <- c(wend, wend.e)
    xmat <- rbind(xmat, xmat.e)

    rm(zy.e, wend.e, xmat.e)
  }

  #
  # Adjust for the offset
  #

  if(any(m$etamap$offsettheta) && !ignore.offset){
    if(any(is.na(theta.offset[m$etamap$offsettheta]))){
      stop("Offset terms without offset coefficients specified!")
    }
    # Compute the offset's effect.
    foffset <- .multiply.with.inf(xmat[,m$etamap$offsetmap,drop=FALSE], 
                                  cbind(ergm.eta(theta.offset,m$etamap)[m$etamap$offsetmap]))
    
    # Remove offset covariate columns.
    xmat <- xmat[,!m$etamap$offsettheta,drop=FALSE] 
    colnames(xmat) <- m$coef.names[!m$etamap$offsettheta]
    # Now, iff a row's offset effect is infinite, then it carries no
    # further information whatsoever, so it should be dropped.
    xmat <- xmat[is.finite(foffset),,drop=FALSE]
    zy <- zy[is.finite(foffset)]
    wend <- wend[is.finite(foffset)]
    foffset <- foffset[is.finite(foffset)]
  }else{
    foffset <- rep(0, length=length(zy))
  }
  
#
# Sample if necessary
#
  if(nrow(xmat) > maxMPLEsamplesize){
   rsample <- sample((1:nrow(xmat))[zy==1], size=min(maxMPLEsamplesize,sum(zy)),
                     replace=FALSE)
   rsample <- c(rsample, 
     sample((1:nrow(xmat))[zy==0], size=min(maxMPLEsamplesize,sum(!zy)),
                     replace=FALSE) )
   tau <- sum(zy*wend)/sum(wend)
   xmat.full <- xmat
   zy.full <- zy
   foffset.full <- foffset
   zy <- zy[rsample]
   wend <- wend[rsample]
   wend <- tau*zy*wend/sum(zy*wend) +
           (1-tau)*(1-zy)*wend/sum((1-zy)*wend)
   wend <- wend*nrow(xmat)/sum(wend)
   xmat <- xmat[rsample,,drop=FALSE]
   foffset <- foffset[rsample]
  }else{
   xmat.full <- NULL
   zy.full <- NULL
   foffset.full <- NULL
  }

  list(xmat=xmat, zy=zy, foffset=foffset, wend=wend, numobs=round(sum(wend)),
       xmat.full=xmat.full, zy.full=zy.full, foffset.full=foffset.full,
       theta.offset=theta.offset, maxMPLEsamplesize=maxMPLEsamplesize)
}
