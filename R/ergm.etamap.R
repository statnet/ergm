#  File ergm/R/ergm.etamap.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###########################################################################
# The <ergm.etamap> function takes a model object and creates a mapping
# from the model parameters, theta, to the canonical eta parameters;
# the mapping is carried out by <ergm.eta>
###############################################################################

ergm.etamap <- function(model) {
  etamap <- list(canonical = NULL, offsetmap=NULL, offset=model$offset,
                 offsettheta=NULL, curved=list(), etalength=0)
  from <- 1
  to <- 1
  a <- 1
  if (is.null(model$terms)) {
    return(etamap)
  }
  for (i in 1:length(model$terms)) {
    j <- model$terms[[i]]$inputs[2]
    if(model$offset[i]){
     etamap$offsetmap <- c(etamap$offsetmap, rep(TRUE,j))
    }else{
     etamap$offsetmap <- c(etamap$offsetmap, rep(FALSE,j))
    }
    mti <- model$terms[[i]]
    if (is.null(mti$params)) { # Not a curved parameter
      etamap$canonical <- c(etamap$canonical, to:(to+j-1))
      from <- from+j
      to <- to+j
      if(model$offset[i]){
       etamap$offsettheta <- c(etamap$offsettheta, rep(TRUE,j))
      }else{
       etamap$offsettheta <- c(etamap$offsettheta, rep(FALSE,j))
      }
    } else { # curved parameter
      k <- length(mti$params)
      etamap$canonical <- c(etamap$canonical, rep(0, k))
      etamap$curved[[a]] <- list(from=from:(from+k-1),
                                 to=to:(to+j-1),
                                 map=mti$map, gradient=mti$gradient,
                                 cov=mti$eta.cov)  #Added by CTB 1/28/06
      from <- from+k
      to <- to+j
      a <- a+1
      if(model$offset[i]){
       etamap$offsettheta <- c(etamap$offsettheta, rep(TRUE,k))
      }else{
       etamap$offsettheta <- c(etamap$offsettheta, rep(FALSE,k))
      }
    }
  }
  etamap$etalength <- to-1
  etamap
} 

