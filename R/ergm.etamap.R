#  File ergm/R/ergm.etamap.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
ergm.etamap <- function(model) {
# Take a model object and create a 'mapping' from the parameter of
# interest, theta, to the canonical parameter, eta.  This 'mapping'
# is a list with three items:
#  The first item, canonical, summarizes the
# one-to-one correspondence between those components of theta that
# are actually canonical and their counterparts in eta.
#  The second item,
# curved, is itself a list, with one element per curved EF term in the
# model.
#  The third item, etalength, is the length of eta.
# The function ergm.eta actually carries out the mapping.
  etamap <- list(canonical = NULL, offsetmap=NULL, offset=model$offset,
                 offsettheta=NULL, curved=list(), etalength=NULL)
  from <- 1
  to <- 1
  a <- 1
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

