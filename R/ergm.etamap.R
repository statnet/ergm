#  File R/ergm.etamap.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2019 Statnet Commons
#######################################################################
###########################################################################
# The <ergm.etamap> function takes a model object and creates a mapping
# from the model parameters, theta, to the canonical eta parameters;
# the mapping is carried out by <ergm.eta>
#
# --PARAMETERS--
#   model: a model object, as returned by <ergm_model>
#
# --RETURNED--
#   etamap: the theta -> eta mapping given by a list of the following:
#     canonical  : a numeric vector whose ith entry specifies whether
#                  the ith component of theta is canonical (via non-
#                  negative integers) or curved (via zeroes)
#     offsetmap  : a logical vector whose ith entry tells whether the
#                  ith coefficient of the canonical parameterization
#                  was "offset", i.e fixed 
#     offset     : a logical vector whose ith entry tells whether the
#                  ith model term was offset/fixed
#     offsettheta: a logical vector whose ith entry tells whether the
#                  ith curved theta coeffient was offset/fixed;
#     curved     : a list with one component per curved EF term in the
#                  model containing
#         from    : the indices of the curved theta parameter that are
#                   to be mapped from
#         to      : the indices of the canonical eta parameters to be
#                   mapped to
#         map     : the map provided by <InitErgmTerm>
#         gradient: the gradient function provided by <InitErgmTerm> 
#         cov     : the eta covariance ??, possibly always NULL (no
#                   <Init> function creates such an item)
#     etalength  : the length of the eta vector
#
###############################################################################

# Note: Documentation for the return data structure is in `ergm.eta.R`.
ergm.etamap <- function(model) {
  etamap <- list(canonical = NULL, offsetmap=NULL, offset=model$offset,
                 offsettheta=NULL, curved=list(), etalength=0L)
  from <- 1L
  to <- 1L
  a <- 1L
  if (is.null(model$terms)) {
    return(etamap)
  }
  for (i in seq_along(model$terms)) {
    mti <- model$terms[[i]]
    if(NVL(mti$offset,FALSE)) model$offset[i] <- TRUE
    j <- length(mti$coef.names)
    if(j==0) next # Auxiliary: no parameters or statistics.
    if(model$offset[i]){
     etamap$offsetmap <- c(etamap$offsetmap, rep(TRUE,j))
    }else{
     etamap$offsetmap <- c(etamap$offsetmap, rep(FALSE,j) | NVL(mti$offsetmap,FALSE))
    }
    if (is.null(mti$params)) { # Not a curved parameter
      etamap$canonical <- c(etamap$canonical, to:(to+j-1L))
      from <- from+j
      to <- to+j
      if(model$offset[i]){
       etamap$offsettheta <- c(etamap$offsettheta, rep(TRUE,j))
      }else{
       etamap$offsettheta <- c(etamap$offsettheta, rep(FALSE,j) | NVL(mti$offsettheta,FALSE))
      }
      
      etamap$mintheta <- c(etamap$mintheta,
                           rep(NVL(mti$minpar, -Inf), length.out=j))
      etamap$maxtheta <- c(etamap$maxtheta,
                           rep(NVL(mti$maxpar, +Inf), length.out=j))
    } else { # curved parameter
      k <- length(mti$params)
      etamap$canonical <- c(etamap$canonical, rep(0L, k))
      etamap$curved[[a]] <- list(from=from+seq_len(k)-1L,
                                 to=to:(to+j-1L),
                                 map=mti$map, gradient=mti$gradient,
                                 cov=mti$eta.cov)  #Added by CTB 1/28/06
      from <- from+k
      to <- to+j
      a <- a+1L
      if(model$offset[i]){
       etamap$offsettheta <- c(etamap$offsettheta, rep(TRUE,k))
      }else{
       etamap$offsettheta <- c(etamap$offsettheta, rep(FALSE,k) | NVL(mti$offsettheta,FALSE))
      }
      
      etamap$mintheta <- c(etamap$mintheta,
                           rep(NVL(mti$minpar, -Inf), length.out=k))
      etamap$maxtheta <- c(etamap$maxtheta,
                           rep(NVL(mti$maxpar, +Inf), length.out=k))
    }
  }
  etamap$etalength <- to-1L
  etamap
} 
