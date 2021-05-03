#  File R/ergm.etamap.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
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
  etamap <- list(canonical = NULL, offsetmap=NULL,
                 offsettheta=NULL, curved=list(), etalength=0L)
  from <- 1L
  to <- 1L
  a <- 1L
  if(length(model$terms)==0L) return(etamap)

  for (i in seq_along(model$terms)) {
    mti <- model$terms[[i]]
    j <- length(mti$coef.names)
    k <- NVL3(mti$params, length(.), j)

    if(j==0L) next # Auxiliary: no parameters or statistics.

    offset <- mti$offset
    etamap$offsettheta <- c(etamap$offsettheta, offset)

    if (is.null(mti$params)) { # Not a curved parameter
      etamap$canonical <- c(etamap$canonical, to:(to+j-1L))
      from <- from+j
      to <- to+j

      etamap$mintheta <- c(etamap$mintheta,
                           rep(NVL(mti$minpar, -Inf), length.out=j))
      etamap$maxtheta <- c(etamap$maxtheta,
                           rep(NVL(mti$maxpar, +Inf), length.out=j))

      etamap$offsetmap <- c(etamap$offsetmap, offset)
    } else { # curved parameter
      etamap$canonical <- c(etamap$canonical, rep(0L, k))
      etamap$curved[[a]] <- list(from=from+seq_len(k)-1L,
                                 to=to:(to+j-1L),
                                 map=mti$map, gradient=mti$gradient,
                                 cov=mti$eta.cov)  #Added by CTB 1/28/06
      from <- from+k
      to <- to+j
      a <- a+1L

      etamap$mintheta <- c(etamap$mintheta,
                           rep(NVL(mti$minpar, -Inf), length.out=k))
      etamap$maxtheta <- c(etamap$maxtheta,
                           rep(NVL(mti$maxpar, +Inf), length.out=k))

      etamap$offsetmap <- c(etamap$offsetmap, rep(all(offset),j)) # In a curved model, only set canonical parameters as offsets if *all* model parameters are offsets.
    }
  }
  etamap$etalength <- to-1L
  etamap
} 
