#  File R/InitErgmProposal.dyadnoise.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################

#' @templateVar name dyadnoiseTNT
#' @aliases InitErgmProposal.dyadnoiseTNT
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitErgmProposal.dyadnoiseTNT<-function(arguments, nw){
  p0to1 <- arguments$constraints$dyadnoise$p01
  p1to0 <- arguments$constraints$dyadnoise$p10
  p1to1 <- 1-p1to0
  p0to0 <- 1-p0to1

  list(name = if(length(p0to1)==1) "dyadnoiseTNT" else "dyadnoisemTNT", inputs=c(
                                                                          deInf(log(p1to0)-log(p0to0)), # Observed 0, State 0
                                                                          deInf(log(p0to0)-log(p1to0)), # Observed 0, State 1
                                                                          deInf(log(p1to1)-log(p0to1)), # Observed 1, State 0
                                                                          deInf(log(p0to1)-log(p1to1)), # Observed 1, State 1
                                                                          to_ergm_Cdouble(nw)),
       bd = ergm_bd_init(arguments, nw))
}

#' @templateVar name dyadnoise
#' @aliases InitErgmProposal.dyadnoise
#' @title TODO
#' @description TODO
#' @template ergmProposal-general
NULL
InitErgmProposal.dyadnoise<-function(arguments, nw){
  p0to1 <- arguments$constraints$dyadnoise$p01
  p1to0 <- arguments$constraints$dyadnoise$p10
  p1to1 <- 1-p1to0
  p0to0 <- 1-p0to1

  list(name = if(length(p0to1)==1) "dyadnoise" else "dyadnoisem", inputs=c(
                                                                                deInf(log(p1to0)-log(p0to0)), # Observed 0, State 0
                                                                                deInf(log(p0to0)-log(p1to0)), # Observed 0, State 1
                                                                                deInf(log(p1to1)-log(p0to1)), # Observed 1, State 0
                                                                                deInf(log(p0to1)-log(p1to1)), # Observed 1, State 1
                                                                    to_ergm_Cdouble(nw)),
       bd = ergm_bd_init(arguments, nw))
}
