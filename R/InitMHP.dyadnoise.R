#  File R/InitMHP.dyadnoise.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
InitMHP.dyadnoiseTNT<-function(arguments, nw){
  p0to1 <- arguments$constraints$dyadnoise$p01
  p1to0 <- arguments$constraints$dyadnoise$p10
  p1to1 <- 1-p1to0
  p0to0 <- 1-p0to1

  deinf <- function(x, replace=1/.Machine$double.eps) ifelse(is.nan(x) | abs(x)<replace, x, sign(x)*replace)
  
  MHproposal <- list(name = if(length(p0to1)==1) "dyadnoiseTNT" else "dyadnoisemTNT", inputs=c(
                                                                                deinf(log(p1to0)-log(p0to0)), # Observed 0, State 0
                                                                                deinf(log(p0to0)-log(p1to0)), # Observed 0, State 1
                                                                                deinf(log(p1to1)-log(p0to1)), # Observed 1, State 0
                                                                                deinf(log(p0to1)-log(p1to1)), # Observed 1, State 1
                                                                                to_ergm_Cdouble(nw)))

  MHproposal                   
}

InitMHP.dyadnoise<-function(arguments, nw){
  p0to1 <- arguments$constraints$dyadnoise$p01
  p1to0 <- arguments$constraints$dyadnoise$p10
  p1to1 <- 1-p1to0
  p0to0 <- 1-p0to1

  deinf <- function(x, replace=1/.Machine$double.eps) ifelse(is.nan(x) | abs(x)<replace, x, sign(x)*replace)
  
  MHproposal <- list(name = if(length(p0to1)==1) "dyadnoise" else "dyadnoisem", inputs=c(
                                                                                deinf(log(p1to0)-log(p0to0)), # Observed 0, State 0
                                                                                deinf(log(p0to0)-log(p1to0)), # Observed 0, State 1
                                                                                deinf(log(p1to1)-log(p0to1)), # Observed 1, State 0
                                                                                deinf(log(p0to1)-log(p1to1)), # Observed 1, State 1
                                                                                to_ergm_Cdouble(nw)))

  MHproposal                   
}
