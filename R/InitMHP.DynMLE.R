#  File ergm/R/InitMHP.DynMLE.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
##########################################2##############################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
############################################################################



InitMHP.formationMLE <- function(arguments, nw) {
  MHproposal <- list(name = "FormationMLE", inputs=ergm.Cprepare.el(arguments$constraints$atleast$nw), package="ergm")
  MHproposal
}

InitMHP.formationMLETNT <- function(arguments, nw) {
  MHproposal <- list(name = "FormationMLETNT", inputs=ergm.Cprepare.el(arguments$constraints$atleast$nw), package="ergm")
  MHproposal
}

InitMHP.dissolutionMLE <- function(arguments, nw) {
  MHproposal <- list(name = "DissolutionMLE", inputs=ergm.Cprepare.el(arguments$constraints$atmost$nw), package="ergm")
  MHproposal
}

InitMHP.formationNonObservedMLE <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are non-ties in y[t-1]
  
  y0<-arguments$costraints$atleast$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  MHproposal <- list(name = "randomtoggleNonObserved", inputs=ergm.Cprepare.el(y.miss-y0), package="ergm")
  MHproposal
}

InitMHP.dissolutionNonObservedMLE <- function(arguments, nw) {
  ## Precalculate toggleable dyads: dyads which
  ## * are unobserved in y[t]
  ## * are ties in y[t-1]

  y0<-arguments$constraints$atmost$nw
  y.miss<-is.na(nw)

  ## Given the list of toggleable dyads, no formation-specific proposal function is needed:
  MHproposal <- list(name = "randomtoggleNonObserved", inputs=ergm.Cprepare.el(y.miss & y0), package="ergm")
  MHproposal
}


