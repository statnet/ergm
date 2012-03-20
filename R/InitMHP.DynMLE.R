#===========================================================================
# The <InitMHP> file contains the following 24 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#        <dissolutionMLE>
#        <formationNonObservedMLE>
#        <dissolutionNonObservedMLE>
#        <formationMLE>       
#============================================================================


##########################################2##############################
# Each of the <InitMHP.X> functions initializes and returns a
# MHproposal list; when appropriate, proposal types are checked against
# covariates and network types for 1 of 2 side effects: to print warning
# messages or to halt execution (only <InitMHP.nobetweengroupties> can
# halt execution)
#
# --PARAMETERS--
#   arguments: is ignored by all but <InitMHP.nobetweengroupties>,
#              where 'arguments' is used to get the nodal attributes
#              via <get.node.attr>
#   nw       : the network given by the model
#   model    : the model for 'nw', as returned by <ergm.getmodel>
#
# --RETURNED--
#   MHproposal: a list containing:
#        name   : the name of the proposal
#        inputs : a vector to be passed to the proposal
#        package: is "ergm"
#
############################################################################

InitMHP.formationMLE <- function(arguments, nw) {
  MHproposal <- list(name = "FormationMLE", inputs=ergm.Cprepare.el(arguments$constraints$atleast$nw), package="ergm")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "atleast",  0, "random", "formationMLE")
#ergm.MHP.table("c", "Bernoulli", "atleast+bd",  0, "random", "formationMLE")


InitMHP.formationMLETNT <- function(arguments, nw) {
  MHproposal <- list(name = "FormationMLETNT", inputs=ergm.Cprepare.el(arguments$constraints$atleast$nw), package="ergm")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "atleast",  1, "TNT", "formationMLETNT")
#ergm.MHP.table("c", "Bernoulli", "atleast+bd",  1, "TNT", "formationMLETNT")

InitMHP.dissolutionMLE <- function(arguments, nw) {
  MHproposal <- list(name = "DissolutionMLE", inputs=ergm.Cprepare.el(arguments$constraints$atmost$nw), package="ergm")
  MHproposal
}
#ergm.MHP.table("c", "Bernoulli", "atmost",  0, "random", "dissolutionMLE")
#ergm.MHP.table("c", "Bernoulli", "atmost+bd",  0, "random", "dissolutionMLE")

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
#ergm.MHP.table("c", "Bernoulli", "atleast+observed",  0, "random", "formationNonObservedMLE")

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
#ergm.MHP.table("c", "Bernoulli", "atmost+observed",  0, "random", "dissolutionNonObservedMLE")

