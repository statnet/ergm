#===========================================================================
# The <InitMHP> file contains the following 24 functions for
# initializing the MHproposal object; each is prepended with 'InitMHP.'
#        <dissolutionMLE>
#        <formationNonObservedMLE>
#        <dissolutionNonObservedMLE>
#        <formationMLE>       
#============================================================================


########################################################################
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



InitMHP.formationMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationMLE", inputs=NULL, package="ergm")
  MHproposal
}
InitMHP.dissolutionMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "DissolutionMLE", inputs=NULL, package="ergm")
  MHproposal
}
InitMHP.formationNonObservedMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationNonObservedMLE", inputs=ergm.Cprepare.miss(nw), package="ergm")
  MHproposal
}
InitMHP.dissolutionNonObservedMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "DissolutionNonObservedMLE", inputs=ergm.Cprepare.miss(nw), package="ergm")
  MHproposal
}


