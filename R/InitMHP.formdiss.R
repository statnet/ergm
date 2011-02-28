#===================================================================
# This file contains the 5 following MHP initializers, each
# prepended with 'InitMHP.'  All of these functions may also be
# found in the <InitMHP> file.
#      <formation>       <formationTNT>
#      <formationMLE>    <dissolution>
#      <dissolutionMLE>
#===================================================================


InitMHP.formation <- function(arguments, nw, model) {
  MHproposal <- list(name = "Formation", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormation"
  }
  MHproposal
}


InitMHP.formationMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationMLE", inputs=NULL, package="ergm")
  MHproposal
}


InitMHP.dissolutionMLE <- function(arguments, nw, model) {
  MHproposal <- list(name = "DissolutionMLE", inputs=NULL, package="ergm")
  MHproposal
}


InitMHP.formationTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationTNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormationTNT"
  }
  MHproposal
}


InitMHP.dissolution <- function(arguments, nw, model) {
  MHproposal <- list(name = "Dissolution", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteDissolution"
  }
  MHproposal
}



