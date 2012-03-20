#===================================================================
# This file contains the 5 following MHP initializers, each
# prepended with 'InitMHP.'  All of these functions may also be
# found in the <InitMHP> file.
#      <formation>       <formationTNT>
#      <dissolution>
#===================================================================

InitMHP.formation <- function(arguments, nw, model) {
  MHproposal <- list(name = "Formation", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormation"
  }
  MHproposal
}
#ergm.MHP.table("f", "Bernoulli", "",  0, "random", "formation")
#ergm.MHP.table("f", "Bernoulli", "bd",  0, "random", "formation")

InitMHP.formationTNT <- function(arguments, nw, model) {
  MHproposal <- list(name = "FormationTNT", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteFormationTNT"
  }
  MHproposal
}
#ergm.MHP.table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
#ergm.MHP.table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")

InitMHP.dissolution <- function(arguments, nw, model) {
  MHproposal <- list(name = "Dissolution", inputs=NULL, package="ergm")
  if(is.bipartite(nw)){
    MHproposal$name <- "BipartiteDissolution"
  }
  MHproposal
}
#ergm.MHP.table("d", "Bernoulli", "",  0, "random", "dissolution")
#ergm.MHP.table("d", "Bernoulli", "bd",  0, "random", "dissolution")
