#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

getMHproposal <- function(name, arguments, nw, model) {
  proposal <- eval(call(paste("InitMHP.", name, sep=""),
                        arguments, nw, model))
#  nwlist <- proposal$nwlist
proposal
}


nativeMHproposals<-
  #         Class Constraint       Weights        MHP
  rbind(I(c("c", "none",          "default",      "TNT")),
          c("c", "none",          "TNT",          "TNT"),
          c("c", "none",          "random",       "randomtoggle"),
          c("c", "none",          "nonobserved",  "randomtoggleNonObserved"),
          c("c", "degrees",       "default",      "CondDegree"),
          c("c", "degrees",       "random",       "CondDegree"),
          c("c", "nodedegrees",   "default",       "CondDegree"),
          c("c", "nodedegrees",   "random",       "CondDegree"),

          c("c", "degreedist",    "default",      "CondDegree"),
          c("c", "degreedist",    "random",       "CondDegree"),

          c("c", "indegrees",     "default",      "CondInDegree"),
          c("c", "indegrees",     "random",       "CondInDegree"),
          c("c", "outdegrees",    "default",      "CondOutDegree"),
          c("c", "outdegrees",    "random",       "CondOutDegree"),
          c("c", "edges",         "default",      "ConstantEdges"),
          c("c", "edges",         "random",       "ConstantEdges"),
          c("c", "hamming",       "default",      "Hamming"),
          c("c", "hamming",       "random",       "Hamming"),
          c("c", "edges+hamming", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "random",       "HammingConstantEdges"),
          c("f", "none",          "default",      "formationTNT"),
          c("f", "none",          "TNT",          "formationTNT"),
          c("f", "none",          "random",       "formation"),
          c("d", "none",          "default",      "dissolution"),
          c("d", "none",          "random",       "dissolution"))
tmp <- nativeMHproposals
nativeMHproposals <- data.frame(I(tmp[,1]), I(tmp[,2]), I(tmp[,3]), I(tmp[,4]))  
colnames(nativeMHproposals)<-c("Class","Constraint","Weights","MHP")

lookupMHproposal <- function(class, constraint, weights){
  if(is.null(constraint)) constraint<-"none"
  # If the proposal is specified explicitly, just use it.
  if(length(grep(":",constraint))) return(gsub(".*:","",constraint))
  
  # Convert vector of constraints to a "standard form".
  constraint<-paste(sort(tolower(constraint)),collapse="+")
  MHP<-with(nativeMHproposals,MHP[Class==class & Constraint==constraint & Weights==weights])
  if(length(MHP)>1) stop("Multiple matching proposals in the lookup table.",
                         "This Should Not Be Happening. Please make a bug report.")
  if(length(MHP)<1){
    constraints<-with(nativeMHproposals,Constraint[Class==class & Weights==weights])
    weightings<-with(nativeMHproposals,Weights[Class==class & Constraint==constraint])
    stop("This combination of model constraint and proposal weighting is not implemented.",
         "Check your arguments for typos.",
         if(length(constraints)) paste("Constraints that go with your selected weighting are as follows: ",
               paste(constraints,collapse=", "),".",sep="")
         else "The supplied weighting is not recognized/implemented.",
         if(length(weightings)) paste("Weightings that go with your selected constraint are as follows: ",
               paste(weightings,collapse=", "),".",sep="")
         else "The supplied constraint is not recognized/implemented."
         )
  }
  MHP
}

forceMH<-function(proposaltype){
  paste(":",proposaltype,sep="")
}

lookupMH.ergm<-function(object,class,constraint=NULL,weights){
  if(is.null(constraint) & !is.null(object$proposaltype))
    object$proposaltype
  else
    lookupMHproposal(class,constraint,weights)
}
