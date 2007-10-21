#  DH:  This file is still a work in progress, but for now it should
#  be set up so as not to break anything!

MHproposals<-
  #         Class Constraint      Weights        MHP
  rbind(I(c("c", "",              "default",      "TNT")),
          c("c", "",              "TNT",          "TNT"),
          c("c", "",              "random",       "randomtoggle"),
          c("c", "bd",            "default",       "TNT"),
          c("c", "bd",            "TNT",           "TNT"),
          c("c", "bd",            "random",       "randomtoggle"),
          c("c", "",              "nonobserved",  "randomtoggleNonObserved"),
          c("c", "degrees",       "default",      "CondDegree"),
          c("c", "degrees",       "random",       "CondDegree"),

          c("c", "degreedist",    "default",      "CondDegreeDist"), # not yet implemented
          c("c", "degreedist",    "random",       "CondDegreeDist"), # not yet implemented 

          c("c", "indegrees",     "default",      "CondInDegree"),
          c("c", "indegrees",     "random",       "CondInDegree"),
          c("c", "outdegrees",    "default",      "CondOutDegree"),
          c("c", "outdegrees",    "random",       "CondOutDegree"),
          c("c", "edges",         "default",      "ConstantEdges"),
          c("c", "edges",         "random",       "ConstantEdges"),
          c("c", "hamming",       "default",      "HammingTNT"),
          c("c", "hamming",       "random",       "HammingTNT"),
          c("c", "edges+hamming", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "random",       "HammingConstantEdges"),
          c("f", "",              "default",      "formationTNT"),
          c("f", "",              "TNT",          "formationTNT"),
          c("f", "",              "random",       "formation"),
          c("f", "bd",            "default",      "formationTNT"),
          c("f", "bd",            "TNT",          "formationTNT"),
          c("f", "bd",            "random",       "formation"),
          c("d", "",              "default",      "dissolution"),
          c("d", "",              "random",       "dissolution"),
          c("d", "bd",            "default",      "dissolution"),
          c("d", "bd",            "random",       "dissolution")
        )
tmp <- MHproposals
MHproposals <- data.frame(I(tmp[,1]), I(tmp[,2]), I(tmp[,3]), I(tmp[,4]))  
colnames(MHproposals)<-c("Class","Constraints","Weights","MHP")


getMHproposal<-function(object, arguments=NULL, nw=NULL, model=NULL, weights="default", ...) UseMethod("getMHproposal")

getMHproposal.NULL<-function(object,...) getMHproposal(~.,...)

getMHproposal.MHproposal<-function(object,...) return(object)

getMHproposal.character <- function(object, arguments, nw, model){
  name<-object
  proposal <- eval(call(paste("InitMHP.", name, sep=""),
                        arguments, nw, model))

  proposal$bd<-ergm.boundDeg(arguments$bd,network.size(nw))

  class(proposal)<-"MHproposal"
  proposal
}

getMHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c") {

  constraints<-object
  ## Construct a list of constraints and arguments from the formula.
  conlist<-list()
  constraints<-try(as.list(attr(terms(constraints),"variables"))[-1],silent=TRUE)
  ## The . in the default formula will break terms(...), signaling no constraints. 
  if(!inherits(constraints,"try-error")){
    for(constraint in constraints){
      if(is.call(constraint)){
        init.call<-list()
        init.call<-list(as.name(paste("InitConstraint.", constraint[[1]], sep = "")),
                        conlist=conlist)
        
        init.call<-c(init.call,as.list(constraint)[-1])
      }else{
        init.call <- list(as.name(paste("InitConstraint.", constraint, sep = "")),conlist=conlist)
      }
      conlist <- eval(as.call(init.call), attr(constraints,".Environment"))
    }
  }
  
  ## Remove constraints implied by other constraints.
  for(constr in names(conlist))
    for(impl in ConstraintImplications[[constr]])
      conlist[[impl]]<-NULL
  
  ## Convert vector of constraints to a "standard form".
  constraints<-paste(sort(tolower(names(conlist))),collapse="+")
  name<-with(MHproposals,MHP[Class==class & Constraints==constraints & Weights==weights])
  if(length(name)>1) stop("Multiple matching proposals in the lookup table.",
                          "This Should Not Be Happening (tm). Unless you have",
                          "been messing with the table, please file a bug report.")
  if(length(name)<1){
    constraints<-with(MHproposals,Constraint[Class==class & Weights==weights])
    weightings<-with(MHproposals,Weights[Class==class & Constraint==constraint])
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
  
  if(is.null(arguments)) arguments<-conlist

  ## Hand it off to the class character method.
  getMHproposal.character(name,arguments,nw,model)
}

getMHproposal.ergm<-function(object,arguments=NULL, nw=NULL, model=NULL,weights=NULL){
  if("proposal" %in% names(object)) object$proposal
  else if(!is.null(object$proposal)){ # (Ab)use list name extension.
    if(missing(arguments)||missing(model)) warning("Possibly guessing proposal type and arguments.")
    getMHproposal(constraints=object$proposal,arguments=arguments,nw=object$network,model=list(),class="c")
  }
  else NULL
}

