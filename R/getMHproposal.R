MHproposals<-
  #         Class Constraints     Reference              Weights        MHP
  rbind(I(c("c", "",              "Bernoulli",  "default",      "TNT")),
          c("c", "",              "Bernoulli",  "TNT",          "TNT"),
          c("c", "",              "Bernoulli",  "random",       "randomtoggle"),
          c("c", "",              "Bernoulli",  "TNT10",        "TNT10"),
          c("c", "bd",            "Bernoulli",  "default",       "TNT"),
          c("c", "bd",            "Bernoulli",  "TNT",           "TNT"),
          c("c", "bd",            "Bernoulli",  "random",       "randomtoggle"),
          c("c", "bd+edges",      "Bernoulli",  "default",      "ConstantEdges"),
          c("c", "bd+edges",      "Bernoulli",      "random",       "ConstantEdges"),          
          c("c", "degrees",       "Bernoulli",       "default",      "CondDegree"),
          c("c", "degrees",       "Bernoulli",       "random",       "CondDegree"),
          c("c", "degreesTetrad", "Bernoulli",       "default",      "CondDegreeTetradToggles"),
          c("c", "degreesTetrad", "Bernoulli",       "random",       "CondDegreeTetradToggles"),
          c("c", "degreesHexad",  "Bernoulli",       "default",      "CondDegreeHexadToggles"),
          c("c", "degreesHexad",  "Bernoulli",       "random",       "CondDegreeHexadToggles"),

          c("c", "degreedist",    "Bernoulli",    "default",      "CondDegreeDist"),
          c("c", "degreedist",    "Bernoulli",    "random",       "CondDegreeDist"), 
          c("c", "indegreedist",  "Bernoulli",  "default",      "CondInDegreeDist"),
          c("c", "indegreedist",  "Bernoulli",  "random",       "CondInDegreeDist"), 
          c("c", "outdegreedist", "Bernoulli",  "default",      "CondOutDegreeDist"),
          c("c", "outdegreedist", "Bernoulli",  "random",       "CondOutDegreeDist"), 
#          c("c", "indegrees",    "Bernoulli",     "default",      "CondInDegree"),
#          c("c", "indegrees",    "Bernoulli",     "random",       "CondInDegree"),
#          c("c", "outdegrees",   "Bernoulli",    "default",      "CondOutDegree"),
#          c("c", "outdegrees",   "Bernoulli",    "random",       "CondOutDegree"),
          c("c", "edges",         "Bernoulli",         "default",      "ConstantEdges"),
          c("c", "edges",         "Bernoulli",         "random",       "ConstantEdges"),
          c("c", "hamming",       "Bernoulli",       "default",      "HammingTNT"),
          c("c", "hamming",       "Bernoulli",       "random",       "HammingTNT"),
          c("c", "edges+hamming", "Bernoulli", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "Bernoulli", "random",       "HammingConstantEdges"),
          c("c", "observed",      "Bernoulli",      "default",      "randomtoggleNonObserved"),
          c("c", "observed",      "Bernoulli",      "random",       "randomtoggleNonObserved"),
          c("f", "",              "Bernoulli",              "default",      "formationTNT"),
          c("f", "",              "Bernoulli",              "TNT",          "formationTNT"),
          c("f", "",              "Bernoulli",              "random",       "formation"),
          c("f", "bd",            "Bernoulli",            "default",      "formationTNT"),
          c("f", "bd",            "Bernoulli",            "TNT",          "formationTNT"),
          c("f", "bd",            "Bernoulli",            "random",       "formation"),
          c("d", "",              "Bernoulli",              "default",      "dissolution"),
          c("d", "",              "Bernoulli",              "random",       "dissolution"),
          c("d", "bd",            "Bernoulli",            "default",      "dissolution"),
          c("d", "bd",            "Bernoulli",            "random",       "dissolution"),
          c("fmle", "",           "Bernoulli",            "default",      "formationMLE"),
          c("fmle", "",           "Bernoulli",            "random",       "formationMLE"),
          c("fmle", "observed",   "Bernoulli",            "default",      "formationNonObservedMLE"),
          c("fmle", "observed",   "Bernoulli",            "random",       "formationNonObservedMLE"),
          c("fmle", "bd",         "Bernoulli",            "default",      "formationMLE"),
          c("fmle", "bd",         "Bernoulli",            "random",       "formationMLE"),
          c("dmle", "",           "Bernoulli",            "default",      "dissolutionMLE"),
          c("dmle", "",           "Bernoulli",            "random",       "dissolutionMLE"),
          c("dmle", "observed",   "Bernoulli",            "default",      "dissolutionNonObservedMLE"),
          c("dmle", "observed",   "Bernoulli",            "random",       "dissolutionNonObservedMLE"),
          c("dmle", "bd",         "Bernoulli",            "default",      "dissolutionMLE"),
          c("dmle", "bd",         "Bernoulli",            "random",       "dissolutionMLE"),
          c("c", "",              "Poisson",  "default",      "Poisson"),
          c("c", "",              "Poisson",  "random",       "Poisson")
        )
MHproposals <- data.frame(I(MHproposals[,1]), I(MHproposals[,2]), 
                          I(MHproposals[,3]), I(MHproposals[,4]),
                          I(MHproposals[,5]))
colnames(MHproposals)<-c("Class","Constraints","Reference","Weights","MHP")

MHproposal<-function(object, ...) UseMethod("MHproposal")

# This could be useful for trapping bugs before they become mysterious
# segfaults.
MHproposal.NULL<-function(object, ...) stop("NULL passed to MHproposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.")

MHproposal.MHproposal<-function(object,...) return(object)

MHproposal.character <- function(object, arguments, nw, model, ...){
  name<-object
  proposal <- eval(call(paste("InitMHP.", name, sep=""),
                        arguments, nw, model))

  proposal$bd<-ergm.bounddeg(arguments$bd,nw)

  class(proposal)<-"MHproposal"
  proposal
}

MHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c", reference="Bernoulli", ...) {

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
      conlist <- try(eval(as.call(init.call), attr(constraints,".Environment")))
      if(inherits(conlist,"try-error")){
        stop(paste("The constraint you have selected ('",constraints,"') does not exist in 'ergm'. Are you sure you have not mistyped it?",sep=""))
      }
    }
  }
  ## Remove constraints implied by other constraints.
  for(constr in names(conlist))
    for(impl in ConstraintImplications[[constr]])
      conlist[[impl]]<-NULL

  ## Convert vector of constraints to a "standard form".
  if(is.null(names(conlist))) {
    constraints <- ""
  } else {
    constraints <- paste(sort(tolower(names(conlist))),collapse="+")
  }
  name<-with(MHproposals,MHP[Class==class & Constraints==constraints & Reference==reference & Weights==weights])
  if(length(name)>1) stop("Multiple matching proposals in the lookup table.",
                          "This Should Not Be Happening (tm). Unless you have",
                          "been messing with the table, please file a bug report.")
  ## TODO: Get intelligent error messages for reference mismatches as well.
  if(length(name)<1){
    constraints<-with(MHproposals,Constraints[Class==class & Weights==weights])
    weightings<-with(MHproposals,Weights[Class==class & Constraints==constraints])
    stop("This combination of model constraint and proposal weighting is not implemented. ",
         "Check your arguments for typos. \n",
         if(length(constraints)) paste("Constraints that go with your selected weighting are as follows: ",
                                       paste(constraints,collapse=", "),".\n",sep="")
         else "The supplied weighting is not recognized/implemented.\n ",
         if(length(weightings)) paste("Weightings that go with your selected constraint are as follows: ",
                                      paste(weightings,collapse=", "),".\n",sep="")
         else "The supplied constraint is not recognized/implemented.\n "
         )
  }
  if(is.null(arguments)) arguments<-conlist
  ## Hand it off to the class character method.
  MHproposal.character(name,arguments,nw,model)
}

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, model=NULL,weights=NULL,class="c", reference="Bernoulli"){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$prop.args
  if(is.null(nw)) nw<-object$network
  if(is.null(weights)) weights<-"default"
  if(is.null(model)){
    model<-if(class %in% c("c","f"))
      ergm.getmodel(object$formula,nw,...)
    else
      ergm.getmodel.dissolve(object$formula,nw,...)
  }  
  MHproposal(constraints,arguments=arguments,nw=nw,model=model,weights=weights,class=class,reference=reference)
}

