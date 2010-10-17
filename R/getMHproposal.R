MHproposals<-
  #         Class Constraints     Family              Weights        MHP
  rbind(I(c("c", "",              "PseudoBernoulli.logit",  "default",      "TNT")),
          c("c", "",              "PseudoBernoulli.logit",  "TNT",          "TNT"),
          c("c", "",              "PseudoBernoulli.logit",  "random",       "randomtoggle"),
          c("c", "",              "PseudoBernoulli.logit",  "TNT10",        "TNT10"),
          c("c", "bd",            "PseudoBernoulli.logit",  "default",       "TNT"),
          c("c", "bd",            "PseudoBernoulli.logit",  "TNT",           "TNT"),
          c("c", "bd",            "PseudoBernoulli.logit",  "random",       "randomtoggle"),
          c("c", "bd+edges",      "PseudoBernoulli.logit",  "default",      "ConstantEdges"),
          c("c", "bd+edges",      "PseudoBernoulli.logit",      "random",       "ConstantEdges"),          
          c("c", "degrees",       "PseudoBernoulli.logit",       "default",      "CondDegree"),
          c("c", "degrees",       "PseudoBernoulli.logit",       "random",       "CondDegree"),
          c("c", "degreesTetrad", "PseudoBernoulli.logit",       "default",      "CondDegreeTetradToggles"),
          c("c", "degreesTetrad", "PseudoBernoulli.logit",       "random",       "CondDegreeTetradToggles"),
          c("c", "degreesHexad",  "PseudoBernoulli.logit",       "default",      "CondDegreeHexadToggles"),
          c("c", "degreesHexad",  "PseudoBernoulli.logit",       "random",       "CondDegreeHexadToggles"),

          c("c", "degreedist",    "PseudoBernoulli.logit",    "default",      "CondDegreeDist"),
          c("c", "degreedist",    "PseudoBernoulli.logit",    "random",       "CondDegreeDist"), 
          c("c", "indegreedist",  "PseudoBernoulli.logit",  "default",      "CondInDegreeDist"),
          c("c", "indegreedist",  "PseudoBernoulli.logit",  "random",       "CondInDegreeDist"), 
          c("c", "outdegreedist", "PseudoBernoulli.logit",  "default",      "CondOutDegreeDist"),
          c("c", "outdegreedist", "PseudoBernoulli.logit",  "random",       "CondOutDegreeDist"), 
#          c("c", "indegrees",    "PseudoBernoulli.logit",     "default",      "CondInDegree"),
#          c("c", "indegrees",    "PseudoBernoulli.logit",     "random",       "CondInDegree"),
#          c("c", "outdegrees",   "PseudoBernoulli.logit",    "default",      "CondOutDegree"),
#          c("c", "outdegrees",   "PseudoBernoulli.logit",    "random",       "CondOutDegree"),
          c("c", "edges",         "PseudoBernoulli.logit",         "default",      "ConstantEdges"),
          c("c", "edges",         "PseudoBernoulli.logit",         "random",       "ConstantEdges"),
          c("c", "hamming",       "PseudoBernoulli.logit",       "default",      "HammingTNT"),
          c("c", "hamming",       "PseudoBernoulli.logit",       "random",       "HammingTNT"),
          c("c", "edges+hamming", "PseudoBernoulli.logit", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "PseudoBernoulli.logit", "random",       "HammingConstantEdges"),
          c("c", "observed",      "PseudoBernoulli.logit",      "default",      "randomtoggleNonObserved"),
          c("c", "observed",      "PseudoBernoulli.logit",      "random",       "randomtoggleNonObserved"),
          c("f", "",              "PseudoBernoulli.logit",              "default",      "formationTNT"),
          c("f", "",              "PseudoBernoulli.logit",              "TNT",          "formationTNT"),
          c("f", "",              "PseudoBernoulli.logit",              "random",       "formation"),
          c("f", "bd",            "PseudoBernoulli.logit",            "default",      "formationTNT"),
          c("f", "bd",            "PseudoBernoulli.logit",            "TNT",          "formationTNT"),
          c("f", "bd",            "PseudoBernoulli.logit",            "random",       "formation"),
          c("d", "",              "PseudoBernoulli.logit",              "default",      "dissolution"),
          c("d", "",              "PseudoBernoulli.logit",              "random",       "dissolution"),
          c("d", "bd",            "PseudoBernoulli.logit",            "default",      "dissolution"),
          c("d", "bd",            "PseudoBernoulli.logit",            "random",       "dissolution"),
          c("fmle", "",           "PseudoBernoulli.logit",            "default",      "formationMLE"),
          c("fmle", "",           "PseudoBernoulli.logit",            "random",       "formationMLE"),
          c("fmle", "observed",   "PseudoBernoulli.logit",            "default",      "formationNonObservedMLE"),
          c("fmle", "observed",   "PseudoBernoulli.logit",            "random",       "formationNonObservedMLE"),
          c("fmle", "bd",         "PseudoBernoulli.logit",            "default",      "formationMLE"),
          c("fmle", "bd",         "PseudoBernoulli.logit",            "random",       "formationMLE"),
          c("dmle", "",           "PseudoBernoulli.logit",            "default",      "dissolutionMLE"),
          c("dmle", "",           "PseudoBernoulli.logit",            "random",       "dissolutionMLE"),
          c("dmle", "observed",   "PseudoBernoulli.logit",            "default",      "dissolutionNonObservedMLE"),
          c("dmle", "observed",   "PseudoBernoulli.logit",            "random",       "dissolutionNonObservedMLE"),
          c("dmle", "bd",         "PseudoBernoulli.logit",            "default",      "dissolutionMLE"),
          c("dmle", "bd",         "PseudoBernoulli.logit",            "random",       "dissolutionMLE"),
          c("c", "",              "PseudoPoisson.log",  "default",      "PseudoPoisson"),
          c("c", "",              "PseudoPoisson.log",  "random",       "PseudoPoisson")
        )
MHproposals <- data.frame(I(MHproposals[,1]), I(MHproposals[,2]), 
                          I(MHproposals[,3]), I(MHproposals[,4]),
                          I(MHproposals[,5]))
colnames(MHproposals)<-c("Class","Constraints","Family","Weights","MHP")

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

MHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c", family="PseudoBernoulli.logit", ...) {

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
  name<-with(MHproposals,MHP[Class==class & Constraints==constraints & Family==family & Weights==weights])
  if(length(name)>1) stop("Multiple matching proposals in the lookup table.",
                          "This Should Not Be Happening (tm). Unless you have",
                          "been messing with the table, please file a bug report.")
  ## TODO: Get intelligent error messages for family mismatches as well.
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

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, model=NULL,weights=NULL,class="c", family="PseudoBernoulli.logit"){
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
  MHproposal(constraints,arguments=arguments,nw=nw,model=model,weights=weights,class=class,family=family)
}

