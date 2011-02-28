#=======================================================================================
# This file contains the following 6 files for creating MHproposal objects
#          <MHproposal>                <MHproposal.character>
#          <MHproposal.NULL>           <MHproposal.formula>
#          <MHproposal.MHproposal>     <MHproposal.ergm>
#=======================================================================================



# the MHproposal look-up table: these are the combinations of constraints, weights,
# and proposals that are supported

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
          c("c", "bd+observed",   "Bernoulli",      "default",      "randomtoggleNonObserved"),
          c("c", "bd+observed",   "Bernoulli",      "random",       "randomtoggleNonObserved"),
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
          c("c", "",              "Poisson",  "random",       "Poisson"),
          c("c", "observed",      "Poisson",  "default",      "PoissonNonObserved"),
          c("c", "observed",      "Poisson",  "random",       "PoissonNonObserved"),
          c("c", "",              "DescRank",  "default",      "DescRank"),
          c("c", "",              "DescRank",  "random",       "DescRank"),
          c("c", "",              "StdNormal",  "default",      "StdNormal"),
          c("c", "",              "StdNormal",  "random",       "StdNormal"),
          c("c", "ranks",         "StdNormal",  "default",      "StdNormalRank"),
          c("c", "ranks",         "StdNormal",  "random",       "StdNormalRank")

        )
MHproposals <- data.frame(I(MHproposals[,1]), I(MHproposals[,2]), 
                          I(MHproposals[,3]), I(MHproposals[,4]),
                          I(MHproposals[,5]))
colnames(MHproposals)<-c("Class","Constraints","Reference","Weights","MHP")



########################################################################################
# The <MHproposal> function initializes and returns an MHproposal object via one of the
# class-specific functions below
#
# --PARAMETERS--
#   (see the class-specific function headers)
#
# --RETURNED--
#   proposal: an MHproposal object as a list containing
#   name   : the name of the proposal
#   args   : NULL (I think - the only non-null value returned by the InitMH
#            is for <nobetweengroupties>, but this isn't included in the 
#            look-up table
#   package: is "ergm"
#   bd     : the list of parameters to bound degree in the fitting process
#            and returned by <ergm.bounddeg>
#
########################################################################################

MHproposal<-function(object, ...) UseMethod("MHproposal")


# This could be useful for trapping bugs before they become mysterious segfaults.
MHproposal.NULL<-function(object, ...) stop("NULL passed to MHproposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.")


MHproposal.MHproposal<-function(object,...) return(object)




########################################################################################
# The <MHproposal.character> function initializes the MHproposal object using the
# <InitMHP.> function that corresponds to the name given in 'object'
#
# --PARAMETERS--
#   object     :  the name of the proposal, one found in the look-up table
#   arguments  :  a list of parameters used by the <Init.MHP> routines possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  the network object orginally given to <ergm> via 'formula'
#   model      :  the initial model object constructed by <ergm>
#
########################################################################################

MHproposal.character <- function(object, arguments, nw, model, ..., response=NULL){
  name<-object
  proposal <- {
    if(is.null(response))
      eval(call(paste("InitMHP", name, sep="."),
                arguments, nw, model))
    else
      eval(call(paste("InitWtMHP", name, sep="."),
                arguments, nw, model, response))
  }

  proposal$bd<-ergm.bounddeg(arguments$bd,nw)

  class(proposal)<-"MHproposal"
  proposal
}





########################################################################################
# The <MHproposal.formula> function verifies that the given constraints exist and
# are supported in conjuction with the given weights and class by a unique MH proposal;
# if so the MHproposal object is created via <MHproposal.character> using the 
# MHP type found in the look-up table above.
#
# --PARAMETERS--
#   object     :  a one-sided formula of constraint terms ( ~ term(s))
#   arguments  :  a list of parameters used by the <Init.MHP> routines  possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  a network object
#   model      :  a model object; default=<ergm.getmodel(object$formula,nw,...)>  
#   constraints:  the constraints as a one sided formula '~ term(s)'
#   weights    :  specifies the method used to allocate probabilities of being proposed
#                 to dyads; options are "TNT", "TNT10", "random", "nonobserved" and
#                 "default"; default="default"
#   class      :  the class of the proposal; choices include "c", "f", "d", "fmle"
#                 and "dmle"; default="c"
#
########################################################################################

MHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c", reference="Bernoulli", response=NULL, ...) {

  constraints<-object
  ## Construct a list of constraints and arguments from the formula.
  conlist<-list()
  constraints<-as.list(attr(terms(constraints,allowDotAsName=TRUE),"variables"))[-1]
  for(constraint in constraints){
    ## The . in the default formula means no constrains.
    ## There may be other constraints in the formula, however.
    if(constraint==".") next
    
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
  MHproposal.character(name,arguments,nw,model,response=response)
}





########################################################################################
# The <MHproposal.ergm> function creates the MHproposal object via <MHproposal.formula>
# after extracting the appropriate parameters from the given ergm object
#
# --PARAMETERS--
#   object     :  an ergm object
#   ...        :  parameters used to create the model via <ergm.getmodel>;
#                 only used if 'model' is not specified; these may include
#                 'silent', 'drop' and 'initialfit'
#   constraints:  the constraints as a one sided formula '~ term(s)';
#                 default=object$constraints
#   arguments  :  a list of parameters used by the <Init.MHP> routines  possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  a network object; default=object.network
#   model      :  a model object; default=<ergm.getmodel(object$formula,nw,...)>
#   weights    :  the proposal weights component of <control.ergm> as either
#                 "TNT", "random", "TNT10", or "default"; default="default"
#                 (these options don't agree with the prop.weights of <control.ergm>)
#   class      :  "c", otherwise execution will halt
#
########################################################################################

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, model=NULL,weights=NULL,class="c", reference="Bernoulli", response=NULL){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$prop.args
  if(is.null(nw)) nw<-object$network
  if(is.null(response)) response<-object$response
  if(is.null(weights)) weights<-"default"
  if(is.null(model)){
    model<-if(class %in% c("c","f"))
      ergm.getmodel(object$formula,nw,response=response,...)
    else
      ergm.getmodel.dissolve(object$formula,nw,response=response,...)
  }  
  MHproposal(constraints,arguments=arguments,nw=nw,model=model,weights=weights,class=class,reference=reference,response=response)
}

