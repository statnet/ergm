#=======================================================================================
# This file contains the following 6 files for creating MHproposal objects
#          <MHproposal>                <MHproposal.character>
#          <MHproposal.NULL>           <MHproposal.formula>
#          <MHproposal.MHproposal>     <MHproposal.ergm>
#=======================================================================================



# the MHproposal look-up table: these are the combinations of constraints, weights,
# and proposals that are supported

MHproposals<-
  #         Class Constraints      Weights        MHP
  rbind(I(c("c", "",              "default",      "TNT")),
          c("c", "",              "TNT",          "TNT"),
          c("c", "",              "random",       "randomtoggle"),
          c("c", "",              "TNT10",        "TNT10"),
          c("c", "bd",            "default",       "TNT"),
          c("c", "bd",            "TNT",           "TNT"),
          c("c", "bd",            "random",       "randomtoggle"),
          c("c", "bd+edges",      "default",      "ConstantEdges"),
          c("c", "bd+edges",      "random",       "ConstantEdges"),          
          c("c", "degrees",       "default",      "CondDegree"),
          c("c", "degrees",       "random",       "CondDegree"),
          c("c", "degreesTetrad",       "default",      "CondDegreeTetradToggles"),
          c("c", "degreesTetrad",       "random",       "CondDegreeTetradToggles"),
          c("c", "degreesHexad",       "default",      "CondDegreeHexadToggles"),
          c("c", "degreesHexad",       "random",       "CondDegreeHexadToggles"),

          c("c", "degreedist",    "default",      "CondDegreeDist"),
          c("c", "degreedist",    "random",       "CondDegreeDist"), 
          c("c", "indegreedist",  "default",      "CondInDegreeDist"),
          c("c", "indegreedist",  "random",       "CondInDegreeDist"), 
          c("c", "outdegreedist",  "default",      "CondOutDegreeDist"),
          c("c", "outdegreedist",  "random",       "CondOutDegreeDist"), 

#          c("c", "indegrees",     "default",      "CondInDegree"),
#          c("c", "indegrees",     "random",       "CondInDegree"),
#          c("c", "outdegrees",    "default",      "CondOutDegree"),
#          c("c", "outdegrees",    "random",       "CondOutDegree"),
          c("c", "edges",         "default",      "ConstantEdges"),
          c("c", "edges",         "random",       "ConstantEdges"),
          c("c", "hamming",       "default",      "HammingTNT"),
          c("c", "hamming",       "random",       "HammingTNT"),
          c("c", "edges+hamming", "default",      "HammingConstantEdges"),
          c("c", "edges+hamming", "random",       "HammingConstantEdges"),
          c("c", "observed",      "default",      "randomtoggleNonObserved"),
          c("c", "observed",      "random",       "randomtoggleNonObserved")
        )
MHproposals <- data.frame(I(MHproposals[,1]), I(MHproposals[,2]), 
                          I(MHproposals[,3]), I(MHproposals[,4]))  
colnames(MHproposals)<-c("Class","Constraints","Weights","MHP")



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

MHproposal.character <- function(object, arguments, nw, model, ...){
  name<-object
  proposal <- eval(call(paste("InitMHP.", name, sep=""),
                        arguments, nw, model))

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
#   class      :  "c", otherwise execution will halt
#
########################################################################################

MHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c", ...) {

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
  name<-with(MHproposals,MHP[Class==class & Constraints==constraints & Weights==weights])
  if(length(name)>1) stop("Multiple matching proposals in the lookup table.",
                          "This Should Not Be Happening (tm). Unless you have",
                          "been messing with the table, please file a bug report.")
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

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, model=NULL,weights=NULL,class="c"){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$prop.args
  if(is.null(nw)) nw<-object$network
  if(is.null(weights)) weights<-"default"
  if(is.null(model)) model<-ergm.getmodel(object$formula,nw,...)
  MHproposal(constraints,arguments=arguments,nw=nw,model=model,weights=weights,class=class)
}

