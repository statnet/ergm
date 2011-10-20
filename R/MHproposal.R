#=======================================================================================
# This file contains the following 6 files for creating MHproposal objects
#          <MHproposal>                <MHproposal.character>
#          <MHproposal.NULL>           <MHproposal.formula>
#          <MHproposal.MHproposal>     <MHproposal.ergm>
#=======================================================================================

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
#   class      :  the class of the proposal; choices include "c", "f", and "d"
#                 default="c"
#
########################################################################################

MHproposal.formula <- function(object, arguments, nw, model, weights="default", class="c", reference="Bernoulli", response=NULL, ...) {
  constraints<-object
  reference<-match.arg(reference,unique(MHproposals$Reference))
  
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
    conlist <- try(eval(as.call(init.call), environment(object)))
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
  MHqualifying<-with(MHproposals,MHproposals[Class==class & Constraints==constraints & Reference==reference & if(is.null(weights) || weights=="default") TRUE else Weights==weights,])

  if(nrow(MHqualifying)<1){
    commonalities<-(MHproposals$Class==class)+(MHproposals$Weights==weights)+(MHproposals$Reference==reference)+(MHproposals$Constraints==constraints)
    stop("The combination of class (",class,"), model constraints (",constraints,"), reference measure (",reference,"), and proposal weighting (",weights,") is not implemented. ", "Check your arguments for typos. ", if(any(commonalities>=3)) paste("Nearest matching proposals: (",paste(apply(MHproposals[commonalities==3,-5],1,paste, sep="), (",collapse=", "),collapse="), ("),")",sep="",".") else "")
  }

  if(nrow(MHqualifying)==1)
    name<-MHqualifying$MHP
  else
    name<-with(MHqualifying,MHP[which.max(Priority)])
  
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
#                 'silent' and 'initialfit'
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

