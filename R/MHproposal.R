#  File ergm/R/MHproposal.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
########################################################################################
# The <MHproposal> function initializes and returns an MHproposal object via one of the
# class-specific functions below
########################################################################################

MHproposal<-function(object, ...) UseMethod("MHproposal")


# This could be useful for trapping bugs before they become mysterious segfaults.
MHproposal.NULL<-function(object, ...) stop("NULL passed to MHproposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.")


MHproposal.MHproposal<-function(object,...) return(object)




########################################################################################
# The <MHproposal.character> function initializes the MHproposal object using the
# <InitMHP.> function that corresponds to the name given in 'object'
########################################################################################

MHproposal.character <- function(object, arguments, nw, ...){
  name<-object
  proposal <- eval(call(paste("InitMHP", name, sep="."), arguments, nw))

  proposal$arguments <- arguments

  proposal$arguments$constraints$bd <- ergm.bounddeg(arguments$bd,nw)
  
  class(proposal)<-"MHproposal"
  proposal
}





########################################################################################
# The <MHproposal.formula> function verifies that the given constraints exist and
# are supported in conjuction with the given weights and class by a unique MH proposal;
# if so the MHproposal object is created via <MHproposal.character> using the 
# MHP type found in the look-up table above.
########################################################################################

MHproposal.formula <- function(object, arguments, nw, weights="default", class="c", ...) {
  constraints<-object

  if("constraints" %in% names(arguments)){
    conlist <- arguments$constraints
  }else{
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
    MHqualifying<-with(MHproposals,MHproposals[Class==class & Constraints==constraints & if(is.null(weights) || weights=="default") TRUE else Weights==weights,])

  if(nrow(MHqualifying)<1){
    commonalities<-(MHproposals$Class==class)+(MHproposals$Weights==weights)+(MHproposals$Constraints==constraints)
    stop("The combination of class (",class,"), model constraints (",constraints,"), and proposal weighting (",weights,") is not implemented. ", "Check your arguments for typos. ", if(any(commonalities>=3)) paste("Nearest matching proposals: (",paste(apply(MHproposals[commonalities==3,-5],1,paste, sep="), (",collapse=", "),collapse="), ("),")",sep="",".") else "")
  }

  if(nrow(MHqualifying)==1)
    name<-MHqualifying$MHP
  else
    name<-with(MHqualifying,MHP[which.max(Priority)])

  arguments$constraints<-conlist
  ## Hand it off to the class character method.
  MHproposal.character(name,arguments,nw)
}





########################################################################################
# The <MHproposal.ergm> function creates the MHproposal object via <MHproposal.formula>
# after extracting the appropriate parameters from the given ergm object
########################################################################################

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, weights=NULL,class="c"){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$prop.args
  if(is.null(nw)) nw<-object$network
  if(is.null(weights)) weights<-"default"
  
  MHproposal(constraints,arguments=arguments,nw=nw,weights=weights,class=class)
}

