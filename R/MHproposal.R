#  File R/MHproposal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

#=======================================================================================
# This file contains the following 6 files for creating MHproposal objects
#          <MHproposal>                <MHproposal.character>
#          <MHproposal.NULL>           <MHproposal.formula>
#          <MHproposal.MHproposal>     <MHproposal.ergm>
#=======================================================================================

# Set up the table mapping constraints, references, etc. to
# MHproposals.  For the moment, there is no way to delete rows, but
# one can always add a row with identical elements but higher
# priority.
ergm.MHP.table <- local({
  MHPs <- data.frame(Class = character(0), Reference = character(0),
                     Constraints = character(0), Priority = numeric(0), Weights = character(0),
                     MHP = character(0), stringsAsFactors=FALSE)
  function(Class, Reference, Constraints, Priority, Weights, MHP) {
    if(!missing(Class))
      MHPs <<- rbind(MHPs,
                     data.frame(Class = Class, Reference = Reference,
                     Constraints = Constraints, Priority = Priority, Weights = Weights,
                     MHP = MHP, stringsAsFactors=FALSE))
    else
      MHPs
  }
})

# Set up the list mapping constraint implications.
ergm.ConstraintImplications <- local({
  implications <- list()
  
  function(implier, implies) {
    if(!missing(implier)){
      new.implies <- c(implications[[implier]], implies)
      implications[implier] <<- if(length(new.implies)) list(sort(unique(new.implies))) else list(NULL)
    }
    else
      implications
  }
})


prune.conlist <- function(conlist){
  ## Remove constraints implied by other constraints.
  for(constr in names(conlist))
    for(impl in ergm.ConstraintImplications()[[constr]])
      conlist[[impl]]<-NULL

  conlist
}


########################################################################################
# The <MHproposal> function initializes and returns an MHproposal object via one of the
# class-specific functions below
#
# --PARAMETERS--
#   (see the class-specific function headers)
#
# --RETURNED--
#   proposal: an MHproposal object as a list containing the following:
#     name   : the C name of the proposal
#     inputs : NULL (I think - the only non-null value returned by the InitMH
#              is for <nobetweengroupties>, but this isn't included in the 
#              look-up table
#     package: shared library name where the proposal can be found (usually "ergm")
#     arguments: list of arguments passed to the InitMHP function; in particular,
#       constraints: list of constraints
#       constraints$bd: the list of parameters to bound degree in the fitting process
#              and returned by <ergm.bounddeg>
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
#
########################################################################################

MHproposal.character <- function(object, arguments, nw, ..., response=NULL, reference=reference){
  name<-object

  arguments$reference <- reference

  f <- locate.InitFunction(name, if(is.null(response)) "InitMHP" else "InitWtMHP", "Metropolis-Hastings proposal")

  proposal <- {
    if(is.null(response))
      eval(as.call(list(f, arguments, nw)))
    else
      eval(as.call(list(f, arguments, nw, response)))
  }

  proposal$arguments <- arguments
  proposal$arguments$reference <- NULL

  proposal$reference <- reference

  proposal$arguments$constraints$bd <- ergm.bounddeg(arguments$constraints$bd,nw)
  # If package not specified, autodetect.
  if(is.null(proposal$pkgname))  proposal$pkgname <- environmentName(environment(f))

  # Add the package to the list of those to be loaded.
  ergm.MCMC.packagenames(proposal$pkgname)
  
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
#   constraints:  the constraints as a one sided formula '~ term(s)'
#   weights    :  specifies the method used to allocate probabilities of being proposed
#                 to dyads; options are "TNT", "TNT10", "random", "nonobserved" and
#                 "default"; default="default"
#   class      :  the class of the proposal; choices include "c", "f", and "d"
#                 default="c"
#
########################################################################################

mk.conlist <- function(object, nw){
  if(is.null(object)) return(NULL)
  ## Construct a list of constraints and arguments from the formula.
  conlist<-list()
  constraints<-as.list(attr(terms(object,allowDotAsName=TRUE),"variables"))[-1]
  for(constraint in constraints){
    ## The . in the default formula means no constrains.
    ## There may be other constraints in the formula, however.
    if(constraint==".") next

    f <- locate.InitFunction(constraint, "InitConstraint", "Sample space constraint")
    
    if(is.call(constraint)){
      init.call<-list(f, conlist=conlist, lhs.nw=nw)
      init.call<-c(init.call,as.list(constraint)[-1])
    }else{
      init.call <- list(f, conlist=conlist, lhs.nw=nw)
    }
    conlist <- eval(as.call(init.call), environment(object))
  }
  conlist <- prune.conlist(conlist)
  class(conlist) <- "conlist"
  conlist
}

MHproposal.formula <- function(object, arguments, nw, weights="default", class="c", reference=~Bernoulli, response=NULL, ...) {
  reference <- reference

  f <- locate.InitFunction(reference[[2]], "InitReference", "Reference distribution") 
  
  if(is.call(reference[[2]])){
    ref.call <- list(f, lhs.nw=nw)
    ref.call <- c(ref.call,as.list(reference[[2]])[-1])
  }else{
    ref.call <- list(f, lhs.nw=nw)
  }
  reference <- eval(as.call(ref.call),environment(reference))

  if("constraints" %in% names(arguments)){
    conlist <- prune.conlist(arguments$constraints)
    class(conlist) <- "conlist"
  }else{
    conlist <- mk.conlist(object, nw)
  }
  
  ## Convert vector of constraints to a "standard form".
  if(is.null(names(conlist))) {
    constraints <- ""
  } else {
    constraints <- paste(sort(tolower(names(conlist))),collapse="+")
  }
  
  MHqualifying<-with(ergm.MHP.table(),ergm.MHP.table()[Class==class & Constraints==constraints & Reference==reference$name & if(is.null(weights) || weights=="default") TRUE else Weights==weights,])

  if(nrow(MHqualifying)<1){
    commonalities<-(ergm.MHP.table()$Class==class)+(ergm.MHP.table()$Weights==weights)+(ergm.MHP.table()$Reference==reference)+(ergm.MHP.table()$Constraints==constraints)
    stop("The combination of class (",class,"), model constraints (",constraints,"), reference measure (",reference,"), and proposal weighting (",weights,") is not implemented. ", "Check your arguments for typos. ", if(any(commonalities>=3)) paste("Nearest matching proposals: (",paste(apply(ergm.MHP.table()[commonalities==3,-5],1,paste, sep="), (",collapse=", "),collapse="), ("),")",sep="",".") else "")
  }

  if(nrow(MHqualifying)==1){
    name<-MHqualifying$MHP
  }else{
    name<-with(MHqualifying,MHP[which.max(Priority)])
  }
    
  arguments$constraints<-conlist
  ## Hand it off to the class character method.
  MHproposal.character(name, arguments, nw, response=response, reference=reference)
}





########################################################################################
# The <MHproposal.ergm> function creates the MHproposal object via <MHproposal.formula>
# after extracting the appropriate parameters from the given ergm object
#
# --PARAMETERS--
#   object     :  an ergm object
#   constraints:  the constraints as a one sided formula '~ term(s)';
#                 default=object$constraints
#   arguments  :  a list of parameters used by the <Init.MHP> routines  possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  a network object; default=object.network
#   weights    :  the proposal weights component of <control.ergm> as either
#                 "TNT", "random", "TNT10", or "default"; default="default"
#                 (these options don't agree with the prop.weights of <control.ergm>)
#   class      :  "c", otherwise execution will halt
#
########################################################################################

MHproposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, weights=NULL,class="c", reference=NULL, response=NULL){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$control$MCMC.prop.args
  if(is.null(weights)) weights<-object$control$MCMC.prop.weights
  if(is.null(nw)) nw<-object$network
  if(is.null(reference)) response<-object$reference
  if(is.null(response)) response<-object$response
  
  MHproposal(constraints,arguments=arguments,nw=nw,weights=weights,class=class,reference=reference,response=response)
}

