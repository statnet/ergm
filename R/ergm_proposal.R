#  File R/ergm_proposal.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2018 Statnet Commons
#######################################################################

#=======================================================================================
# This file contains the following 6 files for creating ergm_proposal objects
#          <ergm_proposal>                <ergm_proposal.character>
#          <ergm_proposal.NULL>           <ergm_proposal.formula>
#          <ergm_proposal.ergm_proposal>     <ergm_proposal.ergm>
#=======================================================================================

# Set up the table mapping constraints, references, etc. to
# ergm_proposals.  For the moment, there is no way to delete rows, but
# one can always add a row with identical elements but higher
# priority.


#' Table mapping reference,constraints, etc. to ERGM Metropolis-Hastings proposals
#' 
#' This is a low-level function not intended to be called directly by
#' end users. For information on Metropolis-Hastings proposal methods,
#' \link{ergm-proposals}. This function sets up the table mapping
#' constraints, references, etc. to `ergm_proposals`. Calling it with
#' arguments adds a new row to this table. Calling it without
#' arguments returns the table so far.
#' 
#' 
#' @param Class default to "c"
#' @param Reference The reference measure used in the model. For the list of
#' reference measures, see \code{\link{ergm-references}}
#' @param Constraints The constraints used in the model. For the list of
#' constraints, see \code{\link{ergm-constraints}}
#' @param Priority On existence of multiple qualifying proposals, specifies the
#' priority (`-1`,`0`,`1`, etc.) of proposals to be used.
#' @param Weights The sampling weights on selecting toggles (random, TNT, etc).
#' @param Proposal The matching proposal from the previous arguments.
#'
#' @note The arguments `Class`, `Reference`, and `Constraints` can
#'   have length greater than 1. If this is the case, the rows added
#'   to the table are a *Cartesian product* of their elements.
#' @export ergm_proposal_table
ergm_proposal_table <- local({
  proposals <- data.frame(Class = character(0), Reference = character(0),
                     Constraints = character(0), Priority = numeric(0), Weights = character(0),
                     Proposal = character(0), stringsAsFactors=FALSE)
  function(Class, Reference, Constraints, Priority, Weights, Proposal) {
    if(!missing(Class)){
      stopifnot(length(Class)>=1,length(Reference)>=1,length(Constraints)>=1,
                length(Priority)==1,length(Weights)==1,length(Proposal)==1)
      newrows <- expand.grid(Class = Class, Reference = Reference,
                             Constraints = Constraints, Priority = Priority, Weights = Weights,
                             Proposal = Proposal, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)
      proposals <<- rbind(proposals,
                     newrows)
    }else proposals
  }
})

#' @describeIn ergm-deprecated Deprecated name for [ergm_proposal_table()].
#' @export ergm.MHP.table
ergm.MHP.table <- function(...){
  .dep_once("ergm_proposal_table()")
  ergm_proposal_table(...)
}

prune.ergm_conlist <- function(conlist){
  ## Remove constraints implied by other constraints.
  for(ed in rev(seq_along(conlist))){
    for(er in rev(seq_along(conlist))){
      if(er==ed) next
      er.con <- NVL(conlist[[er]]$constrain,character(0))
      er.implies <- NVL(conlist[[er]]$implies,character(0))
      ed.con <- NVL(conlist[[ed]]$constrain,character(0))
      ed.impliedby <- NVL(conlist[[ed]]$impliedby,character(0))
      
      if(any(er.con %in% ed.impliedby) || any(ed.con %in% er.implies) || any(er.implies %in% ed.impliedby))
        conlist[[ed]]<-NULL
    }
  }
  conlist
}


########################################################################################
# The <ergm_proposal> function initializes and returns an ergm_proposal object via one of the
# class-specific functions below
#
# --PARAMETERS--
#   (see the class-specific function headers)
#
# --RETURNED--
#   proposal: an ergm_proposal object as a list containing the following:
#     name   : the C name of the proposal
#     inputs : NULL (I think - the only non-null value returned by the InitErgmProposal
#              is for <nobetweengroupties>, but this isn't included in the 
#              look-up table
#     package: shared library name where the proposal can be found (usually "ergm")
#     arguments: list of arguments passed to the InitErgmProposal function; in particular,
#       constraints: list of constraints
#       constraints$bd: the list of parameters to bound degree in the fitting process
#              and returned by <ergm.bounddeg>
#
########################################################################################



#' Functions to initialize the ergm_proposal object
#' 
#' S3 Functions that initialize the Metropolis-Hastings Proposal (ergm_proposal)
#' object using the `InitErgmProposal.*` function that corresponds to the name given in
#' 'object'.  These functions are not generally called directly by the user.
#' See \link{ergm-proposals} for general explanation and lists of available
#' Metropolis-Hastings proposal types.
#' 
#' 
#' @aliases ergm_proposal.NULL ergm_proposal.ergm_proposal
#' @param object Either a character, a \code{\link{formula}} or an
#' \code{\link{ergm}} object.  The \code{\link{formula}} should be of the form
#' \code{y ~ <model terms>}, where \code{y} is a network object or a matrix
#' that can be coerced to a \code{\link[network]{network}} object.
#' @param \dots Further arguments passed to other functions.
#' @return Returns an ergm_proposal object: a list with class `ergm_proposal`
#' containing the following named elements:
#' \item{ name}{the C name of the proposal}
#' \item{inputs}{inputs to be passed to C}
#' \item{package}{shared library name where the proposal
#' can be found (usually `"ergm"`)}
#' \item{arguments}{list of arguments passed to
#' the `InitErgmProposal` function; in particular,
#' \describe{
#' \item{`constraints`}{list of constraints}
#' }
#' }
#' @seealso \code{\link{InitErgmProposal}}
#' @keywords models
#' @export
ergm_proposal<-function(object, ...) UseMethod("ergm_proposal")


#' @noRd
#' @export
# This could be useful for trapping bugs before they become mysterious segfaults.
ergm_proposal.NULL<-function(object, ...) stop("NULL passed to ergm_proposal. This may be due to passing an ergm object from an earlier version. If this is the case, please refit it with the latest version, and try again. If this is not the case, this may be a bug, so please file a bug report.")


#' @noRd
#' @export
ergm_proposal.ergm_proposal<-function(object,...) return(object)




########################################################################################
# The <ergm_proposal.character> function initializes the ergm_proposal object using the
# <InitErgmProposal.> function that corresponds to the name given in 'object'
#
# --PARAMETERS--
#   object     :  the name of the proposal, one found in the look-up table
#   arguments  :  a list of parameters used by the <InitErgmProposal> routines possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  the network object orginally given to <ergm> via 'formula'
#
########################################################################################

#' @describeIn ergm_proposal `object` argument is a character string
#'   giving the \R name of the proposal.
#' @param nw The network object originally given to \code{\link{ergm}}
#'   via 'formula'
#' @param arguments A list of parameters used by the InitErgmProposal routines
#' @template response
#' @param reference One-sided formula whose RHS gives the reference
#'   measure to be used. (Defaults to \code{~Bernoulli}.)
#' @export
ergm_proposal.character <- function(object, arguments, nw, ..., response=NULL, reference=reference){
  name<-object

  arguments$reference <- reference

  # TODO: Remove this around the end of 2018.
  f <- try(
    locate.InitFunction(name, NVL2(response, "InitWtMHP", "InitMHP"), "Metropolis-Hastings proposal"),
    silent=TRUE
  )
  if(!is(f, "try-error")){
    .dep_once(msg="InitWtMHP and InitMHP convention has been replaced by InitWtErgm and InitErgm, respectively.")
  }else f <- locate.InitFunction(name, NVL2(response, "InitWtErgmProposal", "InitErgmProposal"), "Metropolis-Hastings proposal")

  proposal <- NVL3(response,
                   eval(as.call(list(f, arguments, nw, .))),
                   eval(as.call(list(f, arguments, nw))))

  proposal$arguments <- arguments
  proposal$arguments$reference <- NULL

  proposal$reference <- reference

  proposal$arguments$constraints$bd <- ergm.bounddeg(arguments$constraints$bd,nw)
  # If package not specified, autodetect.
  if(is.null(proposal$pkgname))  proposal$pkgname <- environmentName(environment(eval(f)))

  # Add the package to the list of those to be loaded.
  ergm.MCMC.packagenames(proposal$pkgname)
  
  class(proposal)<-"ergm_proposal"
  proposal
}





########################################################################################
# The <ergm_proposal.formula> function verifies that the given constraints exist and
# are supported in conjuction with the given weights and class by a unique proposal;
# if so the ergm_proposal object is created via <ergm_proposal.character> using the 
# proposal type found in the look-up table above.
#
# --PARAMETERS--
#   object     :  a one-sided formula of constraint terms ( ~ term(s))
#   arguments  :  a list of parameters used by the <InitErgmProposal> routines  possibly including
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

ergm_conlist <- function(object, nw){
  if(is.null(object)) return(NULL)
  ## Construct a list of constraints and arguments from the formula.
  conlist<-list()
  constraints<-list_rhs.formula(object)
  consigns <- c(attr(constraints, "sign"), +1)
  constraints<-c(constraints, list(call(".attributes")))
  for(i in seq_along(constraints)){
    constraint <- constraints[[i]]
    consign <- consigns[[i]]
    
    ## The . in the default formula means no constraints.
    ## There may be other constraints in the formula, however.
    if(constraint==".") next

    # TODO: Remove this around the end of 2018.
    f <- try(
      f <- locate.InitFunction(constraint, "InitConstraint", "Sample space constraint"),
      silent=TRUE
    )
    if(!is(f, "try-error")){
      .dep_once(msg="InitConstraint convention has been replaced by InitErgmConstraint.")
    }else f <- locate.InitFunction(constraint, "InitErgmConstraint", "Sample space constraint")
    
    if(is.call(constraint)){
      conname <- as.character(constraint[[1]])
      init.call<-list(f, lhs.nw=nw)
      init.call<-c(init.call,as.list(constraint)[-1])
    }else{
      conname <- as.character(constraint)
      init.call <- list(f, lhs.nw=nw)
    }
    con <- eval(as.call(init.call), environment(object))
    NVL(con$dependence) <- TRUE
    if(con$dependence && consign < 0) stop("Only dyad-independent costraints can have negative signs.")
    con$sign <- consign

    if(is.null(con$constrain)) con$constrain <- conname
    
    conlist[[length(conlist)+1]] <- con
    names(conlist)[length(conlist)] <- conname
  }
  
  conlist <- prune.ergm_conlist(conlist)
  
  class(conlist) <- "ergm_conlist"
  conlist
}

#' @describeIn ergm_proposal `object` argument is an ERGM constraint formula.
#' @param weights Specifies the method used to allocate probabilities of being
#' proposed to dyads; options are "TNT", "TNT10", "random", "nonobserved" and
#' "default"; default="default"
#' @param class The class of the proposal; choices include "c", "f", and "d"
#' default="c".
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being simulated. See the
#' documentation for a similar argument for \code{\link{ergm}} and see
#' [list of implemented constraints][ergm-constraints] for more information.
#' @export
ergm_proposal.formula <- function(object, arguments, nw, weights="default", class="c", reference=~Bernoulli, response=NULL, ...) {
  reference <- reference

  # TODO: Remove this around the end of 2018.
  f <- try(
    locate.InitFunction(reference[[2]], "InitReference", "Reference distribution"),
    silent=TRUE
  )
  if(!is(f, "try-error")){
    .dep_once(msg="InitReference convention has been replaced by InitErgmReference.")
  }else f <- locate.InitFunction(reference[[2]], "InitErgmReference", "Reference distribution") 
  
  if(is.call(reference[[2]])){
    ref.call <- list(f, lhs.nw=nw)
    ref.call <- c(ref.call,as.list(reference[[2]])[-1])
  }else{
    ref.call <- list(f, lhs.nw=nw)
  }
  reference <- eval(as.call(ref.call),environment(reference))

  # TODO: Remove this around the end of 2018.
  if(is.null(reference$init_methods)){
    .dep_once(msg="Initial methods are now specified in the InitErgmReference.* function. Defaulting to 'CD' and 'zeros'.")
    reference$init_methods <- c("CD", "zeros")
  }

  if(length(object)==3){
    lhs <- object[[2]]
    if(is.character(lhs)){
      name <- object[[2]]
    }else{
      name <- try(eval(lhs, envir=environment(object)), silent=TRUE)
      if(is(name, "try-error") || !is.character(name)) stop("Constraints formula must be either one-sided or have a string expression as its LHS.")
    }
    object <- object[-2]
  }else name <- NULL
  
  if("constraints" %in% names(arguments)){
    conlist <- prune.ergm_conlist(arguments$constraints)
    class(conlist) <- "ergm_conlist"
  }else{
    conlist <- ergm_conlist(object, nw)
  }

  if(is.null(name)){ # Unless specified, autodetect.
    ## Convert vector of constraints to a "standard form".

    # Try the specific constraint combination.
    constraints.specific <- tolower(unlist(lapply(conlist, `[[`, "constrain")))
    constraints.specific <- paste(sort(unique(constraints.specific)),collapse="+")

    qualifying.specific <-
      if(all(sapply(conlist, `[[`, "sign")==+1)){ # If all constraints are conjunctive...
        with(ergm_proposal_table(),ergm_proposal_table()[Class==class & Constraints==constraints.specific & Reference==reference$name & if(is.null(weights) || weights=="default") TRUE else Weights==weights,])
      }

    # Try the general dyad-independent constraint combination.
    constraints.general <- tolower(unlist(ifelse(sapply(conlist,`[[`,"dependence"),lapply(conlist, `[[`, "constrain"), ".dyads")))
    constraints.general <- paste(sort(unique(constraints.general)),collapse="+")
    qualifying.general <- with(ergm_proposal_table(),ergm_proposal_table()[Class==class & Constraints==constraints.general & Reference==reference$name & if(is.null(weights) || weights=="default") TRUE else Weights==weights,])

    qualifying <- rbind(qualifying.general, qualifying.specific)
    
    if(nrow(qualifying)<1){
      commonalities<-(ergm_proposal_table()$Class==class)+(ergm_proposal_table()$Weights==weights)+(ergm_proposal_table()$Reference==reference)+(ergm_proposal_table()$Constraints==constraints.specific)
      stop("The combination of class (",class,"), model constraints (",constraints.specific,"), reference measure (",reference,"), proposal weighting (",weights,"), and conjunctions and disjunctions is not implemented. ", "Check your arguments for typos. ", if(any(commonalities>=3)) paste("Nearest matching proposals: (",paste(apply(ergm_proposal_table()[commonalities==3,-5],1,paste, sep="), (",collapse=", "),collapse="), ("),")",sep="",".") else "")
    }
    
    if(nrow(qualifying)==1){
      name<-qualifying$Proposal
    }else{
      name<-with(qualifying,Proposal[which.max(Priority)])
    }
  }
  
  arguments$constraints<-conlist
  ## Hand it off to the class character method.
  ergm_proposal.character(name, arguments, nw, response=response, reference=reference)
}





########################################################################################
# The <ergm_proposal.ergm> function creates the ergm_proposal object via <ergm_proposal.formula>
# after extracting the appropriate parameters from the given ergm object
#
# --PARAMETERS--
#   object     :  an ergm object
#   constraints:  the constraints as a one sided formula '~ term(s)';
#                 default=object$constraints
#   arguments  :  a list of parameters used by the <InitErgmProposal> routines  possibly including
#                  bd: a list of parameters used to bound degree via <ergm.bounddeg>
#   nw         :  a network object; default=object.network
#   weights    :  the proposal weights component of <control.ergm> as either
#                 "TNT", "random", "TNT10", or "default"; default="default"
#                 (these options don't agree with the prop.weights of <control.ergm>)
#   class      :  "c", otherwise execution will halt
#
########################################################################################

#' @describeIn ergm_proposal `object` argument is an [`ergm`] fit whose proposals are extracted which is reproduced as best as possible.
#' @export
ergm_proposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, weights=NULL,class="c", reference=NULL, response=NULL){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$control$MCMC.prop.args
  if(is.null(weights)) weights<-object$control$MCMC.prop.weights
  if(is.null(nw)) nw<-object$network
  if(is.null(reference)) response<-object$reference
  if(is.null(response)) response<-object$response
  
  ergm_proposal(constraints,arguments=arguments,nw=nw,weights=weights,class=class,reference=reference,response=response)
}

#' @describeIn ergm-deprecated Deprecated name of [ergm_proposal()].
#' @export MHproposal
MHproposal <- function(...){
  .dep_once("ergm_proposal()")
  ergm_proposal(...)
}

#' @describeIn ergm-deprecated Deprecated name of [ergm_proposal()].
#' @export MHproposal.character
MHproposal.character <- function(...){
  .dep_once("ergm_proposal()")
  ergm_proposal.character(...)
}

#' @describeIn ergm-deprecated Deprecated name of [ergm_proposal()].
#' @export MHproposal.ergm
MHproposal.ergm <- function(...){
  .dep_once("ergm_proposal()")
  ergm_proposal.ergm(...)
}

#' @describeIn ergm-deprecated Deprecated name of [ergm_proposal()].
#' @export MHproposal.formula
MHproposal.formula <- function(...){
  .dep_once("ergm_proposal()")
  ergm_proposal.formula(...)
}

#' @describeIn ergm-deprecated Deprecated: specify in the `InitErgmReference.*` implementation.
#' @export ergm.init.methods
ergm.init.methods <- function(...){
  .dep_once(msg="Function ergm.init.methods() has been deprecated in favor of specifying init_methods in InitErgmReference.*() functions, and has no effect.")
}

#' @describeIn ergm-deprecated Deprecated: specify in the `InitErgmConstraint.*` implementation.
#' @export ergm.ConstraintImplications
ergm.ConstraintImplications <- function(...){
  .dep_once(msg="Function ergm.ConstraintImplications() has been deprecated in favor of specifying the implications in the InitErgmConstraint.*() functions, and has no effect.")

}
