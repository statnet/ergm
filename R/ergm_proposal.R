#  File R/ergm_proposal.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

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
#' reference measures, see \code{\link{ergmReference}}
#'
#' @param Constraints The constraints used in the model. For the list
#'   of constraints, see \code{\link{ergmConstraint}}. They are
#'   specified as a single string of text, with each contrast prefixed
#'   by either `&` for constraints that the proposal *always* enforces
#'   or `|` for constraints that the proposal *can* enforce if needed.
#'
#' @param Priority On existence of multiple qualifying proposals, specifies the
#' priority (`-1`,`0`,`1`, etc.) of proposals to be used.
#' @param Weights The sampling weights on selecting toggles (random, TNT, etc).
#' @param Proposal The matching proposal from the previous arguments.
#'
#' @note The arguments `Class`, `Reference`, and `Constraints` can
#'   have length greater than 1. If this is the case, the rows added
#'   to the table are a *Cartesian product* of their elements.
#' @keywords internal
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

prune.ergm_conlist <- function(conlist){
  ## Remove constraints implied by other constraints.
  for(ed in rev(seq_along(conlist))){
    for(er in rev(seq_along(conlist))){
      if(er==ed) next
      er.con <- NVL(conlist[[er]]$constrain,character(0))
      er.implies <- NVL(conlist[[er]]$implies,character(0))
      ed.con <- NVL(conlist[[ed]]$constrain,character(0))
      ed.impliedby <- NVL(conlist[[ed]]$impliedby,character(0))
      
      if(any(er.con %in% ed.impliedby) || any(ed.con %in% er.implies) || any(er.implies %in% ed.impliedby)){
        conlist[[ed]]<-NULL
        break
      }
    }
  }

  structure(conlist, class = "ergm_conlist", lhs = attr(conlist, "lhs"))
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
#' \code{\link{ergm}} object.  The \code{\link{formula}} should be of the format documented in the `constraints` argument of [ergm()] and in the [ERGM constraints][ergmConstraint] documentation.
#' @param \dots Further arguments passed to other functions.
#' @return Returns an ergm_proposal object: a list with class `ergm_proposal`
#' containing the following named elements:
#' \item{name}{the C name of the proposal}
#' \item{inputs}{inputs to be passed to C}
#' \item{pkgname}{shared library name where the proposal
#' can be found (usually `"ergm"`)}
#' \item{reference}{the reference distribution}
#' \item{arguments}{list of arguments passed to
#' the `InitErgmProposal` function; in particular,
#' \describe{
#' \item{`constraints`}{list of constraints}
#' \item{uid}{a string generated with the proposal, \UIDalgo; different proposals are, generally, guaranteed to have different strings, but identical proposals are not guaranteed to have the same string}
#' }
#' }
#' @seealso \code{\link{InitErgmProposal}}
#' @keywords models internal
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
#' @template term_options
#' @template reference
#' @export
ergm_proposal.character <- function(object, arguments, nw, ..., reference=ergm_reference(trim_env(~Bernoulli), nw, term.options=term.options, ...), term.options=list()){
  name<-object

  arguments$reference <- reference

  f <- locate_prefixed_function(name, if(is.valued(nw)) "InitWtErgmProposal" else "InitErgmProposal", "Metropolis-Hastings proposal")

  prop.call <-
    if((argnames <- names(formals(eval(f))))[1]=="nw"){
      if(! "..."%in%argnames) stop("New-type InitErgmProposal ", sQuote(format(f)), " must have a ... argument.")
      termCall(f, arguments, nw, term.options, ...)
    }else as.call(list(f, arguments, nw))

  proposal <- eval(prop.call)

  storage.mode(proposal$inputs) <- "double"
  storage.mode(proposal$iinputs) <- "integer"
  proposal$arguments <- arguments
  proposal$arguments$reference <- NULL

  proposal$reference <- reference

  # If package not specified, autodetect.
  if(is.null(proposal$pkgname))  proposal$pkgname <- environmentName(environment(eval(f)))

  # Add the package to the list of those to be loaded.
  ergm.MCMC.packagenames(proposal$pkgname)
  
  class(proposal)<-"ergm_proposal"
  proposal$uid <- .GUID()
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
#                 to dyads; options are "TNT", "StratTNT", "TNT10", "random", "nonobserved" and
#                 "default"; default="default"
#   class      :  the class of the proposal; choices include "c", "f", and "d"
#                 default="c"
#
########################################################################################

ergm_conlist <- function(object, ...) UseMethod("ergm_conlist")
ergm_conlist.ergm_conlist <- function(object, ...) object
ergm_conlist.NULL <- function(object, ...) NULL

ergm_conlist.formula <- function(object, nw, ..., term.options=list()){
  env <- environment(object)
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

    f <- locate_prefixed_function(constraint, "InitErgmConstraint", "Sample space constraint")

    if(names(formals(eval(f)))[1]=="lhs.nw"){
      if(is.call(constraint)){
        conname <- as.character(constraint[[1]])
        init.call<-list(f, lhs.nw=nw)
        init.call<-c(init.call,as.list(constraint)[-1])
      }else{
        conname <- as.character(constraint)
        init.call <- list(f, lhs.nw=nw)
      }
    }else{
      conname <- as.character(if(is.name(constraint)) constraint else constraint[[1]])
      init.call <- termCall(f, constraint, nw, term.options, ..., env=env)
    }

    con <- eval(as.call(init.call), env)
    NVL(con$dependence) <- TRUE
    if(con$dependence && consign < 0) stop("Only dyad-independent costraints can have negative signs.")
    con$sign <- consign
    NVL(con$priority) <- Inf
    NVL(con$constrain) <- conname
    if(!con$dependence && con$priority==Inf){
      con$constrain <- if(con$sign < 0) ".dyads" # If disjunctive, override specific in favour of general.
                       else unique(c(con$constrain,".dyads")) # FIXME: should .dyads go first?
    }
#' @import memoise
    if(is.function(con$free_dyads) && !is.memoised(con$free_dyads)) con$free_dyads <- memoise(con$free_dyads)

    conlist[[length(conlist)+1]] <- con
    names(conlist)[length(conlist)] <- conname
  }

  if (length(object) == 3) {
    lhs <- try(eval_lhs.formula(object), silent = TRUE)
    if (is(lhs, "try-error") || !is.character(lhs)) stop("Constraints formula must be either one-sided or have a string expression as its LHS.")

    attr(conlist, "lhs") <- lhs
  }

  prune.ergm_conlist(conlist)
}

c.ergm_conlist <- function(...){
  o <- NextMethod() %>% prune.ergm_conlist()

  lhss <- lapply(list(...), attr, "lhs") %>% compact()
  if(length(lhss)) attr(o, "lhs") <- ult(lhss)

  o
}

`[.ergm_conlist` <- function(x, ...){
  structure(NextMethod(), class = "ergm_conlist", lhs = attr(x, "lhs"))
}

select_ergm_proposal <- function(conlist, class, ref, name, weights){
  candidates <- ergm_proposal_table()
  candidates <- candidates[candidates$Class==class & candidates$Reference==ref$name & if(is.null(weights) || weights=="default") TRUE else candidates$Weights==weights, , drop=FALSE]
  if(!is.null(name)) candidates <- candidates[candidates$Proposal==name, , drop=FALSE]

  decode_constraints <- function(s){
    # Convert old-style specification to the new-style
    # specification. Note that .dyads is always optional.
    if(nchar(s) && !startsWith(s,"&") && !startsWith(s,"|"))
      s <- strsplit(s, "+", fixed=TRUE)[[1L]] %>% paste0(ifelse(.==".dyads", "|", "&"), ., collapse="")

    # Split on flags, but keep flags.
    s <- strsplit(s, "(?<=.)(?=[&|])", perl=TRUE)[[1L]]

    names <- substr(s, 2, 2147483647L)
    does <- map_lgl(s, startsWith, "&")
    can <- !does

    list(does = names[does],
         can = names[can])
  }

  candidates <- as_tibble(candidates)
  candidates$Constraints <- lapply(candidates$Constraints, decode_constraints)

  # proposals = the proposal table
  # constraints = an ergm_conlist
  score_proposals <- function(proposals, conlist){
    conlist <- conlist[names(conlist)!=".attributes"]
    priorities <- map_dbl(conlist, "priority")
    wanted <- map(conlist, "constrain")
    allwanted <- unlist(wanted)

    add_score <- function(proposal){
      propcon <- proposal$Constraints
      does <- propcon$does
      can <- propcon$can
      knows <- c(does, can)

      if(any(! does%in%allwanted)) return(NULL) # Proposal has an unwanted constraint.

      unmet <- map_lgl(wanted, ~!any(. %in% knows))
      wanted <- wanted[unmet]
      # Penalised proposal score.
      proposal$Unmet <- if(length(wanted)) wanted %>% map_chr(paste0, collapse="/") %>% sQuote %>% paste.and else ""
      proposal$UnmetScore <- sum(priorities[unmet])
      proposal$Score <- proposal$Priority - proposal$UnmetScore
      if(proposal$Score==-Inf) return(NULL)
      proposal
    }
    proposals %>% transpose %>% map(add_score) %>% compact %>% transpose %>% map(simplify_simple,toNA="keep") %>% as_tibble
  }

  qualifying <- score_proposals(candidates, conlist)

  if(nrow(qualifying)<1){
    stop("The combination of class (",class,"), model constraints and hints (",paste.and(sQuote(unique(names(conlist)))),"), reference measure (",deparse(ult(ref$name)),"), proposal weighting (",weights,"), and conjunctions and disjunctions is not implemented. ", "Check your arguments for typos. ")
  }

  proposal <- qualifying[which.max(qualifying$Score),]
  if(proposal$Unmet!="") message("Best valid proposal ", sQuote(proposal$Proposal), " cannot take into account hint(s) ", proposal$Unmet, ".")
  proposal
}

ergm_reference <- function(object, ...) UseMethod("ergm_reference")
ergm_reference.ergm_reference <- function(object, ...) object

ergm_reference.formula <- function(object, nw, ..., term.options=list()) {
    env <- environment(object)

    f <- locate_prefixed_function(object[[2]], "InitErgmReference", "Reference distribution")

    if (names(formals(eval(f)))[1] == "lhs.nw") {
      if (is.call(object[[2]])) {
        ref.call <- list(f, lhs.nw = nw)
        ref.call <- c(ref.call, as.list(object[[2]])[-1])
      }else{
        ref.call <- list(f, lhs.nw = nw)
      }
    }else{
      ref.call <- termCall(f, object[[2]], nw, term.options, ..., env = env)
    }

    structure(eval(as.call(ref.call), env), class = "ergm_reference")
}

#' @describeIn ergm_proposal `object` argument is an ERGM constraint formula; constructs the [`ergm_conlist`] object and hands off to `ergm_proposal.ergm_conlist()`.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being simulated. See the
#' documentation for a similar argument for \code{\link{ergm}} and see
#' [list of implemented constraints][ergmConstraint] for more information.
#' @export
ergm_proposal.formula <- function(object, arguments, nw, hints=trim_env(~sparse), ..., term.options=list()) {
  NVL(hints) <- trim_env(~sparse)

  conlist <- if("constraints" %in% names(arguments))
               prune.ergm_conlist(arguments$constraints)
             else c(ergm_conlist(object, nw, term.options=term.options, ...),
                    ergm_conlist(hints, nw, term.options=term.options, ...))

  ## Hand it off to the class ergm_conlist method.
  ergm_proposal(conlist, arguments, nw, ..., term.options = term.options)
}

#' @describeIn ergm_proposal `object` argument is an ERGM constraint
#'   list; constructs the internal `ergm_reference` object, looks up the
#'   proposal, and hands off to `ergm_proposal.character()`.
#' @param weights Specifies the method used to allocate probabilities
#'   of being proposed to dyads, providing an intermediate method
#'   (between hints and specifying the proposal name directly) for
#'   specifying the proposal; options include "TNT", "StratTNT",
#'   "TNT10", "random", "nonobserved" and "default"; default="default"
#' @param class The class of the proposal; choices include "c", "f",
#'   and "d" default="c".
#' @export
ergm_proposal.ergm_conlist <- function(object, arguments, nw, weights="default", class="c", reference=trim_env(~Bernoulli), ..., term.options=list()) {
  reference <- ergm_reference(reference, nw, term.options=term.options, ...)
  proposal <- select_ergm_proposal(object, class = class, ref = reference, name = attr(object, "lhs"), weights = weights)
  name <- proposal$Proposal
  arguments$constraints <- object
  ## Hand it off to the class character method.
  ergm_proposal(name, arguments, nw, reference = reference, ..., term.options = term.options)
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
#                 "TNT", "StratTNT", "random", "TNT10", or "default"; default="default"
#                 (these options don't agree with the prop.weights of <control.ergm>)
#   class      :  "c", otherwise execution will halt
#
########################################################################################

#' @describeIn ergm_proposal `object` argument is an [`ergm`] fit whose proposals are extracted which is reproduced as best as possible.
#' @export
ergm_proposal.ergm<-function(object,...,constraints=NULL, arguments=NULL, nw=NULL, weights=NULL,class="c", reference=NULL){
  if(is.null(constraints)) constraints<-object$constraints
  if(is.null(arguments)) arguments<-object$control$MCMC.prop.args
  if(is.null(weights)) weights<-object$control$MCMC.prop.weights
  if(is.null(nw)) nw<-object$network
  if(is.null(reference)) reference<-object$reference

  ergm_proposal(constraints,arguments=arguments,nw=nw,weights=weights,class=class,reference=reference, ...)
}

DyadGenType <- list(RandDyadGen=0L, WtRandDyadGen=1L, RLEBDM1DGen=2L, WtRLEBDM1DGen=3L, EdgeListGen=4L, WtEdgeListGen=5L)

#' A helper function to select and construct a dyad generator for C.
#'
#' @param arguments argumements passed to the [`ergm_proposal`].
#' @param nw a [`network`].
#' @param extra_rlebdm an [`rlebdm`] representing any additional constraints.
#' @return A list understood by the C `DyadGen` API.
#' @keywords internal
#' @export
ergm_dyadgen_select <- function(arguments, nw, extra_rlebdm=NULL){
  valued <- is.valued(nw)

  dyadgen <- list(valued=valued)

  r <- as.rlebdm(arguments$constraints)
  if(!is.null(extra_rlebdm)) r <- r & extra_rlebdm

  if(all(r==free_dyads(arguments$constraints$.attributes))){
    dyadgen$type <- if(valued) DyadGenType$WtRandDyadGen else DyadGenType$RandDyadGen
  }else{
    # If the number of selectable edges exceeds the number of runs by
    # some factor, use RLEBDM, otherwise just edgelist.
    #
    # TODO: The exact constant needs to be tuned.
    if(sum(r) > length(r$lengths)*20){
      dyadgen$type <- if(valued) DyadGenType$WtRLEBDM1DGen else DyadGenType$RLEBDM1DGen
      dyadgen$dyads <- to_ergm_Cdouble(r)
    }else{
      dyadgen$type <- if(valued) DyadGenType$WtEdgeListGen else DyadGenType$EdgeListGen
      dyadgen$dyads <- as.integer(to_ergm_Cdouble(as.edgelist(r)))
    }
  }
  dyadgen
}

#' A function to extract a constraint's free_dyads RLEBDM.
#'
#' @param con a list containing constraint information as returned by `InitErgmConstraiont.*()`; if `con$free_dyads` is `NULL` or an [`rlebdm`], returned as is; if a function, the function is called with no arguments and the result returned.
#'
#' @return an [`rlebdm`] or `NULL`.
#' @noRd
free_dyads <- function(con){
  fd <- con$free_dyads
  if(is.null(fd) || is(fd, "rlebdm")) fd
  else if(is.function(fd)) fd()
  else stop("Unsupported free_dyad type; this is probably a programming error.")
}
