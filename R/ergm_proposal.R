#  File R/ergm_proposal.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
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
#' reference measures, see [`ergmReference`]
#'
#' @param Constraints The constraints used in the model. For the list
#'   of constraints, see [`ergmConstraint`]. They are
#'   specified as a single string of text, with each contrast prefixed
#'   by either `&` for constraints that the proposal *always* enforces
#'   or `|` for constraints that the proposal *can* enforce if needed.
#'
#' @param Priority On existence of multiple qualifying proposals, specifies the
#' priority (`-1`,`0`,`1`, etc.) of proposals to be used.
#' @param Weights The sampling weights on selecting toggles (random, TNT, etc).
#' @param Proposal The matching proposal from the previous arguments.
#'
#' @param Package The package in which the proposal is implemented;
#'   it's normally autodetected based on the package to which the
#'   calling function belongs.
#'
#' @note The arguments `Class`, `Reference`, and `Constraints` can
#'   have length greater than 1. If this is the case, the rows added
#'   to the table are a *Cartesian product* of their elements.
#'
#' @details The first time a particular package calls
#'   `ergm_proposal_table()`, it also sets a call-back to remove all
#'   of its proposals from the table should the package be unloaded.
#'
#' @keywords internal
#' @export ergm_proposal_table
ergm_proposal_table <- local({
  proposals <- data.frame(Class = character(0), Reference = character(0),
                     Constraints = character(0), Priority = numeric(0), Weights = character(0),
                     Proposal = character(0), Package = character(0), stringsAsFactors=FALSE)

  unload <- function(pkg_name, ...) proposals <<- proposals[proposals$Package != pkg_name, ]

  function(Class, Reference, Constraints, Priority, Weights, Proposal, Package = NULL) {
    if(!missing(Class)){
      NVL(Package) <- utils::packageName(parent.frame())
      NVL(Package) <- ""

      if(Package != "" && ! Package %in% proposals$Package)
        setHook(packageEvent(Package, "onUnload"), unload)

      stopifnot(length(Class)>=1, length(Reference)>=1, length(Constraints)>=1,
                length(Priority)==1, length(Weights)==1, length(Proposal)==1,
                length(Package)==1)

      newrows <- expand.grid(Class = Class, Reference = Reference, Constraints = Constraints,
                             Priority = Priority, Weights = Weights, Proposal = Proposal,
                             Package = Package, stringsAsFactors=FALSE, KEEP.OUT.ATTRS=FALSE)
      proposals <<- rbind(proposals, newrows)
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

  structure(conlist, class = "ergm_conlist")
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

# If the constraints formula is two-sided, add a term .select(LHS) and remove LHS.
.embed_constraint_lhs <- function(formula){
  if(length(formula) > 2){
    lhs <- try(eval_lhs.formula(formula), silent = TRUE)
    if (is(lhs, "try-error") || !is.character(lhs)) stop("Constraint formula must be either one-sided or have a string expression as its LHS.")
    nonsimp_update.formula(formula, substitute(~. + .select(..), list(..=lhs)))
  }else formula
}

.delete_term <- function(tl, terms) discard(tl, ~any(as.character(.)[1] %in% terms))
.keep_term <- function(tl, terms) keep(tl, ~any(as.character(.)[1] %in% terms))
.delete_constraint <- function(cl, constraints) discard(cl, ~any(.$constrain %in% constraints))
.keep_constraint <- function(cl, constraints) keep(cl, ~any(.$constrain %in% constraints))

#' Functions to initialize the ergm_proposal object
#' 
#' S3 Functions that initialize the Metropolis-Hastings Proposal (ergm_proposal)
#' object using the `InitErgmProposal.*` function that corresponds to the name given in
#' 'object'.  These functions are not generally called directly by the user.
#' See [`ergmProposal`] for general explanation and lists of available
#' Metropolis-Hastings proposal types.
#' 
#' 
#' @aliases ergm_proposal.NULL ergm_proposal.ergm_proposal
#' @param object Either a character, a [`formula`] or an
#' [`ergm`] object.  The [`formula`] should be of the format documented in the `constraints` argument of [ergm()] and in the [ERGM constraints][ergmConstraint] documentation.
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
#' @seealso [`InitErgmProposal`]
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
#' @param nw The network object originally given to [ergm()]
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
      termCall(f, arguments, nw, term.options, reference=reference, ...)
    }else as.call(list(f, arguments, nw))

  proposal <- eval(prop.call)
  if(is.null(proposal) || is.character(proposal)) return(proposal)

  storage.mode(proposal$inputs) <- "double"
  storage.mode(proposal$iinputs) <- "integer"
  proposal$arguments <- arguments
  proposal$arguments$reference <- NULL

  proposal$reference <- reference

  # If package not specified, autodetect.
  if(is.null(proposal$pkgname))  proposal$pkgname <- environmentName(environment(eval(f)))

  # Add the package to the list of those to be loaded.
  ergm.MCMC.packagenames(proposal$pkgname)
  check_ABI(proposal$pkgname)

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

#' @noRd
ergm_conlist <- function(object, ...) UseMethod("ergm_conlist")

#' @noRd
ergm_conlist.ergm_conlist <- function(object, ...) object

#' @noRd
ergm_conlist.NULL <- function(object, ...) NULL

#' @noRd
ergm_conlist.formula <- function(object, nw, ..., term.options=list())
  object %>% .embed_constraint_lhs() %>% list_rhs.formula() %>%
    ergm_conlist(nw, ..., term.options=term.options)

#' @noRd
ergm_conlist.term_list <- function(object, nw, ..., term.options=list()){
  object <-
    if(is(object, "AsIs")) structure(object, class = class(object)[class(object) != "AsIs"])
    else c(object, list(call(".attributes")))
  consigns <- sign(object)
  conenvs <- envir(object)

  conlist <- structure(list(), class = "ergm_conlist")
  for(i in seq_along(object)){
    constraint <- object[[i]]
    consign <- consigns[[i]]
    conenv <- conenvs[[i]]

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
      init.call <- termCall(f, constraint, nw, term.options, ..., env=conenv)
    }

    con <- eval(as.call(init.call), conenv)
    if(is.null(con)) next
    else if (is(con, "ergm_conlist")) {
      conlist <- c(conlist, con)
      next
    }

    NVL(con$priority) <- Inf # Default priority
    if(con$priority < Inf) con$dependence <- FALSE # Hints do not induce dependence in the sample space.
    NVL(con$dependence) <- TRUE # Default dependence
    if(con$dependence && consign < 0) stop("Only dyad-independent constraints can have negative signs.")
    con$sign <- consign
    NVL(con$constrain) <- conname
    conname <- con$constrain[1] # If constraint provides a name, use it.
    if(!con$dependence && con$priority==Inf){
      con$constrain <- if(con$sign < 0) ".dyads" # If disjunctive, override specific in favour of general.
                       else unique(c(con$constrain,".dyads")) # FIXME: should .dyads go first?
    }
#' @import memoise
    if(is.function(con$free_dyads) && !is.memoised(con$free_dyads)) con$free_dyads <- memoise(con$free_dyads)

    conlist[[length(conlist)+1]] <- con
    names(conlist)[length(conlist)] <- conname
  }

  prune.ergm_conlist(conlist)
}

#' @noRd
c.ergm_conlist <- function(...) NextMethod() %>% prune.ergm_conlist()

#' @noRd
`[.ergm_conlist` <- function(x, ...){
  structure(NextMethod(), class = "ergm_conlist")
}

select_ergm_proposals <- function(conlist, class, ref, weights){
  # Extract directly selected proposal, if given, check that it's unique, and discard its constraint and other placeholders.
  name <- conlist %>% .keep_constraint(".select") %>% map_chr("proposal") %>% unique()
  if (length(name) > 1) stop("Error in direct proposal selection: multiple distinct proposals selected: ", paste.and(sQuote(name)), ".", call. = FALSE)
  conlist <- conlist %>% .delete_constraint(c(".", ".select"))

  # Initial narrowing down of the proposal table.
  candidates <- ergm_proposal_table()
  candidates <- candidates[candidates$Class==class & candidates$Reference==ref$name & if(is.null(weights) || weights=="default") TRUE else candidates$Weights==weights, , drop=FALSE]
  if(length(name)) candidates <- candidates[candidates$Proposal==name, , drop=FALSE]

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

  qualifying[order(qualifying$Score, decreasing=TRUE),]
}

ergm_reference <- function(object, ...) UseMethod("ergm_reference")
ergm_reference.ergm_reference <- function(object, ...) object

#' @noRd
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
#' documentation for a similar argument for [ergm()] and see
#' [`ergmConstraint`] for more information.
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

#' @describeIn ergm_proposal `object` argument is a [`term_list`];
#'   same implementation as the `formula` method.
#' @export
ergm_proposal.term_list <- ergm_proposal.formula

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
  proposals <- select_ergm_proposals(object, class = class, ref = reference, weights = weights)

  for(i in seq_len(nrow(proposals))){
    proposal <- proposals[i,]
    name <- proposal$Proposal
    arguments$constraints <- object
    ## Hand it off to the class character method.
    proposal <- ergm_proposal(name, arguments, nw, reference = reference, weights = weights, class = class, ..., term.options = term.options)

    ## Keep trying until some proposal function accepts.
    if(!is.null(proposal) && !is.character(proposal)) break
    if(is.character(proposal)) message(proposal," Falling back.")
  }

  if(proposals[i,]$Unmet!="") message("Best valid proposal ", sQuote(proposals[i,]$Proposal), " cannot take into account hint(s) ", proposals[i,]$Unmet, ".")
  proposal
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
