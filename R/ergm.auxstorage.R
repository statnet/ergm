#  File R/ergm.auxstorage.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
ergm.auxstorage <- function(model, nw,..., extra.aux=list(), term.options=list()){

  aux_list_list <- function(terms, extra=NULL) {
    # As formulas
    aux.reqs <- c(lapply(terms, "[[", "auxiliaries"), extra)
    
    # A nested list: outer list are the model terms requesting
    # auxiliaries and the inner list is the outputs from InitErgmTerm
    # calls of the auxiliaries.
    lapply(aux.reqs, function(aux.req){
      if(is.null(aux.req)) list()
      else{
        if(is(aux.req, "formula")) aux.req <- list_rhs.formula(aux.req)
        lapply(seq_along(aux.req), function(i){
          call.ErgmTerm(aux.req[[i]], attr(aux.req, "env")[[i]], nw, term.options=term.options, ...)
        })
      }
    })
  }

  aux.outlists <- aux_list_list(model$terms, extra.aux)

  # Remove duplicated auxiliaries.
  uniq.aux.outlists <- unique_aux_terms(unlist(aux.outlists, recursive=FALSE))
  prev <- NULL
  aux.aux.outlists <- list()

  # Until we reach a fixed point (which we should, unless there is a circular dependency).
  #
  # TODO: Check for circular dependencies.
  while(!identical(uniq.aux.outlists,prev, ignore.environment=TRUE)){
    prev <- uniq.aux.outlists
    aux.aux.outlists <- aux_list_list(uniq.aux.outlists)
    uniq.aux.outlists <- unique_aux_terms(c(uniq.aux.outlists, unlist(aux.aux.outlists, recursive=FALSE)))
  }

  # uniq.aux.outlists is now a list of unique initialized auxiliaries
  # and auxiliaries' auxiliaries, such that depended-on auxiliaries
  # are always after the dependent auxiliaries.
  #
  # aux.aux.outlists is now a nested list in the same form as aux.outlists, but for auxiliaries.

  # Append the auxiliary terms to the model.
  n.stat.terms <- length(model$terms)
  for(i in seq_along(uniq.aux.outlists)){
    aux.outlist <- uniq.aux.outlists[[i]]
    model <- updatemodel.ErgmTerm(model, aux.outlist)
    attr(model$terms[[n.stat.terms+i]],"aux.slots")[1L] <- i-1L # The storage slot belonging to this auxiliary.
  }

  # Which term is requiring which auxiliary slot? (+1)
  aux.slots <- lapply(aux.outlists, match_aux_terms, uniq.aux.outlists)
  slots.extra.aux <- rep(list(integer(0)), length(extra.aux))
  
  for(i in seq_along(aux.outlists)){
    if(length(aux.outlists[[i]])){
      if(i<=n.stat.terms) # If it's a model term.
        attr(model$terms[[i]],"aux.slots")[seq_len(length(aux.outlists[[i]]))] <- aux.slots[[i]]-1L
      else # If it's some other entity requesting auxiliaries.
        slots.extra.aux[[i-n.stat.terms]] <- aux.slots[[i]]-1L
    }
  }
  names(slots.extra.aux) <- names(extra.aux)

  # Which auxiliary is requiring which auxiliary slot? (+1)
  aux.aux.slots <- lapply(aux.aux.outlists, match_aux_terms, uniq.aux.outlists)
  
  for(i in seq_along(aux.aux.outlists)){
    if(length(aux.aux.outlists[[i]])){
      # 1st slot is the auxiliary's own slot, so its auxiliaries get put into subsequent slots.
      attr(model$terms[[n.stat.terms+i]], "aux.slots")[1L+seq_len(length(aux.aux.outlists[[i]]))] <- aux.aux.slots[[i]]-1L
    }
  }

  # Check that terms and auxiliaries are properly positioned.
  assert_aux_dependencies(model$terms)

  # For those terms exporting dependence = NA, set it based on their auxiliaries.
  model$terms <- set_aux_dependence(model$terms)

  model$slots.extra.aux <- slots.extra.aux
  model
}

unique_aux_terms <- function(terms){
  # Known issue: unique() and match() don't necessarily have the same notion of equality. This can cause problems. Hopefully, assert_aux_dependencies() can catch them early.
  IGNORE <- "call"
  terms.clean <- lapply(terms, function(term) term[! names(term)%in%IGNORE])
  terms[!duplicated(terms.clean, fromLast=TRUE)]
}

match_aux_terms <- function(x, table){
  # Known issue: unique() and match() don't necessarily have the same notion of equality. This can cause problems. Hopefully, assert_aux_dependencies() can catch them early.
  IGNORE <- "call"
  x.clean <- lapply(x, function(term) term[! names(term)%in%IGNORE])
  table.clean <- lapply(table, function(term) term[! names(term)%in%IGNORE])
  match(x.clean, table.clean)
}

assert_aux_dependencies <- function(terms){
  aux <- (lapply(terms, `[[`, "coef.names") %>% lengths)==0
  aux.slots <- lapply(terms, attr, "aux.slots")

  provided <- ifelse(aux, lapply(aux.slots, `[`, 1), NA)
  requested <- ifelse(aux, lapply(aux.slots, `[`, -1), aux.slots)

  if(anyNA(aux.slots, TRUE)) stop("A requested auxiliary is not provided or is positioned before the requester in the term list. This indicates an implementation bug.")

  for(i in seq_along(aux))
    if(any(! requested[[i]] %in% unlist(provided[-seq_len(i)])))
      stop("A requested auxiliary is not provided or is positioned before the requester in the term list. This indicates an implementation bug.")
}

set_aux_dependence <- function(terms){
  is_aux <- lengths(map(terms, "coef.names")) == 0
  aux_slots <- ifelse(is_aux, map(terms, attr, "aux.slots") %>% map_int(`[`, 1L), NA)

  for(i in rev(seq_along(terms))){
    if(is.na(terms[[i]]$dependence)){ # A term is dependent if any of its auxiliaries is dependent.
      terms[[i]]$dependence <- any(map_lgl(attr(terms[[i]], "aux.slots")[if(is_aux[i]) - 1L else TRUE], # Auxiliaries' first slot is their own.
                                           function(slot) terms[[which(aux_slots == slot)]]$dependence))
    }
  }

  terms
}
