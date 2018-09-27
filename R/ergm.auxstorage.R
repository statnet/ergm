#  File R/ergm.auxstorage.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
ergm.auxstorage <- function(model, nw, response=NULL,..., extra.aux=list(), term.options=list()){

  aux_list_list <- function(terms, extra=NULL) {
    # As formulas
    aux.forms <- c(lapply(terms, "[[", "auxiliaries"), extra)
    
    # A nested list: outer list are the model terms requesting
    # auxiliaries and the inner list is the outputs from InitErgmTerm
    # calls of the auxiliaries.
    lapply(aux.forms, function(aux.form){
      if(is.null(aux.form)) list()
      else{
        formula.env <- environment(aux.form)
        lapply(list_rhs.formula(aux.form), function(aux.term){
          call.ErgmTerm(aux.term, formula.env, nw, response=response,term.options=term.options,...)
        })
      }
    })
  }

  aux.outlists <- aux_list_list(model$terms, extra.aux)

  # Remove duplicated auxiliaries.
  uniq.aux.outlists <- unique(unlist(aux.outlists, recursive=FALSE), fromLast=TRUE)
  prev <- 0
  aux.aux.outlists <- list()

  # Until we reach a fixed point (which we should, unless there is a circular dependency.
  #
  # TODO: Check for circular dependencies.
  while(length(uniq.aux.outlists)!=prev){
    prev <- length(uniq.aux.outlists)
    aux.aux.outlists <- aux_list_list(uniq.aux.outlists)
    uniq.aux.outlists <- unique(c(uniq.aux.outlists, unlist(aux.aux.outlists, recursive=FALSE)), fromLast=TRUE)
  }

  # uniq.aux.outlists is now a list of unique initialized auxiliaries
  # and auxiliaries' auxiliaries, such that depended-on auxiliaries
  # are always after the dependent auxiliaries.
  #
  # aux.aux.outlists is now a nested list in the same form as aux.outlists, but for auxiliaries.

  # Initialize the auxiliary model.
  aux.model <- structure(list(formula=NULL, coef.names = NULL,
                          offset = NULL,
                          terms = NULL, networkstats.0 = NULL, etamap = NULL),
                         class = "ergm_model")
  for(i in seq_along(uniq.aux.outlists)){
    aux.outlist <- uniq.aux.outlists[[i]]
    aux.model <- updatemodel.ErgmTerm(aux.model, aux.outlist)
    aux.model$terms[[i]]$inputs[4] <- i-1 # The storage slot belonging to this auxiliary.
  }

  # Which term is requiring which auxiliary slot? (+1)
  aux.slots <- lapply(aux.outlists, match, uniq.aux.outlists)
  slots.extra.aux <- list()
  
  for(i in seq_along(aux.outlists)){
    if(length(aux.outlists[[i]])){
      if(i<=length(model$terms)) # If it's a model term.
        model$terms[[i]]$inputs[3+seq_len(length(aux.outlists[[i]]))] <- aux.slots[[i]]-1
      else # If it's some other entity requesting auxiliaries.
        slots.extra.aux[[i-length(model$terms)]] <- aux.slots[[i]]-1
    }
  }

  # Which auxiliary is requiring which auxiliary slot? (+1)
  aux.aux.slots <- lapply(aux.aux.outlists, match, uniq.aux.outlists)
  
  for(i in seq_along(aux.aux.outlists)){
    if(length(aux.aux.outlists[[i]])){
      # 4th input is auxiliary's own slot, so its auxiliaries get put into subsequent slots.
      aux.model$terms[[i]]$inputs[4+seq_len(length(aux.aux.outlists[[i]]))] <- aux.aux.slots[[i]]-1
    }
  }
  
  model$model.aux <- aux.model
  model$slots.extra.aux <- slots.extra.aux
  model
}
