#  File R/ergm.auxstorage.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
ergm.auxstorage <- function(model, nw, response=NULL,..., extra.aux=list()){
  # As formulas
  aux.forms <- c(lapply(model$terms, "[[", "auxiliaries"), extra.aux)

  # A nested list: outer list are the model terms requesting
  # auxiliaries and the inner list is the outputs from InitErgmTerm
  # calls of the auxiliaries.
  aux.outlists <- lapply(aux.forms, function(aux.form){
    if(is.null(aux.form)) list()
    else{
      formula.env <- environment(aux.form)
      lapply(term.list.formula(aux.form[[length(aux.form)]]), function(aux.term){
        call.ErgmTerm(aux.term, formula.env, nw, response=response,...)
      })
    }
  })

  # Remove duplicated auxiliaries.
  uniq.aux.outlists <- unique(unlist(aux.outlists, recursive=FALSE))

  # Initialize the auxiliary model.
  aux.model <- structure(list(formula=NULL, coef.names = NULL,
                          offset = NULL,
                          terms = NULL, networkstats.0 = NULL, etamap = NULL),
                         class = "ergm.model")
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
 
  model$model.aux <- aux.model
  model$slots.extra.aux <- slots.extra.aux
  model
}
