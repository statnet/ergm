ergm.auxstorage <- function(model, nw, response=NULL,...){
  # As formulas
  aux.forms <- lapply(model$terms, "[[", "auxiliaries")

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
  uniq.aux.outlists <- unique(unlist(term.outlists))

  # Initialize the auxiliary model.
  aux.model <- structure(list(formula=formula, coef.names = NULL,
                          offset = NULL,
                          terms = NULL, networkstats.0 = NULL, etamap = NULL),
                         class = "model.ergm")
  for(i in seq_along(uniq.aux.outlists)){
    aux.outlist <- uniq.aux.outlists[[i]]
    aux.model <- updatemodel.ErgmTerm(aux.model, aux.outlist)
    aux.model[[i]]$inputs[2] <- -(i-1) # The storage slot belonging to this auxiliary (replacing the number of stats.
  }

  # Which term is requiring which auxiliary slot? (+1)
  aux.slots <- lapply(term.outlists, match, uniq.aux.outlists)

  for(i in seq_along(model$terms))
    if(length(aux.outlists[[i]]))
      model$terms[[i]]$inputs[3+seq_len(length(aux.outlists[[i]]))] <- aux.slots[[i]]-1

  model$model.aux <- aux.model
  model
}
