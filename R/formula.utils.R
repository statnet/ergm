#  File R/formula.utils.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
###################################################################
## This file has utilities whose primary purpose is examining or ##
## manipulating ERGM formulas.                                   ##
###################################################################

model.transform.formula <- function(object, theta, recipes, ...){
  ## Recipe syntax:
  ##
  ## Recipes are a named list with the names representing a map from
  ## term name to recipe. Each recipe is a list with instructions
  ## about what should be done with a term with that name.
  ##
  #### General settings
  ##
  ## filter: A function that takes a list containing the term name and
  ## arguments to the term, and returns either TRUE or FALSE. If it
  ## returns FALSE, the term is not processed further and is treated
  ## like a term not in the recipes list. If the function is absent,
  ## TRUE is assumed (i.e., all terms with this name are processed).
  ##
  ## custom: A function that takes two arguments: a list comprising
  ## the term name and arguments to the term and a vector of model
  ## parameters for the term. It must return a list with named
  ## elements: `theta`, a vector of model parameters for the term to
  ## be added to the ouptut parameters and `term`, an unevaluated call
  ## (or a list with the term name being the first element and
  ## arguments being the remaining (named) elements) to be appended to
  ## the rest of the formula. *If present, overrides all settings
  ## listed below. Most of the time, you should use those.*
  ##
  #### Convenience settings
  ##
  ## name: A string or a function. If a string, the term will be
  ## renamed to this. If a function, it must take two arguments: a
  ## string with the current name of the term and a list with the term name
  ## and arguments to the term, and return a string with the new name.
  ##
  ## tocoef: Either a numeric vector of integers or a function. If
  ## numeric, gives the indices of the model parameters for the term
  ## that are to be copied directly the output parameter vector. If a
  ## function, the function must take a vector of model parameters for
  ## the term and a list with the term name and arguments to the term
  ## and return a numeric vector to be appended to the output
  ## parameters.
  ##
  ## toarg: A named list of either numeric vectors of integers or
  ## functions (may be heterogeneous). For each element in the list,
  ## the term's corresponding argument will be set to either the model
  ## parameters indicated by the indices (if numeric) or the return
  ## value (if function). As with the others, the function must take
  ## two arguments: a vector of model parameters for the term and a
  ## list with the term name and arguments to the term.
  ##
  ## constant: A named list of elements of any type. The corresponding
  ## arguments of the term are set to the values in the list. This is
  ## a simple special case of toarg, if it were given a function that
  ## returned a constant value.

  m <- ergm_model(object, ergm.getnetwork(object), ...)
  theta.inds<-cumsum(c(1,nparam(m, byterm=TRUE)))
  terms<-list_rhs.formula(object)
  form<-object
  ## This deletes the formula's RHS, and LHS becomes RHS (for the moment).
  form[[3]]<-NULL
  newtheta<-c()
  for(i in seq_along(terms)){
    if(!is.call(terms[[i]]) ||
       !(as.character(terms[[i]][[1]]) %in% names(recipes)) ||
       (!is.null(recipes[[as.character(terms[[i]][[1]])]]$filter) &&
        !recipes[[as.character(terms[[i]][[1]])]]$filter(as.list(terms[[i]])))){
      ## If it's not a call OR is a call but does not have a recipe OR
      ## does have a recipe, but the filter function says it should be
      ## skipped (e.g. it's already fixed), then just append it.
      form<-append_rhs.formula(form,list(terms[[i]]))
      newtheta<-c(newtheta,theta[theta.inds[i]:(theta.inds[i+1]-1)])
    }else{
      ## Otherwise, it gets complicated...
      recipe<-recipes[[as.character(terms[[i]][[1]])]]
      orig.list<-call.list<-as.list(terms[[i]])

      if(!is.null(recipe$custom)){
        ## Custom recipe
        out<-recipe$custom(orig.list,theta[theta.inds[i]:(theta.inds[i+1]-1)])
        newtheta<-c(newtheta,out$theta)
        form<-append_rhs.formula(form,
                                 if(is.call(out$term)) list(out$term)
                                 else as.call(c(as.name(out$term[[1]]),out$term[-1])))
      }else{

        if("" %in% names(call.list)[-1]) stop("Curved terms must have all their arguments passed by name.")
        
        ## Rename the term.
        if(!is.null(recipe$name)){
          call.list[[1]]<-
            if(is.function(recipe$name)) as.name(recipe$name(as.character(orig.list[[1]]),orig.list))
            else as.name(recipe$name)
        }
        ## Now, go through the arguments to be replaced:
        ## The constants:
        for(name in names(recipe$constant))
          call.list[[name]]<-recipe$constant[[name]]
        
        ## The elements of theta:
        for(name in names(recipe$toarg))
          call.list[[name]]<-
            if(is.function(recipe$toarg[[name]])) recipe$toarg[[name]](theta[theta.inds[i]:(theta.inds[i+1]-1)],orig.list)
            else theta[theta.inds[i]+recipe$toarg[[name]]-1]
        
        ## Now, add the newly rewritten call to the formula.
        form<-append_rhs.formula(form,list(as.call(call.list)))
        
        ## The parts that remain in theta:
        newtheta<-c(newtheta,
                    if(is.function(recipe$tocoef)) recipe$tocoef[[name]](theta[theta.inds[i]:(theta.inds[i+1]-1)],orig.list)
                    else theta[theta.inds[i]+recipe$tocoef-1])
      }
    }
  }
  list(formula=form,theta=newtheta)
}


## Convert a fitted curved ERGM (or its formula + theta) into a linear
## (fixed=TRUE) model with the curved parameters in theta substituted
## into the formula according to a set of recipes. Returns the new
## formula and the appropriate parameter vector.



#' Convert a curved ERGM into a corresponding "fixed" ERGM.
#' 
#' The generic \code{fix.curved} converts an [`ergm`] object or
#' formula of a model with curved terms to the variant in which the curved
#' parameters are fixed. Note that each term has to be treated as a special
#' case.
#' 
#' Some ERGM terms such as [`gwesp`][gwesp-ergmTerm] and [`gwdegree`][gwdegree-ergmTerm] have
#' two forms: a curved form, for which their decay or similar parameters are to
#' be estimated, and whose canonical statistics is a vector of the term's
#' components ([`esp(1)`][esp-ergmTerm], [`esp(2)`][esp-ergmTerm], \dots{} and
#' [`degree(1)`][degree-ergmTerm], [`degree(2)`][degree-ergmTerm], \dots{}, respectively) and
#' a "fixed" form where the decay or similar parameters are fixed, and whose
#' canonical statistic is just the term itself. It is often desirable to fit a
#' model estimating the curved parameters but simulate the "fixed" statistic.
#' 
#' This function thus takes in a fit or a formula and performs this mapping,
#' returning a "fixed" model and parameter specification.  It only works for
#' curved ERGM terms included with the \CRANpkg{ergm}
#' package. It does not work with curved terms not included in ergm.
#' 
#' @param object An [`ergm`] object or an ERGM formula. The curved
#' terms of the given formula (or the formula used in the fit) must have all of
#' their arguments passed by name.
#' @param \dots Unused at this time.
#' @return A list with the following components: \item{formula}{The "fixed"
#' formula.} \item{theta}{The "fixed" parameter vector.}
#' @seealso [ergm()], [simulate.ergm()]
#' @keywords model
#' @examples
#' 
#' \donttest{
#' \dontshow{
#' options(ergm.eval.loglik=FALSE)
#' }
#' data(sampson)
#' gest<-ergm(samplike~edges+gwesp(),
#'            control=control.ergm(MCMLE.maxit=2))
#' summary(gest)
#' # A statistic for esp(1),...,esp(16)
#' simulate(gest,output="stats")
#' 
#' tmp<-fix.curved(gest)
#' tmp
#' # A gwesp() statistic only
#' simulate(tmp$formula, coef=tmp$theta, output="stats") 
#' }
#' 
#' @export fix.curved
fix.curved <- function(object, ...) UseMethod("fix.curved")

#' @rdname fix.curved
#' @export
fix.curved.ergm <- function(object,...){
  fix.curved.formula(object$formula, coef(object), ...)
}

#' @rdname fix.curved
#' @param theta Curved model parameter configuration.
#' @export
fix.curved.formula <- function(object, theta, ...){
  recipes<-list()
  is.fixed.1<-function(a) is.null(a$fixed) || a$fixed==FALSE
  recipes$dgwdsp<-recipes$dgwesp<-recipes$dgwnsp<-recipes$gwdsp<-recipes$gwesp<-recipes$gwnsp<-recipes$gwb1dsp<-recipes$gwb2dsp<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2), constant=list(fixed=TRUE))
  recipes$altkstar<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(lambda=2), constant=list(fixed=TRUE))
  recipes$gwb1degree<-recipes$gwb2degree<-recipes$gwdegree<-recipes$gwidegree<-recipes$gwodegree<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2), constant=list(fixed=TRUE))

  model.transform.formula(object, theta, recipes, ...)
}


#' Convert a curved ERGM into a form suitable as initial values for the same
#' ergm. Deprecated in 4.0.0.
#' 
#' The generic \code{enformulate.curved} converts an [`ergm`] object
#' or formula of a model with curved terms to the variant in which the curved
#' parameters embedded into the formula and are removed from the parameter
#' vector. This is the form that used to be required by [ergm()] calls.
#' 
#' Because of a current kludge in [ergm()], output from one run
#' cannot be directly passed as initial values (\code{control.ergm(init=)}) for
#' the next run if any of the terms are curved. One workaround is to embed the
#' curved parameters into the formula (while keeping \code{fixed=FALSE}) and
#' remove them from \code{control.ergm(init=)}.
#' 
#' This function automates this process for curved ERGM terms included with the
#' \CRANpkg{ergm} package. It does not work with curved
#' terms not included in ergm.
#' 
#' @param object An [`ergm`] object or an ERGM formula. The curved
#' terms of the given formula (or the formula used in the fit) must have all of
#' their arguments passed by name.
#' @param \dots Unused at this time.
#' @return A list with the following components: \item{formula}{The formula
#' with curved parameter estimates incorporated.} \item{theta}{The coefficient
#' vector with curved parameter estimates removed.}
#' @seealso [ergm()], [simulate.ergm()]
#' @keywords model
#' @name enformulate.curved-deprecated
## #' @examples
## #' 
## #' \donttest{
## #' \dontshow{
## #' options(ergm.eval.loglik=FALSE)
## #' }
## #' data(sampson)
## #' gest<-ergm(samplike~edges+gwesp(decay=.5, fixed=FALSE), 
## #'     control=control.ergm(MCMLE.maxit=1))
## #' # Error:
## #' gest2<-try(ergm(gest$formula, control=control.ergm(init=coef(gest), MCMLE.maxit=1)))
## #' print(gest2)
## #' 
## #' # Works:
## #' tmp<-enformulate.curved(gest)
## #' tmp
## #' gest2<-try(ergm(tmp$formula, control=control.ergm(init=tmp$theta, MCMLE.maxit=1)))
## #' summary(gest2)
## #' }
#' 
#' @export enformulate.curved
enformulate.curved <- function(object, ...) UseMethod("enformulate.curved")

#' @rdname enformulate.curved-deprecated
#' @export
enformulate.curved.ergm <- function(object,...){
  .Deprecated(msg="enformulate.curved() family of functions has been obviated by native handling of curved ERGM terms.")
  fix.curved.formula(object$formula, coef(object), ...)
}

#' @rdname enformulate.curved-deprecated
#' @param theta Curved model parameter configuration.
#' @export
enformulate.curved.formula <- function(object, theta, ...){
  recipes<-list()
  is.fixed.1<-function(a) is.null(a$fixed) || a$fixed==FALSE
  recipes$dgwdsp<-recipes$dgwesp<-recipes$dgwnsp<-recipes$gwdsp<-recipes$gwesp<-recipes$gwnsp<-recipes$gwb1dsp<-recipes$gwb2dsp<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2))
  recipes$altkstar<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(lambda=2))
  recipes$gwb1degree<-recipes$gwb2degree<-recipes$gwdegree<-recipes$gwidegree<-recipes$gwodegree<-
    list(filter=is.fixed.1, tocoef=1, toarg=list(decay=2))

  model.transform.formula(object, theta, recipes, ...)
}

#' @describeIn ergm-deprecated \code{offset.info.formula} returns the offset
#'   vectors associated with a formula.
#' @export offset.info.formula
offset.info.formula <- function(object, ...){
  .Deprecated()
  nw <- ergm.getnetwork(object)
  m<-ergm_model(object, nw, ...)
  with(m$etamap, list(term=offset, theta=offsettheta,eta=offsetmap))
}
