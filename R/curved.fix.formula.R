## Converts a fitted curved ERGM (or its formula + theta) into a
## linear (fixed=TRUE) model with the curved parameters in theta
## substituted into the formula according to a set of recipes. Returns
## the new formula and the appropriate parameter vector.

curved.fix <- function(object, ...) UseMethod("curved.fix")

curved.fix.ergm <- function(object,...){
  curved.fix.formula(object$formula, coef(object), response=object$response, ...)
}

curved.fix.formula <- function(object, theta, response=NULL, ...){
  ## Construct the recipes.
  recipes<-list()
  is.fixed.1<-function(a) is.null(a$fixed) || a$fixed==FALSE
  recipes$gwdsp<-recipes$gwesp<-recipes$gwnsp<-
    list(filter=is.fixed.1, coef=1, theta=list(alpha=2), constant=list(fixed=TRUE))
  recipes$altkstar<-
    list(filter=is.fixed.1, coef=1, theta=list(lambda=2), constant=list(fixed=TRUE))
  recipes$gwb1degree<-recipes$gwb2degree<-recipes$gwdegree<-recipes$gwidegree<-recipes$gwodegree<-
    list(filter=is.fixed.1, coef=1, theta=list(decay=2), constant=list(fixed=TRUE))

  m <- ergm.getmodel(object, ergm.getnetwork(object), drop=FALSE, response=response)
  theta.inds<-cumsum(c(1,theta.sublength.model(m)))
  terms<-term.list.formula(object[[3]])
  form<-object
  ## This deletes the formula's RHS, and LHS becomes RHS (for the moment).
  form[[3]]<-NULL
  newtheta<-c()
  for(i in seq_along(terms)){
    if(!is.call(terms[[i]]) ||
       !(as.character(terms[[i]][[1]]) %in% names(recipes)) ||
       (!is.null(recipes[[as.character(terms[[i]][[1]])]]$filter) &&
        !recipes[[as.character(terms[[i]][[1]])]]$filter(as.list(terms[[i]])[-1]))){
      ## If it's not a call OR is a call but does not have a recipe OR
      ## does have a recipe, but the filter function says it should be
      ## skipped (e.g. it's already fixed), then just append it.
      form<-append.rhs.formula(form,list(terms[[i]]))
      newtheta<-c(newtheta,theta[theta.inds[i]:(theta.inds[i+1]-1)])
    }else{
      ## Otherwise, it gets complicated...
      recipe<-recipes[[as.character(terms[[i]][[1]])]]
      call.list<-as.list(terms[[i]])

      if("" %in% names(call.list)[-1]) stop("Curved terms must have all their arguments passed by name.")
      
      ## Rename the term.
      if(!is.null(recipe$name)) call.list[[1]]<-as.name(recipe$name)
      ## Now, go through the arguments to be replaced:
      ## The constants:
      for(name in names(recipe$constant))
        call.list[[name]]<-recipe$constant[[name]]
      ## The elements of theta
      for(name in names(recipe$theta)){
        if(is.numeric(recipe$theta[[name]]))
          ## Position in theta
          call.list[[name]]<-theta[theta.inds[i]+recipe$theta[[name]]-1]
        else
          ## Function of theta
          call.list[[name]]<-recipe$theta[[name]](theta[theta.inds[i]:(theta.inds[i+1]-1)])
      }

      ## Now, add the newly rewritten call to the formula.
      form<-append.rhs.formula(form,list(as.call(call.list)))
      
      ## The parts that remain in theta.
      newtheta<-c(newtheta,theta[theta.inds[i]+recipe$coef-1])
    }
  }
  list(formula=form,theta=newtheta)
}
