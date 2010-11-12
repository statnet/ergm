InitWtErgmTerm.sum<-function(nw, arglist, drop=TRUE, response=NULL, ...) {
  if(is.null(response)) stop('Statistic "sum" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
##  if(drop) { # Check for zero statistics, print -Inf messages if applicable
##    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
##    if (extremewarnings(obsstats)) {
##      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
##    } 
##  }
  list(name="sum",
       coef.names="sum",
       inputs=NULL,
       dependence=FALSE)
}

InitWtErgmTerm.nonzero<-function(nw, arglist, drop=TRUE, response=NULL, ...) {
  if(is.null(response)) stop('Statistic "nonzero" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
##  if(drop) { # Check for zero statistics, print -Inf messages if applicable
##    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
##    if (extremewarnings(obsstats)) {
##      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
##    } 
##  }

  list(name="nonzero",
       coef.names="nonzero",
       inputs=NULL,
       dependence=FALSE)
}

InitWtErgmTerm.mutualweights<-function (nw, arglist, drop=TRUE, response=NULL, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("min"),
                      required = c(FALSE))
  
##  if(drop) { # Check for zero statistics, print -Inf messages if applicable
##    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
##    if (extremewarnings(obsstats)) {
##      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
##    } 
##  }
  form<-match.arg(a$form,c("min","nabsdiff"))
  list(name=switch(form,min="mutualweights_min",nabsdiff="mutualweights_nabsdiff"),
       coef.names=switch(form,min="mutualweights.min",nabsdiff="mutualweights.nabsdiff"),
       inputs=NULL,
       dependence=FALSE)
  
}

InitWtErgmTerm.transitiveweights<-function (nw, arglist, drop=TRUE, response=NULL, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("max"),
                      required = c(FALSE))
  
##  if(drop) { # Check for zero statistics, print -Inf messages if applicable
##    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
##    if (extremewarnings(obsstats)) {
##      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
##    } 
##  }
  form<-match.arg(a$form,c("max","sum"))
  list(name=switch(form,max="transitiveweights_max",sum="transitiveweights_sum"),
       coef.names=switch(form,max="transitiveweights.max",sum="transitiveweights.sum"),
       inputs=NULL,
       dependence=FALSE)
  
}

InitWtErgmTerm.cyclicweights<-function (nw, arglist, drop=TRUE, response=NULL, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,  bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("max"),
                      required = c(FALSE))
  
##  if(drop) { # Check for zero statistics, print -Inf messages if applicable
##    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
##    if (extremewarnings(obsstats)) {
##      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
##    } 
##  }
  form<-match.arg(a$form,c("max","sum"))
  list(name=switch(form,max="cyclicweights_max",sum="cyclicweights_sum"),
       coef.names=switch(form,max="cyclicweights.max",sum="cyclicweights.sum"),
       inputs=NULL,
       dependence=FALSE)
  
}
