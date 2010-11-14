InitWtErgmTerm.CMP<-function(nw, arglist, response, drop=TRUE, ...) {
  if(is.null(response)) stop('Conway-Maxwell-Poisson form requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=NULL,maxval=0)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="nsumlogfactorial",
       coef.names="CMP",
       inputs=NULL,
       dependence=FALSE)
}

InitWtErgmTerm.sum<-function(nw, arglist, response, drop=TRUE, ...) {
  if(is.null(response)) stop('Statistic "sum" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="sum",
       coef.names="sum",
       inputs=NULL,
       dependence=FALSE)
}

InitWtErgmTerm.nodefactor<-function (nw, arglist, response, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments

  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      print("Warning:  nodefactor term deleted because it contributes no statistics")
      return()
    }
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    ew <- extremewarnings(obsstats)
    u <- u[!ew]
    ui <- ui[!ew]
  }
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file R/InitErgmTerm.R
  list(name="nodefactor_wt",                                        #required
       coef.names = paste("nodefactor", paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  

InitWtErgmTerm.nonzero<-function(nw, arglist, response, drop=TRUE, ...) {
  if(is.null(response)) stop('Statistic "nonzero" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }

  list(name="nonzero",
       coef.names="nonzero",
       inputs=NULL,
       dependence=FALSE)
}

InitWtErgmTerm.mutual<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("min"),
                      required = c(FALSE))

  form<-match.arg(a$form,c("min","nabsdiff"))
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,
                        minval=switch(form,min=0,nabsdiff=NULL),
                        maxval=switch(form,min=NULL,nabsdiff=0))
        ) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name=switch(form,min="mutual_wt_min",nabsdiff="mutual_wt_nabsdiff"),
       coef.names=switch(form,min="mutual.min",nabsdiff="mutual.nabsdiff"),
       inputs=NULL,
       dependence=FALSE)
  
}

InitWtErgmTerm.transitiveweights<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("max"),
                      required = c(FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  form<-match.arg(a$form,c("max","sum"))
  list(name=switch(form,max="transitiveweights_max",sum="transitiveweights_sum"),
       coef.names=switch(form,max="transitiveweights.max",sum="transitiveweights.sum"),
       inputs=NULL,
       dependence=FALSE)
  
}

InitWtErgmTerm.cyclicweights<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,  bipartite=NULL,
                      varnames = c("form"),
                      vartypes = c("character"),
                      defaultvalues = list("max"),
                      required = c(FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  form<-match.arg(a$form,c("max","sum"))
  list(name=switch(form,max="cyclicweights_max",sum="cyclicweights_sum"),
       coef.names=switch(form,max="cyclicweights.max",sum="cyclicweights.sum"),
       inputs=NULL,
       dependence=FALSE)
  
}
