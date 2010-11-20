InitWtErgmTerm.CMP<-function(nw, arglist, response, drop=TRUE, ...) {
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

InitWtErgmTerm.ininterval<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("lower","upper","open"),
                      vartypes = c("numeric","numeric","logical"),
                      defaultvalues = list(-Inf,+Inf,c(TRUE,TRUE)),
                      required = c(FALSE,FALSE,FALSE))

  a$open<-rep(a$open,length.out=2)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="ininterval",
       coef.names=paste("ininterval",if(a$open[0]) "(" else "[", a$lower,",",a$upper, if(a$open[1]) ")" else "]",sep=""),
       inputs=c(a$lower,a$upper,a$open),
       dependence=FALSE)
}


InitWtErgmTerm.atleast<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="atleast",
       coef.names=paste("atleast",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE)
}


InitWtErgmTerm.greaterthan<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="atleast",
       coef.names=paste("greaterthan",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE)
}



InitWtErgmTerm.sum<-function(nw, arglist, response, drop=TRUE, ...) {
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
                      varnames = c("attrname", "base","form"),
                      vartypes = c("character", "numeric","character"),
                      defaultvalues = list(NULL, 1, "sum"),
                      required = c(TRUE, FALSE, FALSE))
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
  form<-match.arg(a$form,c("sum","nonzero"))
  list(name=paste("nodefactor",form,sep="_"),                                        #required
       coef.names = paste("nodefactor",form, paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  

InitWtErgmTerm.nonzero<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
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
                      varnames = c("form","threshold"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("min",0),
                      required = c(FALSE,FALSE))

  form <- match.arg(a$form,c("min","nabsdiff","threshold"))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,
                        minval=switch(form,min=0,nabsdiff=NULL,threshold=0),
                        maxval=switch(form,min=NULL,nabsdiff=0,threshold=NULL))
        ) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    }
  }
  list(name=switch(form,min="mutual_wt_min",nabsdiff="mutual_wt_nabsdiff",threshold="mutual_wt_threshold"),
       coef.names=switch(form,min="mutual.min",nabsdiff="mutual.nabsdiff",threshold=paste("mutual",a$threshold,sep=".")),
       inputs=if(form=="threshold") threshold,
       dependence=FALSE)
}

InitWtErgmTerm.transitiveties<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="transitiveweights_threshold",
       coef.names="transitiveties",
       inputs=c(0),
       dependence=FALSE)
  
}

InitWtErgmTerm.transitiveweights<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("form","threshold"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("max",0),
                      required = c(FALSE,FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  form<-match.arg(a$form,c("max","sum","threshold"))
  list(name=switch(form,max="transitiveweights_max",sum="transitiveweights_sum",threshold="transitiveweights_threshold"),
       coef.names=switch(form,max="transitiveweights.max",sum="transitiveweights.sum",threshold=paste("transitiveweights",a$threshold,sep=".")),
       inputs=if(form=="threshold") threshold,
       dependence=FALSE)
  
}

InitWtErgmTerm.cyclicties<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats,minval=0,maxval=network.dyadcount(nw,TRUE))) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  list(name="cyclicweights_threshold",
       coef.names="cyclicties",
       inputs=c(0),
       dependence=FALSE)
  
}

InitWtErgmTerm.cyclicweights<-function (nw, arglist, response, drop=TRUE, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("form","threshold"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("max",0),
                      required = c(FALSE,FALSE))
  
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }
  form<-match.arg(a$form,c("max","sum","threshold"))
  list(name=switch(form,max="cyclicweights_max",sum="cyclicweights_sum",threshold="cyclicweights_threshold"),
       coef.names=switch(form,max="cyclicweights.max",sum="cyclicweights.sum",threshold=paste("cyclicweights",a$threshold,sep=".")),
       inputs=if(form=="threshold") threshold,
       dependence=FALSE)
  
}
