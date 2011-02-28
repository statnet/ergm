InitWtErgmTerm.absdiff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow","form"),
                      vartypes = c("character","numeric","character"),
                      defaultvalues = list(NULL,1,"sum"),
                      required = c(TRUE,FALSE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  ### Construct the list to return
  form<-match.arg(a$form,c("sum","nonzero"))
  list(name=paste("absdiff",form,sep="_"),                                     #name: required
       coef.names = paste(paste("absdiff",if(a$pow!=1) a$pow else "",sep=""), form,  a$attrname, sep="."), #coef.names: required
       inputs = c(a$pow,nodecov),  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

InitWtErgmTerm.absdiffcat <- function(nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","base","form"),
                      vartypes = c("character","numeric","character"),
                      defaultvalues = list(NULL,NULL,"sum"),
                      required = c(TRUE,FALSE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=TRUE)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(a$base)) u <- u[-(a$base)]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  u2 <- u[!is.na(u)]
  ### Construct the list to return
  inputs <- c(u2, NAsubstitute, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(u2)+1 # See comment at top of file
  form<-match.arg(a$form,c("sum","nonzero"))
  list(name=paste("absdiffcat",form,sep="_"),                                  #name: required
       coef.names = paste("absdiff",form, a$attrname, u, sep="."), #coef.names: required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
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
                      varnames = c("pow"),
                      vartypes = c("numeric"),
                      defaultvalues = list(1),
                      required = c(FALSE))
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, response=response, ...)
    if (extremewarnings(obsstats)) {
      return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
    } 
  }

  if(a$pow==1){
    list(name="sum",
         coef.names="sum",
         inputs=NULL,
         dependence=FALSE)
  }else{
    list(name="sum_pow",
         coef.names=paste("sum",a$pow,sep=""),
         inputs=a$pow,
         dependence=FALSE)
  }
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

InitWtErgmTerm.nodeofactor<-function (nw, arglist, response, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
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
      print("Warning:  nodeofactor term deleted because it contributes no statistics")
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
  list(name=paste("nodeofactor",form,sep="_"),                                        #required
       coef.names = paste("nodeofactor",form, paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

InitWtErgmTerm.nodeifactor<-function (nw, arglist, response, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
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
      print("Warning:  nodeifactor term deleted because it contributes no statistics")
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
  list(name=paste("nodeifactor",form,sep="_"),                                        #required
       coef.names = paste("nodeifactor",form, paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

InitWtErgmTerm.nodecov<-InitWtErgmTerm.nodemain<-function (nw, arglist, response, drop=TRUE, ...) {
  a <- ergm.ErgmTerm(nw, arglist,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  inputs <- f(get.node.attr(nw, attrname, "nodecov", numeric=TRUE))
  attr(inputs, "ParamsBeforeCov")<-length(inputs)
  form<-match.arg(a$form,c("sum","nonzero"))
  list(name=paste("nodecov",form,sep="_"), soname="ergm",
       coef.names = paste("nodecov",form,f.name,attrname,sep="."),
       inputs=inputs,
       dependence=FALSE)
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
       dependence=TRUE)
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
       inputs=c(a$threshold),
       dependence=TRUE)
  
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
       dependence=TRUE)
  
}

InitWtErgmTerm.cyclicalties<-function (nw, arglist, response, drop=TRUE, ...) {
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
       coef.names="cyclicalties",
       inputs=c(a$threshold),
       dependence=TRUE)
  
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
       dependence=TRUE)
  
}
