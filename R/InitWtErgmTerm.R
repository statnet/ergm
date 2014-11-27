#  File R/InitWtErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
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


InitWtErgmTerm.atleast<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="atleast",
       coef.names=paste("atleast",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
}

InitWtErgmTerm.edgecov <- function(nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname", "form"),
                      vartypes = c("matrix,network,character", "character", "character"),
                      defaultvalues = list(NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix.network(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), 
              as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), sep = ".")
  }

  form<-match.arg(a$form,c("sum","nonzero"))
  
  inputs <- c(as.double(xm))
  list(name=paste("edgecov",form,sep="_"), coef.names = paste(cn,form,sep="."), inputs = inputs, dependence=FALSE)
}


InitWtErgmTerm.ininterval<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("lower","upper","open"),
                      vartypes = c("numeric","numeric","logical"),
                      defaultvalues = list(-Inf,+Inf,c(TRUE,TRUE)),
                      required = c(FALSE,FALSE,FALSE))

  a$open<-rep(a$open,length.out=2)
  list(name="ininterval",
       coef.names=paste("ininterval",if(a$open[0]) "(" else "[", a$lower,",",a$upper, if(a$open[1]) ")" else "]",sep=""),
       inputs=c(a$lower,a$upper,a$open),
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
}

InitWtErgmTerm.greaterthan<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="atleast",
       coef.names=paste("greaterthan",a$threshold,sep="."),
       inputs=a$threshold,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
}



InitWtErgmTerm.sum<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("pow"),
                      vartypes = c("numeric"),
                      defaultvalues = list(1),
                      required = c(FALSE))
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

InitWtErgmTerm.nodecovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite = FALSE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodecovar",
       coef.names = "nodecovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodesqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite = FALSE, nonnegative=TRUE, response=response,
                      varnames = c("center"),
                      vartypes = c("logical"),
                      defaultvalues = list(TRUE),
                      required = c(TRUE))
  ### Process the arguments

  name <- "nodesqrtcovar"

  if(a$center) name <- paste(name,"centered",sep="_")

  coef.name <- gsub("_",".",name)
  
  list(name=name,
       coef.names = coef.name,
       dependence = TRUE,
       minval = if(a$center) NULL else 0
       )
}

InitWtErgmTerm.nodeosqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE, nonnegative=TRUE, response=response,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  ### Process the arguments

  list(name="nodeosqrtcovar",
       coef.names = "nodeosqrtcovar",
       dependence = TRUE,
       # arithmetic mean >= geometric mean
       minval = 0
       )
}

InitWtErgmTerm.nodeisqrtcovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE, nonnegative=TRUE, response=response,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(FALSE),
                      required = c(FALSE))
  ### Process the arguments

  list(name="nodeisqrtcovar",
       coef.names = "nodeisqrtcovar",
       dependence = TRUE,
       # arithmetic mean >= geometric mean
       minval = 0
       )
}

InitWtErgmTerm.nodefactor<-function (nw, arglist, response, ...) {
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

InitWtErgmTerm.nodeocovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodeocovar",
       coef.names = "nodeocovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeofactor<-function (nw, arglist, response, ...) {
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

InitWtErgmTerm.nodeicovar<-function (nw, arglist, response, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed = TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = NULL,
                      required = NULL)
  ### Process the arguments

  list(name="nodeicovar",
       coef.names = "nodeicovar",
       dependence = TRUE
       )
}

InitWtErgmTerm.nodeifactor<-function (nw, arglist, response, ...) {
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

InitWtErgmTerm.nodematch<-InitWtErgmTerm.match<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "diff", "keep", "form"),
                      vartypes = c("character", "logical", "numeric", "character"),
                      defaultvalues = list(NULL, FALSE, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE, FALSE))
  ### Process the arguments
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
  u <- sort(unique(nodecov))
  if (!is.null(a$keep)) {
    u <- u[a$keep]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov==(length(u)+1)
  nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
  ui <- seq(along=u)

  form<-match.arg(a$form,c("sum","nonzero"))

  ### Construct the list to return
  if (a$diff) {
    coef.names <- paste("nodematch", form, paste(a$attrname,collapse="."), u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", form, paste(a$attrname,collapse="."), sep=".")
    inputs <- nodecov
  }
  list(name=paste("nodematch",form,sep="_"),                                 #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

InitWtErgmTerm.nodemix<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "base", "form"),
                      vartypes = c("character", "numeric", "character"),
                      defaultvalues = list(NULL, NULL, "sum"),
                      required = c(TRUE, FALSE, FALSE))

  form<-match.arg(a$form,c("sum","nonzero"))

  ### Process the arguments
  if (is.bipartite(nw) && is.directed(nw)) {
    stop("Directed bipartite networks are not currently possible")
  }
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
    
  if (is.bipartite(nw)) {
    #  So undirected network storage but directed mixing
    nb1 <- get.network.attribute(nw, "bipartite")       
    #  Recode nodecov to numeric (but retain original sorted names in "namescov")
    b1namescov <- sort(unique(nodecov[1:nb1]))
    b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
    namescov <- c(b1namescov, b2namescov)
    b1nodecov <- match(nodecov[1:nb1],b1namescov)
    b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    nodecov <- c(b1nodecov, b2nodecov + nr)
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}    
    if (!is.null(a$base) && !identical(a$base,0)) {
      u <- u[-a$base,]
    }
    name <- "mix"
    cn <- paste("mix", form, paste(a$attrname,collapse="."), apply(matrix(namescov[u],ncol=2),
                                       1,paste,collapse="."), sep=".")
    inputs <- c(u[,1], u[,2], nodecov)
    attr(inputs, "ParamsBeforeCov") <- NROW(u)
  } else {# So one mode, but could be directed or undirected
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    #   Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    ucount<-sapply(ui,function(x){sum(nodecov==x,na.rm=TRUE)}) #Count cases
    uui <- matrix(1:length(ui)^2,length(ui),length(ui))  #Create int tables
    urm <- t(sapply(ui,rep,length(ui)))   #This is the reverse of what you'd
    ucm <- sapply(ui,rep,length(ui))      #expect for r/c, but it's correct
    uun <- outer(u,u,paste,sep=".")
    if (!is.directed(nw)) {
      uui <- uui[upper.tri(uui,diag=TRUE)]
      urm <- urm[upper.tri(urm,diag=TRUE)]  
      ucm <- ucm[upper.tri(ucm,diag=TRUE)]
      uun <- uun[upper.tri(uun,diag=TRUE)]
    }
    if (!is.null(a$base) && !identical(a$base,0)) {
      urm <- as.vector(urm)[-a$base]
      ucm <- as.vector(ucm)[-a$base]
      uun <- as.vector(uun)[-a$base]
    }
    name <- "nodemix"
    cn <- paste("mix", form, paste(a$attrname,collapse="."), uun, sep=".")
    inputs <- c(urm, ucm, nodecov)
    #attr(inputs, "ParamsBeforeCov") <- 2*length(uui)
    attr(inputs, "ParamsBeforeCov") <- 2*length(uun)
  }
  ### Construct the list to return
  list(name = paste(name, form, sep="_"), coef.names = cn, # required
       inputs = inputs, 
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

InitWtErgmTerm.nodecov<-InitWtErgmTerm.nodemain<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  form <- match.arg(a$form,c("sum","nonzero"))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste("nodecov",form,f.name,attrname,sep=".")
  nodecov <- f(get.node.attr(nw, attrname, "nodecov", numeric=TRUE))
  list(name=paste("nodecov",form,sep="_"), soname="ergm",
       coef.names=coef.names,
       inputs=c(nodecov),
       dependence=FALSE)
}

InitWtErgmTerm.nodeicov<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  form<-match.arg(a$form,c("sum","nonzero"))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste("nodeicov",form,f.name,attrname,sep=".")
  nodecov <- f(get.node.attr(nw, attrname, "nodeicov", numeric=TRUE))
  list(name=paste("nodeicov",form,sep="_"), soname="ergm",
       coef.names=coef.names,
       inputs=c(nodecov),
       dependence=FALSE)
}


InitWtErgmTerm.nodeocov<-function (nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                     varnames = c("attrname","transform","transformname","form"),
                     vartypes = c("character","function","character","character"),
                     defaultvalues = list(NULL,identity,"","sum"),
                     required = c(TRUE,FALSE,FALSE,FALSE))
  form<-match.arg(a$form,c("sum","nonzero"))
  attrname<-a$attrname
  f<-a$transform
  f.name<-a$transformname
  coef.names <- paste("nodeocov",form,f.name,attrname,sep=".")
  nodecov <- f(get.node.attr(nw, attrname, "nodeicov", numeric=TRUE))
  list(name=paste("nodeocov",form,sep="_"), soname="ergm",
       coef.names=coef.names,
       inputs=c(nodecov),
       dependence=FALSE)
}

InitWtErgmTerm.edges<-InitWtErgmTerm.nonzero<-function(nw, arglist, response, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nonzero",
       coef.names="nonzero",
       inputs=NULL,
       dependence=FALSE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
}

InitWtErgmTerm.mutual<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("form","threshold"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list("min",0),
                      required = c(FALSE,FALSE))

  form <- match.arg(a$form,c("min","nabsdiff","threshold","product","geometric"))
  
  list(name=switch(form,min="mutual_wt_min",nabsdiff="mutual_wt_nabsdiff",threshold="mutual_wt_threshold",product="mutual_wt_product", geometric="mutual_wt_geom_mean"),
       coef.names=switch(form,min="mutual.min",nabsdiff="mutual.nabsdiff",threshold=paste("mutual",a$threshold,sep="."), product="mutual.product",geometric="mutual.geom.mean"),
       inputs=if(form=="threshold") a$threshold,
       dependence=TRUE,
       minval=switch(form,min=NULL,nabsdiff=NULL,threshold=0,product=NULL,geometric=0),
       maxval=switch(form,min=NULL,nabsdiff=0,threshold=NULL,product=NULL,geometric=NULL)
       )
}

InitWtErgmTerm.transitiveties<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  
  list(name="transitiveweights_threshold",
       coef.names="transitiveties",
       inputs=c(a$threshold),
       dependence=TRUE,
       minval=0, maxval=network.dyadcount(nw,FALSE))  
}

InitWtErgmTerm.transitiveweights<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE), response=response)
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max","sum")
  combine<-match.arg(a$combine,combines)
  affects<-c("min","geomean")
  affect<-match.arg(a$affect,affects)

  list(name="transitiveweights",
       coef.names=paste("transitiveweights",twopath,combine,affect,sep="."),
       inputs=c(
         which(twopaths==twopath),
         which(combines==combine),
         which(affects==affect)
         ),
       dependence=TRUE,
       minval = 0)
}

InitWtErgmTerm.cyclicalties<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL,
                      varnames = c("threshold"),
                      vartypes = c("numeric"),
                      defaultvalues = list(0),
                      required = c(FALSE))
  list(name="cyclicweights_threshold",
       coef.names="cyclicalties",
       inputs=c(a$threshold),
       dependence=TRUE,
       minval=0, maxval=network.dyadcount(nw,FALSE))
  
}

InitWtErgmTerm.cyclicalweights<-function (nw, arglist, response, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, bipartite=NULL, nonnegative=TRUE,
                      varnames = c("twopath","combine","affect"),
                      vartypes = c("character","character","character"),
                      defaultvalues = list("min","max","min"),
                      required = c(FALSE,FALSE), response=response)
  twopaths<-c("min","geomean")
  twopath<-match.arg(a$twopath,twopaths)
  combines<-c("max","sum")
  combine<-match.arg(a$combine,combines)
  affects<-c("min","geomean")
  affect<-match.arg(a$affect,affects)

  list(name="cyclicalweights",
       coef.names=paste("cyclicalweights",twopath,combine,affect,sep="."),
       inputs=c(
         which(twopaths==twopath),
         which(combines==combine),
         which(affects==affect)
         ),
       dependence=TRUE,
       minval = 0)
}

