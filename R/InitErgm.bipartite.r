#  See InitErgm.r for a general explanation 
#  of InitErgm functions

# DH:  These still need to be translated into the new
# InitErgm format!

ergm.checkbipartite <- function(fname, nw.bipartiteflag, requirement,
                               extramessage="") {
  if (!nw.bipartiteflag && requirement)
    stop(paste(fname, "model term may not be used with an non-bipartite network.",
               extramessage), call.=FALSE)
  if (nw.bipartiteflag && !requirement)
    stop(paste(fname, "model term may not be used with a bipartite network.",
               extramessage), call.=FALSE)
}

##
## Fixed effects
##
#InitErgm.event<-function(g, model, drop=TRUE, ...)
#{
#    nevents <- is.bipartite(g)
#    if (!nevents)
#      stop("The event term is for bipartite graphs.",
#           call.=FALSE)
#    if (nargs()!=3)
#        stop(paste("event model term expected zero argument, got ", 
#            nargs()-2, sep=""), call.=FALSE)
##
##   Check for degeneracy
## 
#    nactors <- get.network.attribute(g,"bipartite")
#    nevents <- network.size(g)-nactors
#    xnames <- network.vertex.names(g)
#    d <- nactors + (2:nevents)
#    if(is.null(xnames)){
#     dnames <- paste("event",d-nactors,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    if(drop){
#     degrees <- summary(g ~ event, drop=FALSE)
#     if(any(degrees==0)){
#      cat(paste("Warning: There are no actors for the event", 
#          dnames[degrees==0],"\n"))
#      dropterms <- dnames[degrees==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     if(any(degrees==nactors)){
#      cat(paste("Warning: The event", 
#          dnames[degrees],"is common to all actors\n"))
#      dropterms <- dnames[degrees==nactors]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     d <- d[degrees!=nactors & degrees!=0] 
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
##
#    if(is.null(xnames)){
#     dnames <- paste("event",d-nactors,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="event", soname="statnet",
#                                          inputs=c(0, ld, ld+1, nactors, d))
##   model$coef.names<-c(model$coef.names,dnames)
#    model$coef.names<-c(model$coef.names,paste("event.",dnames,sep=""))
#    model
#}
#
#InitErgm.actor<-function(g, model, drop=TRUE, ...)
#{
#    nevents <- is.bipartite(g)
#    if (!nevents)
#      stop("The actor term is for bipartite graphs.",
#           call.=FALSE)
#    if (nargs()!=3)
#        stop(paste("actor model term expected zero argument, got ", 
#            nargs()-2, sep=""), call.=FALSE)
##
##   Check for degeneracy
## 
#    nactors <- get.network.attribute(g,"bipartite")
#    nevents <- network.size(g)-nactors
#    d <- 2:nactors
#    xnames <- network.vertex.names(g)
#    if(is.null(xnames)){
#     dnames <- paste("actor",d,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    if(drop){
#     degrees <- summary(g ~ actor, drop=FALSE)
#     if(any(degrees==0)){
#      cat(paste("Warning: There are no events for the actor", 
#          dnames[degrees==0],"\n"))
#      dropterms <- dnames[degrees==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     if(any(degrees==nevents)){
#      cat(paste("Warning: The actor", 
#          dnames[degrees==nevents],"participated in all events.\n"))
#      dropterms <- dnames[degrees==nevents]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     d <- d[degrees!=nevents & degrees!=0] 
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
##
#    if(is.null(xnames)){
#     dnames <- paste("actor",d,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="actor", soname="statnet",
#                                          inputs=c(0, ld, ld+1, nactors, d))
##   model$coef.names<-c(model$coef.names,dnames)
#    model$coef.names<-c(model$coef.names,paste("actor.",dnames,sep=""))
#    model
#}
#
InitErgm.gwevent706<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("biduration", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwevent", arglist,
    varnames = c("alpha"),
    vartypes = c("numeric"),
    defaultvalues = list(0.5),
    required = c(FALSE))
  attach(a)
  alpha<-a$alpha
  nactors <- get.network.attribute(nw,"bipartite")
  nevents <- network.size(nw)-nactors
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="gwevent",
                                    soname="statnet",
                                    inputs=c(0, 1, 2, nactors, alpha))
  m$coef.names<-c(m$coef.names,"gwevent")
  m
}

InitErgm.gwactor706<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("biduration", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwactor", arglist,
    varnames = c("alpha"),
    vartypes = c("numeric"),
    defaultvalues = list(0.5),
    required = c(FALSE))
  attach(a)
  alpha<-a$alpha
  nactors <- get.network.attribute(nw,"bipartite")
  nevents <- network.size(nw)-nactors
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="gwactor",
                                    soname="statnet",
                                    inputs=c(0, 1, 2, nactors, alpha))
  m$coef.names<-c(m$coef.names,"gwactor")
  m
}

InitErgm.r0a<-function(nw, m, arglist, ...) {
  ergm.checkdirected("r0a", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("r0a", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="r0a", soname="statnet",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"r0a")
  m
}

InitErgm.r0e<-function(nw, m, arglist, ...) {
  ergm.checkdirected("r0e", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("r0e", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="r0e", soname="statnet",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"r0e")
  m
}
#InitErgm.esa<-function(g, model, d, drop=TRUE, ...)
#{
#    if (nargs()!=4)
#        stop(paste("esa() model term expected 1 argument, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#    nevents <- is.bipartite(g)
#    if (!nevents)
#      stop("The esa term is for bipartite graphs.",
#           call.=FALSE)
#    nactors <- get.network.attribute(g,"bipartite")
#    if (is.directed(g))
#      stop("the esa() term is not allowed with a directed graph",
#           call.=FALSE)
##
##   Check for degeneracy
## 
#    if(drop){
#     mesa <- paste("c(",paste(d,collapse=","),")",sep="")
#     mesa <- summary(
#       as.formula(paste('g ~ esa(',mesa,')',sep="")),
#       drop=FALSE)
#     if(any(mesa==0)){
#      cat(paste("Warning: There are no esa", d[mesa==0],"dyads.\n"))
#      dropterms <- paste("esa", d[mesa==0],sep="")
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      d <- d[mesa!=0] 
#     }
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="esa", soname="statnet",
#                                          inputs=c(0, ld, ld+1, nactors, d))
#    model$coef.names<-c(model$coef.names,paste("esa",d,sep=""))
#    model
#}
#
#InitErgm.ase<-function(g, model, d, drop=TRUE, ...)
#{
#    if (nargs()!=4)
#        stop(paste("ase() model term expected 1 argument, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#    nevents <- is.bipartite(g)
#    if (!nevents)
#      stop("The ase term is for bipartite graphs.",
#           call.=FALSE)
#    nactors <- get.network.attribute(g,"bipartite")
#    nevents <- network.size(g)-nactors
#    if (is.directed(g))
#      stop("the ase() term is not allowed with a directed graph",
#           call.=FALSE)
##
##   Check for degeneracy
## 
#    if(drop){
#     mase <- paste("c(",paste(d,collapse=","),")",sep="")
#     mase <- summary(
#       as.formula(paste('g ~ ase(',mase,')',sep="")),
#       drop=FALSE)
#     if(any(mase==0)){
#      cat(paste("Warning: There are no ase", d[mase==0],"dyads.\n"))
#      dropterms <- paste("ase", d[mase==0],sep="")
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      d <- d[mase!=0] 
#     }
#    }
#    ld<-length(d)
#    if(ld==0){return(m)}
#    termnumber<-1+length(m$terms)
##  No covariates here, so input component 1 is arbitrary
#    m$terms[[termnumber]] <- list(name="ase", soname="statnet",
#                                          inputs=c(0, ld, ld+1, nactors, d))
#    m$coef.names<-c(m$coef.names,paste("ase",d,sep=""))
#    m
#}
#
#InitErgm.bimixall<-function (nw, m, arglist, drop=TRUE, ...) {
#  a <- ergm.checkargs("bimix", arglist,
#    varnames = c("attrname","contrast"),
#    vartypes = c("character","logical"),
#    defaultvalues = list(NULL,FALSE),
#    required = c(TRUE, FALSE))
#  attach(a)
#  attrname<-a$attrname
#  nactors <- get.network.attribute(nw,"bipartite")
#  nevents <- network.size(nw)-nactors
##
#  nodecov <- get.node.attr(nw, attrname, "bimixall")
#  mixmat <- mixingmatrix(nw,attrname)
#  mixmat <- mixmat[-nrow(mixmat),-nrow(mixmat)]
#  u <- cbind(as.vector(row(mixmat)), as.vector(col(mixmat)))
#  if(any(is.na(nodecov))){u<-rbind(u,NA)}
##
##   Recode to numeric if necessary
##
#  nodecov <- match(nodecov,sort(unique(nodecov)))
#  if (length(nodecov)==1)
#        stop ("Argument to bimixall() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#  if(drop){
#   mesa <- summary(
#     as.formula(paste('nw ~ bimixall("',attrname,'")',sep="")),
#     drop=FALSE)
#   if(any(mesa==0)){
#    cat(paste("Warning: There are no counts for the bimixall", 
#      mesa[mesa==0],"\n"))
#    dropterms <- mesa[mesa==0]
#    cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#   }
#   u <- u[mesa>0,] 
#  }
#  if(contrast){
#   cat(paste("Using contrasted version: bimixall.",attrname,".",
#     paste(c(u[1,2],u[1,1]),collapse="")," has been dropped.\n",
#     sep=""))
#   u <- u[-1,] 
#  }
#  termnumber<-1+length(m$terms)  
#  #  Number of input parameters before covariates equals number of
#  #  unique elements in nodecov, namely length(u), so that's what
#  #  input component 1 equals
#  m$terms[[termnumber]] <- list(name="bimix", soname="statnet",
#      inputs=c(nrow(u), nrow(u), length(nodecov)+length(u),
#        u[,1], u[,2],nodecov),
#      dependence=FALSE)
#  m$coef.names<-c(m$coef.names,paste("bimixall",
#      attrname, apply(cbind(u[,2],u[,1]),1,paste,collapse=""), sep="."))
#  m
#}
#
#InitErgm.bimixconddeg<-function (g, model, attrname, drop=TRUE, ...) {
#    if (!is.bipartite(g))
#      stop("The bimixconddeg term is for bipartite graphs.",
#           call.=FALSE)
#    if (!(nargs() %in% 4:5))
#        stop(paste("bimixconddeg() expected 1 arguments, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#
#    nactors <- get.network.attribute(g,"bipartite")
#    nevents <- network.size(g)-nactors
##
#    nodecov <- get.node.attr(g, attrname, "bimixconddeg")
#    mixmat <- mixingmatrix(g,attrname)
#    mixmat <- mixmat[-nrow(mixmat),-nrow(mixmat)]
#    mixmat <- mixmat[-1,-1,drop=FALSE]
#    u <- cbind(as.vector(row(mixmat)+1), as.vector(col(mixmat)+1))
#    if(any(is.na(nodecov))){u<-rbind(u,NA)}
##
##   Recode to numeric if necessary
##
#    nodecov <- match(nodecov,sort(unique(nodecov)))
#    if (length(nodecov)==1)
#        stop ("Argument to bimixconddeg() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#    if(drop){
#     mesa <- summary(
#       as.formula(paste('g ~ bimixconddeg("',attrname,'")',sep="")),
#       drop=FALSE)
#     if(any(mesa==0)){
#      cat(paste("Warning: There are no counts for the bimixconddeg", 
#          mesa[mesa==0],"\n"))
#      dropterms <- mesa[mesa==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     u <- u[mesa>0,] 
#    }
##   u <- u[-1,] 
#    termnumber<-1+length(model$terms)  
#    #  Number of input parameters before covariates equals number of
#    #  unique elements in nodecov, namely length(u), so that's what
#    #  input component 1 equals
#    model$terms[[termnumber]] <- list(name="bimix", soname="statnet",
#        inputs=c(nrow(u), nrow(u), length(nodecov)+length(u),
#          u[,1], u[,2],nodecov))
#    model$coef.names<-c(model$coef.names, paste("bimixconddeg",
#        attrname, apply(cbind(u[,2],u[,1]),1,paste,collapse=""), sep="."))
#    model
#}
#
#InitErgm.bimix<-function (g, model, attrname, drop=TRUE, ...) {
#    if (!is.bipartite(g))
#      stop("The bimix term is for bipartite graphs.",
#           call.=FALSE)
#    if (!(nargs() %in% 4:5))
#        stop(paste("bimix() expected 1 arguments, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#
#    nactors <- get.network.attribute(g,"bipartite")
#    nevents <- network.size(g)-nactors
##
#    nodecov <- get.node.attr(g, attrname, "bimix")
#    mixmat <- mixingmatrix(g,attrname)
#    mixmat <- mixmat[-nrow(mixmat),-nrow(mixmat)]
##   mixmat <- mixmat[-1,-1,drop=FALSE]
#    u <- cbind(as.vector(row(mixmat)), as.vector(col(mixmat)))
#    if(any(is.na(nodecov))){u<-rbind(u,NA)}
##
##   Drop the [1,1] cell
##
#    u <- u[-1,] 
##
##   Recode to numeric if necessary
##
#    nodecov <- match(nodecov,sort(unique(nodecov)))
#    if (length(nodecov)==1)
#        stop ("Argument to bimix() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#    if(drop){
#     mesa <- summary(
#       as.formula(paste('g ~ bimix("',attrname,'")',sep="")),
#       drop=FALSE)
#     if(any(mesa==0)){
#      cat(paste("Warning: There are no counts for the bimix", 
#          mesa[mesa==0],"\n"))
#      dropterms <- mesa[mesa==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     u <- u[mesa>0,] 
#    }
#    termnumber<-1+length(model$terms)  
#    #  Number of input parameters before covariates equals number of
#    #  unique elements in nodecov, namely length(u), so that's what
#    #  input component 1 equals
#    model$terms[[termnumber]] <- list(name="bimix", soname="statnet",
#        inputs=c(nrow(u), nrow(u), length(nodecov)+length(u),
#          u[,1], u[,2],nodecov))
#    model$coef.names<-c(model$coef.names, paste("bimix",
#        attrname, apply(cbind(u[,2],u[,1]),1,paste,collapse=""), sep="."))
#    model
#}
##
## End untranslated Inits
##

InitErgm.eventfactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("eventfactor", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("eventfactor", arglist,
    varnames = c("attrname","contrast"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL, TRUE),
    required = c(TRUE,FALSE))
  attach(a)
  attrname<-a$attrname
  nodecov <- get.node.attr(nw, attrname, "eventfactor")
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop){
    if (!is.directed(nw)){
      nfc <- tapply(tabulate(as.matrix.network.edgelist(nw),
                             nbins=network.size(nw)),
                    nodecov,sum)
    }else{
      nfc <- tapply(tabulate(as.matrix.network.edgelist(nw),
                             nbins=network.size(nw)),
                    nodecov,sum)
    }
    if(any(nfc==0)){
      dropterms <- paste(paste("eventfactor",attrname,sep="."),u[nfc==0],sep="")
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  if(contrast){
   ui <- ui[-1]
   u <- u[-1]
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to eventfactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="eventfactor", soname="statnet",
                                inputs=c(lu, lu, lu+length(nodecov),
                                         ui, nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("eventfactor",
                                      attrname, paste(u), sep="."))
  m
}

InitErgm.actorfactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("actorfactor", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("actorfactor", arglist,
    varnames = c("attrname","contrast"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,TRUE),
    required = c(TRUE,FALSE))
  attach(a)
  attrname<-a$attrname
  nodecov <- get.node.attr(nw, attrname, "actorfactor")
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop){
    if (!is.directed(nw)){
      nfc <- tapply(tabulate(as.matrix.network.edgelist(nw),
                             nbins=network.size(nw)),
                    nodecov,sum)
    }else{
      nfc <- tapply(tabulate(as.matrix.network.edgelist(nw),
                             nbins=network.size(nw)),
                    nodecov,sum)
    }
    if(any(nfc==0)){
      dropterms <- paste(paste("actorfactor",attrname,sep="."),u[nfc==0],sep="")
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  if(contrast){
   ui <- ui[-1]
   u <- u[-1]
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to actorfactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="actorfactor", soname="statnet",
                                inputs=c(lu, lu, lu+length(nodecov),
                                         ui, nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("actorfactor",
                                      attrname, paste(u), sep="."))
  m
}
InitErgm.biduration<-function (nw, m, arglist, ...) {
  ergm.checkdirected("biduration", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("biduration", arglist,
    varnames = c("form", "dissolve", "x"),
    vartypes = c("matrixnetwork", "matrixnetwork", "matrixnetwork"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(TRUE, TRUE, FALSE))
  attach(a)
  x<-a$x;form<-a$form;dissolve<-a$dissolve
# m$coef.names<-c(m$coef.names, paste("biduration.",x,sep=""))
  #Coerce x to an adjacency matrix
  if(is.null(x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")
    x<-"self"
  }else{
    if(is.network(x)){
      xm<-as.matrix.network(x,matrix.type="edgelist")
      x<-paste(quote(x))
    }else if(is.character(x)){
      xm<-get.network.attribute(nw,x)
      xm<-as.matrix.network(xm,matrix.type="edgelist")
    }else{
      xm<-as.matrix(x)
      x<-paste(quote(x))
    }
  }
  #Check for symmetry
  if (is.null(xm) || ncol(xm)!=2){
    stop("biduration requires the edgelist of the base network")
  }
  nactors <- get.network.attribute(nw,"bipartite")
  nevents <- network.size(nw)-nactors
  #Coerce form to an adjacency matrix
  if(is.network(form)){
    formm<-as.matrix.network(form,matrix.type="adjacency")
    form<-paste(quote(form))
  }else if(is.character(form)){
    formm<-get.network.attribute(nw,form)
  }else{
    formm<-as.matrix(form)
    form<-paste(quote(form))
  }
  #Check for matrix
  if (is.null(formm) || dim(formm)!=c(nactors, nevents)){
    stop("biduration requires a matrix of formation rates")
  }
  #Coerce dissolve to an adjacency matrix
  if(is.network(dissolve)){
    dissolvem<-as.matrix.network(dissolve,matrix.type="adjacency")
    dissolve<-paste(quote(dissolve))
  }else if(is.character(dissolve)){
    dissolvem<-get.network.attribute(nw,dissolve)
  }else{
    dissolvem<-as.matrix(dissolve)
    dissolve<-paste(quote(dissolve))
  }
  #Check for matrix
  if (is.null(dissolvem) || dim(dissolvem)!=c(nactors, nevents)){
    stop("biduration requires a matrix of dissolution rates")
  }
  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_biduration ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "biduration", soname="statnet",
                                inputs = c(1, 1, 
                        2 + 2*nrow(xm)+nrow(formm)*ncol(formm),
                        nrow(xm), nrow(formm), as.double(c(xm, formm)))
                               )
  #Update the coefficient name list, adding biduration.nameofx
  model$coef.names<-c(model$coef.names, paste(c("duration","formation"),x,sep="."))
  #Return the updated model list
  model
}

InitErgm.biendure<-function (g, model, drop=TRUE, ...) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 3))
    stop(paste("biendure() model term expected 0 arguments, got ", 
                                nargs() - 3, sep = ""), call. = FALSE)
  xm<-as.matrix.network(g,matrix.type="edgelist")
  nactors <- get.network.attribute(g,"bipartite")
  nevents <- network.size(g)-nactors
  #Check for symmetry
  if (is.null(xm) || ncol(xm)!=2){
    stop("biendure requires the edgelist of the base graph")
  }
  #Update the term number
  termnumber <- 1 + length(model$terms)
  #Update the terms list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_biendure ignores the value of inp->attrib).
  model$terms[[termnumber]] <- list(name = "endure",
                                        soname="statnet",
             inputs = c(1, 1, 
                        2 + 2*nrow(xm),
                        nrow(xm), nactors, as.double(c(xm))
                       ) )
  #Update the coefficient name list, adding dyadcov.nameofx
  model$coef.names<-c(model$coef.names, paste("endure"))
  #Return the updated model list
  model
}

InitErgm.bichange<-function (g, model, form=NULL, x=NULL, drop=TRUE, ...) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 3:5))
    stop(paste("bichange() model term expected 0 or 1 arguments, got ", 
                                nargs() - 3, sep = ""), call. = FALSE)
  if (nargs()==4){drop <- x; x <- NULL}
  #Coerce x to an adjacency matrix
  if(is.null(x)){
   xm<-as.matrix.network(g,matrix.type="edgelist")
   x<-"self"
  }else{
   if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist")
    x<-paste(quote(x))
   }else if(is.character(x)){
#   xm<-as.matrix.network(g,matrix.type="edgelist",x)
    xm<-get.network.attribute(g,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
   }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
   }
  }
  nactors <- get.network.attribute(g,"bipartite")
  nevents <- network.size(g)-nactors
  #Check for symmetry
  if (is.null(xm) || ncol(xm)!=2){
    stop("bichange requires the edgelist of the base graph")
  }
  #Coerce form to an adjacency matrix
  if(is.network(form)){
    formm<-as.matrix.network(form,matrix.type="adjacency")
    form<-paste(quote(form))
  }else if(is.character(form)){
    formm<- g %n% form
  }else{
   if(is.null(x)){
    formm <- matrix(1, nrow=nactors, ncol=nevents)
   }else{
    formm<-as.matrix(form)
    form<-paste(quote(form))
   }
  }
  #Check for matrix
  if (is.null(formm) || dim(formm)!=c(nactors, nevents)){
    stop("bichange requires a matrix of formation rates")
  }
  #Update the term number
  termnumber <- 1 + length(model$terms)
  #Update the terms list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_bichange ignores the value of inp->attrib).
  model$terms[[termnumber]] <- list(name = "bichange",
                                        soname="statnet",
             inputs = c(1, 1, 
                        2 + 2*nrow(xm)+nrow(formm)*ncol(formm),
                        nrow(xm), nrow(formm), as.double(c(xm, formm))
                       ) )
  #Update the coefficient name list, adding dyadcov.nameofx
  model$coef.names<-c(model$coef.names, paste(c("change"),x,sep="."))
  #Return the updated model list
  model
}

InitErgm.adegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("adegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("adegree", arglist,
    varnames = c("d", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  d<-a$d
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "adegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to adegree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      adegreeattr <- summary(
       as.formula(paste('nw ~ adegree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(adegreeattr)){
        dropterms <- paste("adeg", du[1,adegreeattr], ".", attrname,
                           u[du[2,adegreeattr]], sep="")
        cat("Warning: These adegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
        du <- matrix(du[,!adegreeattr], nrow=2)
      }
    }
  }else{
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      madegree <- paste("c(",paste(d,collapse=","),")",sep="")
      madegree <- summary(
       as.formula(paste('nw ~ adegree(',madegree,')',sep="")),
       drop=FALSE) == 0
      if(any(madegree)){
       cat(paste("Warning: There are no order", d[madegree],"adegrees.\n"))
       dropterms <- paste("adegree", d[madegree],sep="")
       cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
       d <- d[!madegree] 
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="adegree_by_attr", soname="statnet",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_adegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("adeg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="adegree", soname="statnet",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("adegree",d,sep=""))
  }
  m
}

InitErgm.edegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("edegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("edegree", arglist,
    varnames = c("d", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  d<-a$d
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "edegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to edegree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      edegreeattr <- summary(
       as.formula(paste('nw ~ edegree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(edegreeattr)){
        dropterms <- paste("edeg", du[1,edegreeattr], ".", attrname,
                           u[du[2,edegreeattr]], sep="")
        cat("Warning: These edegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
        du <- matrix(du[,!edegreeattr], nrow=2)
      }
    }
  }else{
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      medegree <- paste("c(",paste(d,collapse=","),")",sep="")
      medegree <- summary(
       as.formula(paste('nw ~ edegree(',medegree,')',sep="")),
       drop=FALSE) == 0
      if(any(medegree)){
       cat(paste("Warning: There are no order", d[medegree],"edegrees.\n"))
       dropterms <- paste("edegree", d[medegree],sep="")
       cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
       d <- d[!medegree] 
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="edegree_by_attr", soname="statnet",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_edegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("edeg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="edegree", soname="statnet",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("edegree",d,sep=""))
  }
  m
}


InitErgm.gwadegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkbipartite("gwadegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwadegree", arglist,
    varnames = c("decay", "attrname", "fixed"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(0, NULL, TRUE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  nactors <- get.network.attribute(nw,"bipartite")
  d <- 1:(network.size(nw) - nactors)
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwadegree term is not yet able to handle a",
           "nonfixed decay term.") # with an attribute.")
#    }
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^
                  {i-1}*(1+i-exp(-x[2])))
            )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="adegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwadegree=NULL,
                                    gwadegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwadegree#",d,sep=""))
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "gwadegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to gwadegree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(nrow(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="gwadegree_by_attr", soname="statnet",
                                  inputs=c(0, nrow(du), 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    # See comment in d_gwadegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("gwadeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwadegree", soname="statnet",
                                       inputs=c(0, 1, 1, decay),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("gwadeg",decay,sep=""))
  }
  m
}

InitErgm.gwedegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkbipartite("gwedegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwedegree", arglist,
    varnames = c("decay", "attrname", "fixed"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(0, NULL, TRUE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  nactors <- get.network.attribute(nw,"bipartite")
  d <- 1:nactors
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwedegree term is not yet able to handle a",
           "nonfixed decay term.") # with an attribute.")
#    }
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^
                  {i-1}*(1+i-exp(-x[2])))
            )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="adegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwedegree=NULL,
                                    gwedegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwedegree#",d,sep=""))
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "gwedegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to gwedegree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(nrow(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="gwedegree_by_attr", soname="statnet",
                                  inputs=c(0, nrow(du), 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    # See comment in d_gwedegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("gwedeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwedegree", soname="statnet",
                                       inputs=c(0, 1, 1, decay),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("gwedeg",decay,sep=""))
  }
  m
}

InitErgm.biduration<-function (nw, m, arglist, ...) {
  ergm.checkdirected("biduration", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("biduration", arglist,
    varnames = c("dissolve", "x"),
    vartypes = c("matrixnetwork", "matrixnetwork"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  x<-a$x;dissolve<-a$dissolve
# m$coef.names<-c(m$coef.names, paste("biduration.",x,sep=""))
  #Coerce x to an adjacency matrix
  if(is.null(x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")
    x<-"self"
  }else{
    if(is.network(x)){
      xm<-as.matrix.network(x,matrix.type="edgelist")
      x<-paste(quote(x))
    }else if(is.character(x)){
      xm<-get.network.attribute(nw,x)
      xm<-as.matrix.network(xm,matrix.type="edgelist")
    }else{
      xm<-as.matrix(x)
      x<-paste(quote(x))
    }
  }
  #Check for symmetry
  if (is.null(xm) || ncol(xm)!=2){
    stop("biduration requires the edgelist of the base network")
  }
  nactors <- get.network.attribute(nw,"bipartite")
  nevents <- network.size(nw)-nactors
  #Coerce dissolve to an adjacency matrix
  if(is.network(dissolve)){
    dissolvem<-as.matrix.network(dissolve,matrix.type="adjacency")
    dissolve<-paste(quote(dissolve))
  }else if(is.character(dissolve)){
    dissolvem<-get.network.attribute(nw,dissolve)
  }else{
    dissolvem<-as.matrix(dissolve)
    dissolve<-paste(quote(dissolve))
  }
  #Check for matrix
  if (is.null(dissolvem) || dim(dissolvem)!=c(nactors, nevents)){
    stop("biduration requires a matrix of dissolution rates")
  }
  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_biduration ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "biduration", soname="statnet",
                                inputs = c(1, 1, 
                        2 + 2*nrow(xm)+nrow(dissolvem)*ncol(dissolvem),
                        nrow(xm), nrow(dissolvem), as.double(c(xm, dissolvem)))
                               )
  #Update the coefficient name list, adding biduration.nameofx
# m$coef.names<-c(m$coef.names, paste(c("duration","formation"),x,sep="."))
  m$coef.names<-c(m$coef.names, paste(c("duration"),x,sep="."))
  m
}

InitErgm.hammingmix<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingmix", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingmix", arglist=arglist,
    varnames = c("attrname","x","omitwhich"),
    vartypes = c("character","matrixnetwork","numeric"),
    defaultvalues = list(NULL,nw,0),
    required = c(TRUE,FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  x<-a$x
  omitwhich<-a$omitwhich
  contrast<-a$contrast
  drop<-a$drop
  contrast<-TRUE
  drop<-TRUE
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist",attrname)
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    stop("hammingmix() requires an edgelist")
  }
    nodecov <- get.node.attr(nw, attrname, "hammingmix")
    mixmat <- mixingmatrix(nw,attrname)
    mixmat <- mixmat[-nrow(mixmat),-nrow(mixmat)]
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
#
#   Recode to numeric if necessary
#
    namescov <- sort(unique(nodecov))
    nodecov <- match(nodecov,namescov)
    if (length(nodecov)==1)
        stop ("Argument to hammingmix() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#    if(drop){
#     ematch <- mixmat[u]
#     mu <- ematch==0
#     mu[is.na(mu)] <- FALSE
#     if(any(mu)){
#      dropterms <- paste(paste("hammingmix",attrname,sep="."),
#        apply(u,1,paste,collapse="")[mu],sep="")
#      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      u <- u[!mu,]
#     }
#    }
  if(contrast){
   u <- u[-1,]
  }
  if(all(omitwhich!=0)){
   u <- u[-omitwhich,]
  }
  termnumber<-1+length(m$terms)
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui), so that's what
  #  input component 1 equals
  m$terms[[termnumber]] <- list(name="hammingmix", soname="statnet",
    inputs=c(nrow(u), nrow(u), nrow(xm)*2+length(nodecov)+length(u)+1,
            nrow(xm),as.integer(xm), u[,1], u[,2],nodecov),
            dependence=FALSE)
  m$coef.names<-c(m$coef.names,
       paste("hammingmix",attrname, apply(matrix(namescov[u],ncol=2),1,paste,collapse="."), sep="."))
  m
}

InitErgm.hammingfixmix<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingfixmix", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingfixmix", arglist=arglist,
    varnames = c("attrname","x","omitwhich"),
    vartypes = c("character","matrixnetwork","numeric"),
    defaultvalues = list(NULL,nw,0),
    required = c(TRUE,FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  x<-a$x
  omitwhich<-a$omitwhich
  contrast<-a$contrast
  drop<-a$drop
  contrast<-TRUE
  drop<-TRUE
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist",attrname)
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    stop("hammingfixmix() requires an edgelist")
  }
    nodecov <- get.node.attr(nw, attrname, "hammingfixmix")
    mixmat <- mixingmatrix(nw,attrname)
    mixmat <- mixmat[-nrow(mixmat),-nrow(mixmat)]
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
#
#   Recode to numeric if necessary
#
    namescov <- sort(unique(nodecov))
    nodecov <- match(nodecov,namescov)
    if (length(nodecov)==1)
        stop ("Argument to hammingfixmix() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#    if(drop){
#     ematch <- mixmat[u]
#     mu <- ematch==0
#     mu[is.na(mu)] <- FALSE
#     if(any(mu)){
#      dropterms <- paste(paste("hammingfixmix",attrname,sep="."),
#        apply(u,1,paste,collapse="")[mu],sep="")
#      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      u <- u[!mu,]
#     }
#    }
  if(contrast){
   u <- u[-1,]
  }
  if(all(omitwhich!=0)){
   u <- u[-omitwhich,]
  }
  termnumber<-1+length(m$terms)
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui), so that's what
  #  input component 1 equals
  m$terms[[termnumber]] <- list(name="hammingfixmix", soname="statnet",
    inputs=c(1, 1, nrow(xm)*2+length(nodecov)+1,
            nrow(xm),as.integer(xm), nodecov),
            dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("hammingfixmix",attrname, sep="."))
  m
}

InitErgm.hammingdyadcov<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingdyadcov", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingdyadcov", arglist=arglist,
    varnames = c("cov","x","covattrname"),
    vartypes = c("matrixnetwork","matrixnetwork","character"),
    defaultvalues = list(NULL,nw,NULL),
    required = c(TRUE,FALSE,FALSE))
  attach(a)
  covattrname<-a$covattrname
  x<-a$x
  cov<-a$cov
#
# Extract dyadic covariate
#
  if(is.network(cov)){
    covm<-as.matrix.network(cov,matrix.type="adjacency",covattrname)
    cov<-paste(quote(cov))
  }else if(is.character(cov)){
    covm<-get.network.attribute(nw,cov)
    covm<-as.matrix.network(covm,matrix.type="adjacency")
  }else{
    covm<-as.matrix(cov)
    cov<-paste(quote(cov))
  }
  if (is.null(covm) || !is.matrix(covm) || nrow(covm)!=get.network.attribute(nw,"bipartite")){
    stop("hammingdyadcov() requires an proper dyadic covariate")
  }
#
# Extract reference network as an edgelist
#
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist")
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else if(is.null(x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    stop("hammingdyadcov() requires an proper network as its reference")
  }
##
##   Check for degeneracy
##
#    if(drop){
#     ematch <- mixmat[u]
#     mu <- ematch==0
#     mu[is.na(mu)] <- FALSE
#     if(any(mu)){
#      dropterms <- paste(paste("hammingdyadcov",attrname,sep="."),
#        apply(u,1,paste,collapse="")[mu],sep="")
#      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      u <- u[!mu,]
#     }
#    }
  termnumber<-1+length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "hammingdyadcov", soname="statnet",
                                 inputs = c(1, 1,
                                   1+2*nrow(xm)+nrow(covm)*ncol(covm),
                                   nrow(xm), as.integer(xm),
                                   as.double(covm)),
                                 dependence=TRUE)
   if(!is.null(covattrname)){
     cn<-paste("hammingdyadcov", as.character(sys.call(0)[[4]][2]),
               as.character(covattrname), sep = ".")
   }else{
     cn<-paste("hammingdyadcov", as.character(sys.call(0)[[4]][2]), sep = ".")
   }
   m$coef.names <- c(m$coef.names, cn)
   m
}
