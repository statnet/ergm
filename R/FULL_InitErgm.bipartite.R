#  See InitErgm.R for a general explanation 
#  of InitErgm functions

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
#    model$terms[[termnumber]] <- list(name="event", soname="ergm",
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
#    model$terms[[termnumber]] <- list(name="actor", soname="ergm",
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
                                    soname="ergm",
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
                                    soname="ergm",
                                    inputs=c(0, 1, 2, nactors, alpha))
  m$coef.names<-c(m$coef.names,"gwactor")
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
#    model$terms[[termnumber]] <- list(name="esa", soname="ergm",
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
#    m$terms[[termnumber]] <- list(name="ase", soname="ergm",
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
#  mixmat <- mixingmatrix(nw,attrname)$mat
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
#  m$terms[[termnumber]] <- list(name="bimix", soname="ergm",
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
#    mixmat <- mixingmatrix(g,attrname)$mat
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
#    model$terms[[termnumber]] <- list(name="bimix", soname="ergm",
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
#    mixmat <- mixingmatrix(g,attrname)$mat
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
#    model$terms[[termnumber]] <- list(name="bimix", soname="ergm",
#        inputs=c(nrow(u), nrow(u), length(nodecov)+length(u),
#          u[,1], u[,2],nodecov))
#    model$coef.names<-c(model$coef.names, paste("bimix",
#        attrname, apply(cbind(u[,2],u[,1]),1,paste,collapse=""), sep="."))
#    model
#}
##
## End untranslated Inits
##

## Commented out because there is another version of InitErgm.biduration below
#InitErgm.biduration<-function (nw, m, arglist, ...) {
#  ergm.checkdirected("biduration", is.bipartite(nw), requirement=TRUE)
#  a <- ergm.checkargs("biduration", arglist,
#    varnames = c("form", "dissolve", "x"),
#    vartypes = c("matrixnetwork", "matrixnetwork", "matrixnetwork"),
#    defaultvalues = list(NULL, NULL, NULL),
#    required = c(TRUE, TRUE, FALSE))
#  attach(a)
#  x<-a$x;form<-a$form;dissolve<-a$dissolve
## m$coef.names<-c(m$coef.names, paste("biduration.",x,sep=""))
#  #Coerce x to an adjacency matrix
#  if(is.null(x)){
#    xm<-as.matrix.network(nw,matrix.type="edgelist")
#    x<-"self"
#  }else{
#    if(is.network(x)){
#      xm<-as.matrix.network(x,matrix.type="edgelist")
#      x<-paste(quote(x))
#    }else if(is.character(x)){
#      xm<-get.network.attribute(nw,x)
#      xm<-as.matrix.network(xm,matrix.type="edgelist")
#    }else{
#      xm<-as.matrix(x)
#      x<-paste(quote(x))
#    }
#  }
#  #Check for symmetry
#  if (is.null(xm) || ncol(xm)!=2){
#    stop("biduration requires the edgelist of the base network")
#  }
#  nactors <- get.network.attribute(nw,"bipartite")
#  nevents <- network.size(nw)-nactors
#  #Coerce form to an adjacency matrix
#  if(is.network(form)){
#    formm<-as.matrix.network(form,matrix.type="adjacency")
#    form<-paste(quote(form))
#  }else if(is.character(form)){
#    formm<-get.network.attribute(nw,form)
#  }else{
#    formm<-as.matrix(form)
#    form<-paste(quote(form))
#  }
#  #Check for matrix
#  if (is.null(formm) || dim(formm)!=c(nactors, nevents)){
#    stop("biduration requires a matrix of formation rates")
#  }
#  #Coerce dissolve to an adjacency matrix
#  if(is.network(dissolve)){
#    dissolvem<-as.matrix.network(dissolve,matrix.type="adjacency")
#    dissolve<-paste(quote(dissolve))
#  }else if(is.character(dissolve)){
#    dissolvem<-get.network.attribute(nw,dissolve)
#  }else{
#    dissolvem<-as.matrix(dissolve)
#    dissolve<-paste(quote(dissolve))
#  }
#  #Check for matrix
#  if (is.null(dissolvem) || dim(dissolvem)!=c(nactors, nevents)){
#    stop("biduration requires a matrix of dissolution rates")
#  }
#  termnumber <- 1 + length(m$terms)
## There is 1 input parameter before the covariate vector, so input
## component 1 is set to 1 (although in this case, input component 1
## is actually arbitrary since d_biduration ignores the value of inp->attrib).
#  m$terms[[termnumber]] <- list(name = "biduration", soname="ergm",
#                                inputs = c(1, 1, 
#                        2 + 2*nrow(xm)+nrow(formm)*ncol(formm),
#                        nrow(xm), nrow(formm), as.double(c(xm, formm)))
#                               )
#  #Update the coefficient name list, adding biduration.nameofx
#  model$coef.names<-c(model$coef.names, paste(c("duration","formation"),x,sep="."))
#  #Return the updated model list
#  model
#}

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
                                        soname="ergm",
             inputs = c(1, 1, 
                        2 + 2*nrow(xm),
                        nrow(xm), nactors, as.double(c(xm))
                       ) )
  #Update the coefficient name list, adding dyadcov.nameofx
  model$coef.names<-c(model$coef.names, paste("endure"))
  #Return the updated model list
  model
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
  m$terms[[termnumber]] <- list(name = "biduration", soname="ergm",
                                inputs = c(1, 1, 
                        2 + 2*nrow(xm)+nrow(dissolvem)*ncol(dissolvem),
                        nrow(xm), nrow(dissolvem), as.double(c(xm, dissolvem)))
                               )
  #Update the coefficient name list, adding biduration.nameofx
# m$coef.names<-c(m$coef.names, paste(c("duration","formation"),x,sep="."))
  m$coef.names<-c(m$coef.names, paste(c("duration"),x,sep="."))
  m
}


