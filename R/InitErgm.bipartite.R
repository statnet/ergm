#  File ergm/R/InitErgm.bipartite.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
#  See InitErgm.R for a general explanation 
#  of InitErgm functions

###################################### InitErgm TERMS:  A
##########################################################
#InitErgm.b1<-function(g, model, drop=TRUE, ...)
#{
#    nb2 <- is.bipartite(g)
#    if (!nb2)
#      stop("The b1 term is for bipartite graphs.",
#           call.=FALSE)
#    if (nargs()!=3)
#        stop(paste("b1 model term expected zero argument, got ", 
#            nargs()-2, sep=""), call.=FALSE)
##
##   Check for degeneracy
## 
#    nb1 <- get.network.attribute(g,"bipartite")
#    nb2 <- network.size(g)-nb1
#    d <- 2:nb1
#    xnames <- network.vertex.names(g)
#    if(is.null(xnames)){
#     dnames <- paste("b1",d,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    if(drop){
#     degrees <- summary(g ~ b1, drop=FALSE)
#     if(any(degrees==0)){
#      cat(paste("Warning: There are no b2s for the b1", 
#          dnames[degrees==0],"\n"))
#      dropterms <- dnames[degrees==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     if(any(degrees==nb2)){
#      cat(paste("Warning: The b1", 
#          dnames[degrees==nb2],"participated in all b2s.\n"))
#      dropterms <- dnames[degrees==nb2]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     d <- d[degrees!=nb2 & degrees!=0] 
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
##
#    if(is.null(xnames)){
#     dnames <- paste("b1",d,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="b1", soname="ergm",
#                                          inputs=c(0, ld, ld+1, nb1, d))
##   model$coef.names<-c(model$coef.names,dnames)
#    model$coef.names<-c(model$coef.names,paste("b1.",dnames,sep=""))
#    model
#}
#

#########################################################
InitErgm.b1kappa<-function(nw, m, arglist, ...) {
  ergm.checkdirected("b1kappa", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("b1kappa", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="b1kappa", soname="ergm",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"b1kappa")
  m
}

##########################################################
InitErgm.b1share<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b1share", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b1share", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b1share", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  d <- a$d
  nb1 <- get.network.attribute(nw,"bipartite")
  nb2 <- network.size(nw)-nb1
#
# Check for degeneracy
# 
  if(drop){
   mb1share <- paste("c(",paste(d,collapse=","),")",sep="")
   mb1share <- summary(
     as.formula(paste('nw ~ b1share(',mb1share,')',sep="")),
     drop=FALSE)
   if(any(mb1share==0)){
     cat(paste("Warning: There are no type 1 nodes sharing", d[mb1share==0], "type 2 nodes;\n",
       " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
    d <- d[mb1share!=0] 
   }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
#No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="b1share", soname="ergm",
                                inputs=c(0, ld, ld+1, nb1, d))
  m$coef.names<-c(m$coef.names,paste("b1share",d,sep=""))
  m
}

##########################################################
InitErgm.b2share<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b2share", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b2share", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b2share", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  d <- a$d
  nb1 <- get.network.attribute(nw,"bipartite")
  nb2 <- network.size(nw)-nb1
#
# Check for degeneracy
# 
  if(drop){
   mb2share <- paste("c(",paste(d,collapse=","),")",sep="")
   mb2share <- summary(
     as.formula(paste('nw ~ b2share(',mb2share,')',sep="")),
     drop=FALSE)
   if(any(mb2share==0)){
     cat(paste("Warning: There are no type 1 nodes sharing", d[mb2share==0], "type 2 nodes;\n",
       " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
    d <- d[mb2share!=0] 
   }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
#No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="b2share", soname="ergm",
                                inputs=c(0, ld, ld+1, nb1, d))
  m$coef.names<-c(m$coef.names,paste("b2share",d,sep=""))
  m
}

###################################### InitErgm TERMS:  B
#########################################################
InitErgm.bichange<-function (g, model, form=NULL, x=NULL, drop=TRUE, ...) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 3:5))
    stop(paste("bichange() model term expected 0 or 1 arguments, got ", 
                                nargs() - 3, sep = ""), call. = FALSE)
  if (nargs()==4){drop <- x; x <- NULL}
  #Coerce x to an adjacency matrix
  if(is.null(x)){
   xm<-as.edgelist(g)
   x<-"self"
  }else{
   if(is.network(x)){
    xm<-as.edgelist(x)
    x<-paste(quote(x))
   }else if(is.character(x)){
#   xm<-as.edgelist(g,x)
    xm<-get.network.attribute(g,x)
    xm<-as.edgelist(xm)
   }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
   }
  }
  nb1 <- get.network.attribute(g,"bipartite")
  nb2 <- network.size(g)-nb1
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
    formm <- matrix(1, nrow=nb1, ncol=nb2)
   }else{
    formm<-as.matrix(form)
    form<-paste(quote(form))
   }
  }
  #Check for matrix
  if (is.null(formm) || dim(formm)!=c(nb1, nb2)){
    stop("bichange requires a matrix of formation rates")
  }
  #Update the term number
  termnumber <- 1 + length(model$terms)
  #Update the terms list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_bichange ignores the value of inp->attrib).
  model$terms[[termnumber]] <- list(name = "bichange",
                                        soname="ergm",
             inputs = c(1, 1, 
                        2 + 2*nrow(xm)+nrow(formm)*ncol(formm),
                        nrow(xm), nrow(formm), as.double(c(xm, formm))
                       ) )
  #Update the coefficient name list, adding dyadcov.nameofx
  model$coef.names<-c(model$coef.names, paste(c("change"),x,sep="."))
  #Return the updated model list
  model
}

##########################################################
#InitErgm.bimixall<-function (nw, m, arglist, drop=TRUE, ...) {
#  a <- ergm.checkargs("bimix", arglist,
#    varnames = c("attrname","contrast"),
#    vartypes = c("character","logical"),
#    defaultvalues = list(NULL,FALSE),
#    required = c(TRUE, FALSE))
#  attach(a)
#  attrname<-a$attrname
#  nb1 <- get.network.attribute(nw,"bipartite")
#  nb2 <- network.size(nw)-nb1
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
#   mb2sb1 <- summary(
#     as.formula(paste('nw ~ bimixall("',attrname,'")',sep="")),
#     drop=FALSE)
#   if(any(mb2sb1==0)){
#    cat(paste("Warning: There are no counts for the bimixall", 
#      mb2sb1[mb2sb1==0],"\n"))
#    dropterms <- mb2sb1[mb2sb1==0]
#    cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#   }
#   u <- u[mb2sb1>0,] 
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
##########################################################
#InitErgm.bimixconddeg<-function (g, model, attrname, drop=TRUE, ...) {
#    if (!is.bipartite(g))
#      stop("The bimixconddeg term is for bipartite graphs.",
#           call.=FALSE)
#    if (!(nargs() %in% 4:5))
#        stop(paste("bimixconddeg() expected 1 arguments, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#
#    nb1 <- get.network.attribute(g,"bipartite")
#    nb2 <- network.size(g)-nb1
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
#     mb2sb1 <- summary(
#       as.formula(paste('g ~ bimixconddeg("',attrname,'")',sep="")),
#       drop=FALSE)
#     if(any(mb2sb1==0)){
#      cat(paste("Warning: There are no counts for the bimixconddeg", 
#          mb2sb1[mb2sb1==0],"\n"))
#      dropterms <- mb2sb1[mb2sb1==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     u <- u[mb2sb1>0,] 
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
##########################################################
#InitErgm.bimix<-function (g, model, attrname, drop=TRUE, ...) {
#    if (!is.bipartite(g))
#      stop("The bimix term is for bipartite graphs.",
#           call.=FALSE)
#    if (!(nargs() %in% 4:5))
#        stop(paste("bimix() expected 1 arguments, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#
#    nb1 <- get.network.attribute(g,"bipartite")
#    nb2 <- network.size(g)-nb1
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
#     mb2sb1 <- summary(
#       as.formula(paste('g ~ bimix("',attrname,'")',sep="")),
#       drop=FALSE)
#     if(any(mb2sb1==0)){
#      cat(paste("Warning: There are no counts for the bimix", 
#          mb2sb1[mb2sb1==0],"\n"))
#      dropterms <- mb2sb1[mb2sb1==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     u <- u[mb2sb1>0,] 
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
##########################################################
#InitErgm.biduration<-function (nw, m, arglist, ...) {
#  ergm.checkdirected("biduration", is.directed(nw), requirement=TRUE)
#  a <- ergm.checkargs("biduration", arglist,
#    varnames = c("form", "dissolve", "x"),
#    vartypes = c("matrix,network", "matrix,network", "matrix,network"),
#    defaultvalues = list(NULL, NULL, NULL),
#    required = c(TRUE, TRUE, FALSE))
#  attach(a)
#  x<-a$x;form<-a$form;dissolve<-a$dissolve
## m$coef.names<-c(m$coef.names, paste("biduration.",x,sep=""))
#  #Coerce x to an adjacency matrix
#  if(is.null(x)){
#    xm<-as.matrix.network(nw)
#    x<-"self"
#  }else{
#    if(is.network(x)){
#      xm<-as.matrix.network(x)
#      x<-paste(quote(x))
#    }else if(is.character(x)){
#      xm<-get.network.attribute(nw,x)
#      xm<-as.matrix.network(xm)
#    }else{
#      xm<-as.matrix(x)
#      x<-paste(quote(x))
#    }
#  }
#  #Check for symmetry
#  if (is.null(xm) || ncol(xm)!=2){
#    stop("biduration requires the edgelist of the base network")
#  }
#  nb1 <- get.network.attribute(nw,"bipartite")
#  nb2 <- network.size(nw)-nb1
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
#  if (is.null(formm) || dim(formm)!=c(nb1, nb2)){
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
#  if (is.null(dissolvem) || dim(dissolvem)!=c(nb1, nb2)){
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

#########################################################
InitErgm.biendure<-function (g, model, drop=TRUE, ...) 
{
  #Check to ensure that we got the right number of arguments
  if (!(nargs() %in% 3))
    stop(paste("biendure() model term expected 0 arguments, got ", 
                                nargs() - 3, sep = ""), call. = FALSE)
  xm<-as.edgelist(g)
  nb1 <- get.network.attribute(g,"bipartite")
  nb2 <- network.size(g)-nb1
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
                        nrow(xm), nb1, as.double(c(xm))
                       ) )
  #Update the coefficient name list, adding dyadcov.nameofx
  model$coef.names<-c(model$coef.names, paste("endure"))
  #Return the updated model list
  model
}


#########################################################
InitErgm.biduration<-function (nw, m, arglist, ...) {
  ergm.checkdirected("biduration", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("biduration", arglist,
    varnames = c("dissolve", "x"),
    vartypes = c("matrix,network", "matrix,network"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  x<-a$x;dissolve<-a$dissolve
# m$coef.names<-c(m$coef.names, paste("biduration.",x,sep=""))
  #Coerce x to an adjacency matrix
  if(is.null(x)){
    xm<-as.edgelist(nw)
    x<-"self"
  }else{
    if(is.network(x)){
      xm<-as.edgelist(x)
      x<-paste(quote(x))
    }else if(is.character(x)){
      xm<-get.network.attribute(nw,x)
      xm<-as.edgelist(xm)
    }else{
      xm<-as.matrix(x)
      x<-paste(quote(x))
    }
  }
  #Check for symmetry
  if (is.null(xm) || ncol(xm)!=2){
    stop("biduration requires the edgelist of the base network")
  }
  nb1 <- get.network.attribute(nw,"bipartite")
  nb2 <- network.size(nw)-nb1
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
  if (is.null(dissolvem) || dim(dissolvem)!=c(nb1, nb2)){
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

###################################### InitErgm TERMS:  E
#########################################################
InitErgm.b2kappa<-function(nw, m, arglist, ...) {
  ergm.checkdirected("b2kappa", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("b2kappa", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="b2kappa", soname="ergm",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"b2kappa")
  m
}

##########################################################
#InitErgm.b2sb1<-function(g, model, d, drop=TRUE, ...)
#{
#    if (nargs()!=4)
#        stop(paste("b2sb1() model term expected 1 argument, got ", 
#            nargs()-3, sep=""), call.=FALSE)
#    nb2 <- is.bipartite(g)
#    if (!nb2)
#      stop("The b2sb1 term is for bipartite graphs.",
#           call.=FALSE)
#    nb1 <- get.network.attribute(g,"bipartite")
#    if (is.directed(g))
#      stop("the b2sb1() term is not allowed with a directed graph",
#           call.=FALSE)
##
##   Check for degeneracy
## 
#    if(drop){
#     mb2sb1 <- paste("c(",paste(d,collapse=","),")",sep="")
#     mb2sb1 <- summary(
#       as.formula(paste('g ~ b2sb1(',mb2sb1,')',sep="")),
#       drop=FALSE)
#     if(any(mb2sb1==0)){
#      cat(paste("Warning: There are no b2sb1", d[mb2sb1==0],"dyads.\n"))
#      dropterms <- paste("b2sb1", d[mb2sb1==0],sep="")
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      d <- d[mb2sb1!=0] 
#     }
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="b2sb1", soname="ergm",
#                                          inputs=c(0, ld, ld+1, nb1, d))
#    model$coef.names<-c(model$coef.names,paste("b2sb1",d,sep=""))
#    model
#}
#

##########################################################
#InitErgm.b2<-function(g, model, drop=TRUE, ...)
#{
#    nb2 <- is.bipartite(g)
#    if (!nb2)
#      stop("The b2 term is for bipartite graphs.",
#           call.=FALSE)
#    if (nargs()!=3)
#        stop(paste("b2 model term expected zero argument, got ", 
#            nargs()-2, sep=""), call.=FALSE)
##
##   Check for degeneracy
## 
#    nb1 <- get.network.attribute(g,"bipartite")
#    nb2 <- network.size(g)-nb1
#    xnames <- network.vertex.names(g)
#    d <- nb1 + (2:nb2)
#    if(is.null(xnames)){
#     dnames <- paste("b2",d-nb1,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    if(drop){
#     degrees <- summary(g ~ b2, drop=FALSE)
#     if(any(degrees==0)){
#      cat(paste("Warning: There are no b1s for the b2", 
#          dnames[degrees==0],"\n"))
#      dropterms <- dnames[degrees==0]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     if(any(degrees==nb1)){
#      cat(paste("Warning: The b2", 
#          dnames[degrees],"is common to all b1s\n"))
#      dropterms <- dnames[degrees==nb1]
#      cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
#     }
#     d <- d[degrees!=nb1 & degrees!=0] 
#    }
#    ld<-length(d)
#    if(ld==0){return(model)}
##
#    if(is.null(xnames)){
#     dnames <- paste("b2",d-nb1,sep="")
#    }else{
#     dnames <- xnames[d]
#    }
#    termnumber<-1+length(model$terms)
##  No covariates here, so input component 1 is arbitrary
#    model$terms[[termnumber]] <- list(name="b2", soname="ergm",
#                                          inputs=c(0, ld, ld+1, nb1, d))
##   model$coef.names<-c(model$coef.names,dnames)
#    model$coef.names<-c(model$coef.names,paste("b2.",dnames,sep=""))
#    model
#}
#

###################################### InitErgm TERMS:  G
#########################################################
InitErgm.gwb2706<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("biduration", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb2", arglist,
    varnames = c("alpha"),
    vartypes = c("numeric"),
    defaultvalues = list(0.5),
    required = c(FALSE))
  attach(a)
  alpha<-a$alpha
  nb1 <- get.network.attribute(nw,"bipartite")
  nb2 <- network.size(nw)-nb1
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="gwb2",
                                    soname="ergm",
                                    inputs=c(0, 1, 2, nb1, alpha))
  m$coef.names<-c(m$coef.names,"gwb2")
  m
}

#########################################################
InitErgm.gwb1706<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("biduration", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb1", arglist,
    varnames = c("alpha"),
    vartypes = c("numeric"),
    defaultvalues = list(0.5),
    required = c(FALSE))
  attach(a)
  alpha<-a$alpha
  nb1 <- get.network.attribute(nw,"bipartite")
  nb2 <- network.size(nw)-nb1
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="gwb1",
                                    soname="ergm",
                                    inputs=c(0, 1, 2, nb1, alpha))
  m$coef.names<-c(m$coef.names,"gwb1")
  m
}

###################################### InitErgm TERMS:  M
#########################################################
InitErgm.monopolymixmat<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("monopolymixmat", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("monopolymixmat", arglist,
                   varnames = NULL,
                   vartypes = NULL,
                   defaultvalues = list(),
                   required = NULL)
  termnumber<-1+length(m$terms)
  m$coef.names<-c(m$coef.names, c("monoFmonoM", "monoFpolyM", "polyFmonoM"))
  m$terms[[termnumber]] <- list(name = "monopolymixmat", soname="ergm",
                                inputs = c(0, 3, 0))
  m
}


#########################################################
InitErgm.gwb1share<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwb1share", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("gwb1share", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb1share", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  termnumber<-1+length(m$terms)
  dname <- "b1share"
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map<- function(x,n,...) {
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwb1share=NULL,gwb1share.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwb1share#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwb1share") # must match params$gwb1share above
  }else{ # fixed == TRUE
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,paste("gwb1share.fixed.",alpha,sep=""))
  }
  m
}

#########################################################
InitErgm.gwb2share<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwb2share", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("gwb2share", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb2share", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  termnumber<-1+length(m$terms)
  dname <- "b2share"
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map<- function(x,n,...) {
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwb2share=NULL,gwb2share.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwb2share#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwb2share") # must match params$gwb2share above
  }else{ # fixed == TRUE
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,paste("gwb2share.fixed.",alpha,sep=""))
  }
  m
}
