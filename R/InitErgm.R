# Upon encountering a model term such as [name](args), ergm and related
# routines will call a function of the form InitErgm.[name].  The specific
# function call will be of the form
#   InitErgm.[name](network, model, arguments, ...)
# where network is the network object, model is a model object that should be
# updated and then returned, arguments is the list (if any) of arguments
# passed to the model term by the user, and ... includes any arguments
# passed to the InitErgm function from within the program.
#
# Arguments of
# the latter type include such items as drop (a logical flag telling whether
# degenerate terms should be dropped) and expanded (a logical flag used
# by curved exponential family terms).  Because such arguments are usually
# passed to ALL InitErgm functions, regardless of whether they are used,
# it is important that each InitErgm function declaration include the
# dot-dot-dot (...) argument.  Finally, such arguments are not guaranteed
# to be passed when the InitErgm function is called, so any InitErgm function
# requiring such an argument should supply a default value.
#
# An example:  If drop=TRUE is passed from inside the program,
# then the statement
#     ergm(nw ~ triangle + kstar (2:4) + nodematch("sex"))
# results in the following function calls:
#     InitErgm.triangle (nw, model, list(), drop=TRUE)
#     InitErgm.kstar (nw, model, list(2:4), drop=TRUE)
#     InitErgm.nodematch (nw, model, list("sex"), drop=TRUE)
#
# Each InitErgm.[name] function should check its argument list for errors, 
# then set termnumber to 1+length(model$terms).
# Next, it should add the names of the statistics that
# will be computed to the vector model$coef.names.  These names must be
# concatenated onto model$coef.names in the same order they will be produced
# by the changestat function.
# Finally, it should create 
# model$terms[[termnumber]] , a list with the following elements, some
# required and some optional:
#
# Required arguments of model$terms[[termnumber]]
# -----------------------------------------------
#    name: This is the (text) name of the term.  It is expected that there
#          is a C function called d_[name].
#  soname: This is the (text) name of the package containing the C function
#          called d_[name].
#  inputs: This is a (numeric) vector with at least 3 elements, as described
#          below:
#    Element 1 -- For functions that require a vector of covariates, either
#                 nodal or dyadic, this optional value is the number of
#                 input parameters BEFORE the beginning of the covariate
#                 vector.  For instance, if there are no input parameters
#                 passed before the covariate vector, this value should be
#                 set to zero.  The changestat function in C will be passed a
#                 pointer to the start of this vector of covariates, though
#                 the changestat function may choose to ignore this pointer,
#                 in which case the value of element 1 is arbitrary.
#    Element 2 -- The number of change statistics returned by the function.
#    Element 3 -- The total number of input parameters and covariates
#                 to be passed to the function.  If there are no nodal or 
#                 dyadic covariates, the value of element 1 is arbitrary.
#   Element 4+ -- The input parameters to be passed to the function.
#                 For example, if element 3 equals 3, then elements
#                 4, 5, 6 are the parameters to be passed.  No 4th element
#                 is necessary if element 3==0.  If there are nodal or
#                 dyadic covariates, they should be appended after any other
#                 input parameters (and element 1 may then be set to the
#                 number of other input parameters excluding the covariates).
#
# Optional arguments of model$terms[[termnumber]]
# -----------------------------------------------
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  Default value:  TRUE
#        params: For curved exponential family models, this argument must be
#                a list:  Each item in the list should be named with the
#                corresponding parameter name (one or more of these will
#                probably coincide with the coef.names used when
#                initialfit=TRUE; the initial values of such parameter values
#                will be set by MPLE, so their values in params are ignored.)
#                Any parameter not having its initial value set by MPLE
#                should be given its initial value in this params list.
#           eta: A function that gives the map from theta (the canonical
#                parameters associated with the statistics for this term)
#                to eta (the corresponding curved parameters).  The length
#                of eta is the same as the length of the params list above.
#                This function takes two args:  theta and length(eta).
#      gradient: A function that gives the gradient of the eta map above.
#                If theta has length p and eta has length q, then gradient
#                should return a p by q matrix.
#                This function takes two args:  theta and length(eta).
#  emptynwstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  Example:  If the degree0 term is
#                among the statistics, this argument is necessary because
#                degree0 = number of nodes for the empty network.


###################################### InitErgm TERMS:  A
#########################################################
InitErgm.absdiff<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("absdiff", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  termnumber<-1+length(m$terms)
  m$coef.names<-c(m$coef.names,paste("absdiff",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "absdiff", numeric=TRUE)
  m$terms[[termnumber]] <- list(name="absdiff", soname="ergm",
                                inputs=c(0,1,length(nodecov),nodecov),
                                dependence=FALSE)
  m
}

#########################################################
InitErgm.absdiffcat<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("absdiffcat", arglist,
                   varnames = c("attrname","base"),
                   vartypes = c("character","numeric"),
                   defaultvalues = list(NULL,NULL),
                   required = c(TRUE,FALSE))
  # base:  If not NULL, indices of non-zero absdiff categories to delete
  # (ordered from 1=smallest to largest).
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nodecov <- get.node.attr(nw, attrname, "absdiffcat")  
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=TRUE)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(base)) u <- u[-base]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  termnumber<-1+length(m$terms)  
  u2 <- u[!is.na(u)]
  m$terms[[termnumber]] <- list(name="absdiffcat", soname="ergm",
                                inputs=c(length(u2)+1, length(u),
                                         length(u2)+1+length(nodecov),
                                         u2, NAsubstitute, nodecov),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,paste("absdiff", attrname, u, sep="."))
  m
}

#########################################################
InitErgm.altkstar<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("altkstar", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("altkstar", arglist,
    varnames = c("lambda","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(1, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  lambda<-a$lambda;fixed<-a$fixed
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(x[2]*((1-1/x[2])^i + i) - 1)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(x[2]*((1-1/x[2])^i + i) - 1,
            x[1]*(i - 1 + (x[2]*x[2]-x[2]+i)*((1-1/x[2])^(i-1))/(x[2]*x[2]) )
           )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altkstar=NULL,
                                    altkstar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altkstar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altkstar", soname="ergm",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altkstar")
  }
  m
}

#########################################################
InitErgm.asymmetric<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("asymmetric", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("asymmetric", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    nasymmetric <- summary(as.formula('nw ~ asymmetric'), drop=FALSE)
    if(nasymmetric==0){
      cat(" ")
      cat(paste("Warning: There are no asymmetric ties;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'asymmetric' term has been dropped.\n"))
      return(m)
    }
    if(nasymmetric==network.dyadcount(nw)){
      cat(" ")
      cat(paste("Warning: All dyads have asymmetric ties!\n",
                 " the corresponding coefficient has been fixed at its MLE of infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'asymmetric' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="asymmetric", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"asymmetric")
  m
}

###################################### InitErgm TERMS:  B
#########################################################
InitErgm.b1concurrent<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b1concurrent", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b1concurrent", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b1concurrent", arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  attach(a)
  attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")       
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b1concurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b1concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    if(drop){ #   Check for degeneracy
      b1concurrentattr <- paste('nw ~ b1concurrent("',attrname,'")',sep="")
      b1concurrentattr <- summary(as.formula(b1concurrentattr),
                                 drop=FALSE) == 0
      if(any(b1concurrentattr)){
        cat(" ")
        cat(paste("Warning: There are no b1concurrent", ".", attrname,
                           u[b1concurrentattr],
                  "b1s;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("b1concurrent", ".", attrname,
                           u[b1concurrentattr], sep="")
        u <- u[-b1concurrentattr]
      }
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb1concurrent <- summary(
                          as.formula(paste('nw ~ b1concurrent',sep="")),
                          drop=FALSE) == 0
      if(any(mb1concurrent)){
        cat(paste("Warning: There are no concurrent b1s.\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b1concurrent_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           u, nodecov),
                                  dependence=TRUE)
    # See comment in d_b1concurrent_by_attr function
    m$coef.names<-c(m$coef.names, paste("b1concurrent",".", attrname,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b1concurrent", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b1concurrent",sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.b1degree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b1degree", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b1degree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b1degree", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b1degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b1degree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      b1degreeattr <- summary(as.formula
                             (paste('nw ~ b1degree(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(b1degreeattr)){
        cat(" ")
        cat(paste("Warning: There are no degree", du[1,b1degreeattr], ".",
                   attrname, u[du[2,b1degreeattr]],
                  "b1s;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        du <- matrix(du[,!b1degreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nb1]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb1degree <- paste("c(",paste(d,collapse=","),")",sep="")
      mb1degree <- summary(
                          as.formula(paste('nw ~ b1degree(',mb1degree,')',sep="")),
                          drop=FALSE) == 0
      if(any(mb1degree)){
        cat(" ")
        cat(paste("Warning: There are no degree", d[mb1degree],
                  "b1s;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        d <- d[!mb1degree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nb1
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b1degree_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_b1degree_by_attr function
    m$coef.names<-c(m$coef.names, paste("adeg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b1degree", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b1degree",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.b1factor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b1factor", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b1factor", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b1factor", arglist,
    varnames = c("attrname", "base"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, 1),
    required = c(TRUE, FALSE))                                    
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, attrname, "b1factor")[1:nb1]
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop){
    nfc <- summary(as.formula(paste('nw ~ b1factor("',attrname,'",base=0)',sep="")),drop=FALSE) == 0
    if(any(nfc)){
      dropterms <- paste(paste("b1factor",attrname,sep="."),u[nfc],sep=".")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      u<-u[!nfc]
      ui<-ui[!nfc]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to b1factor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)  
  if(base[1]==0){
   m$terms[[termnumber]] <- list(name="b1factor", soname="ergm",
                                 inputs=c(lu, 
                                          lu, 
                                          lu+length(nodecov),
                                          ui, nodecov), dependence=FALSE)
   m$coef.names<-c(m$coef.names, paste("b1factor",
                                       attrname, paste(u), sep="."))
  }else{
   m$terms[[termnumber]] <- list(name="b1factor", soname="ergm",
                                 inputs=c(lu-length(base), 
                                          lu-length(base), 
                                          lu-length(base)+length(nodecov),
                                          ui[-base], nodecov), dependence=FALSE)
   m$coef.names<-c(m$coef.names, paste("b1factor",
                                       attrname, paste(u[-base]), sep="."))
  }
  m
}

#########################################################
InitErgm.b1star<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b1star", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b1star", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b1star", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)                           
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b1star")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b1star() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      b1starattr <- summary(as.formula
                             (paste('nw ~ b1star(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(b1starattr)){
        cat(" ")
        cat(paste("Warning: There are no b1 stars", du[1,b1starattr], ".",
                   attrname, u[du[2,b1starattr]],
                  ";\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        du <- matrix(du[,!b1starattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nb1]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb1star <- paste("c(",paste(d,collapse=","),")",sep="")
      mb1star <- summary(
                          as.formula(paste('nw ~ b1star(',mb1star,')',sep="")),
                          drop=FALSE) == 0
      if(any(mb1star)){
        cat(" ")
        cat(paste("Warning: There are no b1 stars", d[mb1star],
                  "b1s;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        d <- d[!mb1star] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nb1
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}            
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ostar", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("b1star", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ostar", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b1star",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.b2concurrent<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b2concurrent", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b2concurrent", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b2concurrent", arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  attach(a)
  attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b2concurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b2concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    if(drop){ #   Check for degeneracy
      b2concurrentattr <- paste('nw ~ b2concurrent("',attrname,'")',sep="")
      b2concurrentattr <- summary(as.formula(b2concurrentattr),
                                 drop=FALSE) == 0
      if(any(b2concurrentattr)){
        cat(" ")
        cat(paste("Warning: There are no b2concurrent", ".", attrname,
                           u[b2concurrentattr],
                  "b2s;\n",
                 " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        u <- u[-b2concurrentattr]
      }
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb2concurrent <- summary(
                          as.formula(paste('nw ~ b2concurrent',sep="")),
                          drop=FALSE) == 0
      if(any(mb2concurrent)){
        cat(paste("Warning: There are no concurrent b2s\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b2concurrent_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           u, nodecov),
                                  dependence=TRUE)
    # See comment in d_b2concurrent_by_attr function
    m$coef.names<-c(m$coef.names, paste("b2concurrent",".", attrname,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b2concurrent", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b2concurrent",sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.b2degree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b2degree", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b2degree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b2degree", arglist,
    varnames = c("d", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b2degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to b2degree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      b2degreeattr <- summary(
       as.formula(paste('nw ~ b2degree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(b2degreeattr)){
        dropterms <- paste("b2deg", du[1,b2degreeattr], ".", attrname,
                           u[du[2,b2degreeattr]], sep="")
        cat("Warning: These b2degree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        du <- matrix(du[,!b2degreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[(1+nb1):n]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb2degree <- paste("c(",paste(d,collapse=","),")",sep="")
      mb2degree <- summary(
       as.formula(paste('nw ~ b2degree(',mb2degree,')',sep="")),
       drop=FALSE) == 0
      if(any(mb2degree)){
        cat(" ")
        cat(paste("Warning: There are no degree", d[mb2degree],
                  "b2s;\n",
                  " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
       dropterms <- paste("b2degree", d[mb2degree],sep="")
       d <- d[!mb2degree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))             
      emptynwstats[d==0] <- n-nb1
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b2degree_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_b2degree_by_attr function
    m$coef.names<-c(m$coef.names, paste("b2deg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="b2degree", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b2degree",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.b2factor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b2factor", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b2factor", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b2factor", arglist,
    varnames = c("attrname", "base"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, 1),
    required = c(TRUE, FALSE))
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nb1 <- get.network.attribute(nw, "bipartite")
  nodecov <- get.node.attr(nw, attrname, "b2factor")[(nb1+1):network.size(nw)]
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop){
    nfc <- summary(as.formula(paste('nw ~ b2factor("',attrname,'",base=0)',sep="")),drop=FALSE) == 0
    if(any(nfc)){
      dropterms <- paste(paste("b2factor",attrname,sep="."),u[nfc],sep=".")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      u<-u[!nfc]
      ui<-ui[!nfc]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to b2factor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)  
  if(base[1]==0){
   m$terms[[termnumber]] <- list(name="b2factor", soname="ergm",
                                 inputs=c(lu, 
                                          lu, 
                                          lu+length(nodecov),
                                          ui, nodecov), dependence=FALSE)
   m$coef.names<-c(m$coef.names, paste("b2factor",
                                       attrname, paste(u), sep="."))
  }else{
   m$terms[[termnumber]] <- list(name="b2factor", soname="ergm",
                                 inputs=c(lu-length(base), 
                                          lu-length(base), 
                                          lu-length(base)+length(nodecov),
                                          ui[-base], nodecov), dependence=FALSE)
   m$coef.names<-c(m$coef.names, paste("b2factor",
                                       attrname, paste(u[-base]), sep="."))
  }
  m
}

#########################################################
InitErgm.b2star<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("b2star", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("b2star", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("b2star", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nb1 <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "b2star")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to b2star() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      b2starattr <- summary(as.formula
                             (paste('nw ~ b2star(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(b2starattr)){
        dropterms <- paste("b2star", du[1,b2starattr], ".", attrname,
                           u[du[2,b2starattr]], sep="")
        cat("Warning: These b2star terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        du <- matrix(du[,!b2starattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nb1]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mb2star <- paste("c(",paste(d,collapse=","),")",sep="")
      mb2star <- summary(
                          as.formula(paste('nw ~ b2star(',mb2star,')',sep="")),
                          drop=FALSE) == 0
      if(any(mb2star)){
        cat(paste("Warning: There are no order", d[mb2star],"b2stars.\n"))
        dropterms <- paste("b2star", d[mb2star],sep="")
        cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
        d <- d[!mb2star] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nb1
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="istar", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("b2star", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="istar", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("b2star",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.balance<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("balance", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
   nodecov <- get.node.attr(nw, attrname, "balance")
   u<-sort(unique(nodecov))
   if(any(is.na(nodecov))){u<-c(u,NA)}
#
#  Recode to numeric if necessary
#
   nodecov <- match(nodecov,u,nomatch=length(u)+1)
   ui <- seq(along=u)

   if (length(u)==1)
         stop ("Attribute given to balance() has only one value", call.=FALSE)
#
#  Check for degeneracy
#
   if(drop){
      triattr <- summary(
       as.formula(paste('nw ~ balance(','"',attrname,'",diff=',diff,')',sep="")),
       drop=FALSE) == 0
      if(diff){
       if(any(triattr)){
        dropterms <- paste(paste("balance",attrname,sep="."),u[triattr],sep="")
      cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        u <- u[!triattr] 
        ui <- ui[!triattr] 
       }
      }else{
       if(triattr){
         dropterms <- paste(paste("balance",attrname,sep="."),sep="")
      cat(" ")
         cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#        cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                  "have been dropped.\n"))
       }
      }
     }
     if (!diff) {
#     No parameters before covariates here, so input element 1 equals 0
      m$terms[[termnumber]] <- list(name="balance", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names,paste("balance",attrname,sep="."))
     } else {
      #  Number of input parameters before covariates equals number of
      #  unique elements in nodecov, namely length(u), so that's what
      #  input element 1 equals
      m$terms[[termnumber]] <- list(name="balance", soname="ergm",
          inputs=c(length(ui), length(ui), length(ui)+length(nodecov),
                   ui, nodecov),
          dependence=TRUE)
      m$coef.names<-c(m$coef.names,paste("balance",
          attrname, u, sep="."))
     }
  }else{
#  No attributes (or diff)
#
#   Check for degeneracy
#   Can't do this as starts an infinite loop
#
#   if(drop){
#    triattr <- summary(as.formula('nw ~ balance'), drop=FALSE) == 0
#    if(triattr){
#       cat(paste("Warning: There are no balanced triads;\n",
#       cat(paste("To avoid degeneracy the balance term has been dropped.\n"))
#    }
#   }
#   No covariates, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="balance", soname="ergm",
                                  inputs=c(0,1,0),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,"balance")
   }
   m
}


###################################### InitErgm TERMS:  C
#########################################################
InitErgm.concurrent<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("concurrent", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("concurrent", arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  attach(a)
  attrname <- a$attrname
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "concurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    if(drop){ #   Check for degeneracy
      concurrentattr <- summary(as.formula
                             (paste('nw ~ concurrent(',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(concurrentattr)){
        dropterms <- paste("concurrent", ".", attrname,
                           u[concurrentattr], sep="")
      cat(" ")
        cat("Warning: These concurrent terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        u <- u[-concurrentattr]
      }
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mconcurrent <- summary(
                          as.formula(paste('nw ~ concurrent',sep="")),
                          drop=FALSE) == 0
      if(any(mconcurrent)){
      cat(" ")
        cat(paste("Warning: There are no concurrent b1s;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        return(m)
      }
    }                                                         
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="concurrent_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           u, nodecov),
                                  dependence=TRUE)
    # See comment in d_concurrent_by_attr function
    m$coef.names<-c(m$coef.names, paste("concurrent",".", attrname,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="concurrent", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("concurrent",sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.ctriple<-InitErgm.ctriad<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ctriple", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ctriple", arglist,
    varnames = c("attrname","diff"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,FALSE),
    required = c(FALSE,FALSE))
  attach(a)
  attrname <- a$attrname; diff <- a$diff;
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    nodecov <- get.node.attr(nw, attrname, "ctriple")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ctriple() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ ctriple(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("ctriple",attrname,sep="."),
                             u[triattr],sep="")
      cat(" ")
          cat(paste("Warning: The count of",
                paste(dropterms,collapse=" and, "),
                    "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=""))
#         cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("ctriple",attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="ctriple", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("ctriple",attrname,sep="."))
    } else {
      #  Number of input parameters before covariates equals number of
      #  unique elements in nodecov, namely length(u), so that's what
      #  input element 1 equals
      m$terms[[termnumber]] <- list(name="ctriple", soname="ergm",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("ctriple", attrname, u, sep="."))
    }
  }else{
#    No attributes (or diff)
#    No covariates, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ctriple", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"ctriple")
  }
  m
}

#########################################################
InitErgm.cycle<-function(nw, m, arglist, drop=TRUE, ...)
{
  #Verify arguments
  a <- ergm.checkargs("cycle", arglist,
    varnames = "k",
    vartypes = "numeric",
    defaultvalues = list(NULL),
    required = TRUE)
  attach(a)
  k<-a$k
  #Check for degeneracy
  if(drop){
    mcycle <- paste("c(",paste(k,collapse=","),")",sep="")
    mcycle <- summary(as.formula(paste('nw ~ cycle(',mcycle,')',sep="")),
      drop=FALSE) == 0
    if(any(mcycle)){
      cat(" ")
      cat(paste("Warning: There are no order", k[mcycle],"cycles;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("cycle", k[mcycle],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      k <- k[!mcycle] 
    }
  }
  #Set things up
  lk<-length(k)             #Find the number of terms remaining
  if(lk==0){return(m)}        #Return if no terms left
  mk<-max(k)                #Get maximum cycle length
  usestats<-(2:mk)%in%k     #Which stats are being used?
  direct<-is.directed(nw)   #Is the graph directed?
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="cycle", soname="ergm",
                                inputs=c(0, lk, mk+1, direct, mk, usestats))
  m$coef.names<-c(m$coef.names,paste("cycle",k,sep=""))
  m
}

###################################### InitErgm TERMS:  D
#########################################################
InitErgm.degree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("degree", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("degree", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "degree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degree() has only one value", call.=FALSE)
  }
  if(!is.null(attrname) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      degreeattr <- summary(
       as.formula(paste('nw ~ degree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(degreeattr)){
        dropterms <- paste("deg", du[1,degreeattr], ".", attrname,
                           u[du[2,degreeattr]], sep="")
        cat(" ")
        cat("Warning: These degree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        du <- matrix(du[,!degreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(drop){
      tmp <- paste("c(",paste(d,collapse=","),")",sep="")
      if(!homophily) {
        mdegree <- summary(as.formula(paste('nw ~ degree(',tmp,')',
                                            sep="")), drop=FALSE) == 0
      } else {
        mdegree <- summary(as.formula(paste('nw ~ degree(',tmp,',"',attrname,
                                                         '", TRUE)', sep="")), 
                                             drop = FALSE) == 0
      }
      if(any(mdegree)){
      cat(" ")
        cat("Warning: These degree terms have extreme counts and will be dropped:\n")
        cat(d[mdegree], "\n", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        d <- d[!mdegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  termnumber<-1+length(m$terms)
  if(is.null(attrname)) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("degree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="degree_w_homophily", soname="ergm",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    # See comment in d_degree_w_homophily function
    m$coef.names<-c(m$coef.names,paste("deg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="degree_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_degree_by_attr function
    m$coef.names<-c(m$coef.names, paste("deg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }
  if (!is.null(emptynwstats))
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.density<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("density", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="density", soname="ergm",
                                inputs=c(0, 1, 0),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,"density")
  m
}

#########################################################
InitErgm.dsp<-function(nw, m, arglist, drop=TRUE, ...) {
# ergm.checkdirected("dsp", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("dsp", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  if(drop){
    mdsp <- paste("c(",paste(d,collapse=","),")",sep="")
    mdsp <- summary(as.formula(paste('nw ~ dsp(',mdsp,')',sep="")),
                    drop=FALSE)
    if(any(mdsp==0)){
      cat(" ")
      cat(paste("Warning: There are no dsp", d[mdsp==0],"dyads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("dsp", d[mdsp==0],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- d[mdsp!=0] 
    }
  }
  if (any(d==0)) {
    emptynwstats <- rep(0, length(d))
    if(is.bipartite(nw)){
      nb1 <- get.network.attribute(nw, "bipartite")
      nb2 <- network.size(nw) - nb1
      emptynwstats[d==0] <- nb1*(nb1-1)/2 + nb2*(nb2-1)/2
    }else{
      emptynwstats[d==0] <- network.dyadcount(nw)
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
  m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("dsp",d,sep=""))
  m
}

#########################################################
InitErgm.dyadcov<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("dyadcov", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrixnetwork","character"),
    defaultvalues = list(NULL,NULL),
    required = c(TRUE,FALSE))
  attach(a)
  x<-a$x;attrname<-a$attrname
  #Coerce x to an adjacency matrix
  if(is.network(x))
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
  else if(is.character(x))
#   xm<-as.matrix.network(nw,matrix.type="adjacency",x)
    xm<-get.network.attribute(nw,x)
  else
    xm<-as.matrix(x)

  if(is.directed(nw)){
   #Check for symmetry
   # DH:  Since nw is directed, why are we testing for symmetry here?  
   if (any(xm[upper.tri(xm)]!=t(xm)[upper.tri(xm)])){
     xm[lower.tri(xm)]<-t(xm)[lower.tri(xm)]
     warning("asymmetric covariate in dyadcov; using upper triangle only")
   }
   #Update the term number
   termnumber <- 1 + length(m$terms)
   #Update the terms list, adding the vectorized adjacency matrix

#  There is 1 input parameter before the covariate vector, so input
#  element 1 is set to 1 (although in this case, input element 1
#  is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "dyadcov",  soname="ergm",
#                                inputs = c(1, 3, 1+NROW(xm)*NROW(xm),
                                 inputs = c(1, 3, 1+length(xm),
                                   NCOL(xm), as.double(xm)),
                                 dependence=FALSE)
   if(!is.null(attrname))
     cn<-paste("dyadcov", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
   else
     cn<-paste("dyadcov", as.character(sys.call(0)[[4]][2]), sep = ".")
   m$coef.names <- c(m$coef.names, paste(cn, c("mutual","utri","ltri"),
                                         sep=".") )
  }else{
#  So it is undirected
   termnumber <- 1 + length(m$terms)
#  There is 1 input parameter before the covariate vector, so input
#  element 1 is set to 1 (although in this case, input element 1
#  is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "dyadcov", soname="ergm", 
                                 inputs = c(1, 1, 1+length(xm),
                                   NCOL(xm), as.double(xm)),
                                 dependence=FALSE)
   if(!is.null(attrname))
     cn<-paste("dyadcov", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
   else
     cn<-paste("dyadcov", as.character(sys.call(0)[[4]][2]), sep = ".")
   m$coef.names <- c(m$coef.names, cn)
  }
  m
}

###################################### InitErgm TERMS:  E
#########################################################
InitErgm.edgecov<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("edgecov", arglist,
    varnames = c("x", "attrname"),
    vartypes = c("matrixnetwork", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  x<-a$x;attrname<-a$attrname
  #Coerce x to an adjacency matrix
  if(is.network(x))
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
  else if(is.character(x))
#   xm<-as.matrix.network(nw,matrix.type="adjacency",x)
    xm<-get.network.attribute(nw,x)
  else
    xm<-as.matrix(x)
  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1 (although in this case, input element 1
# is actually arbitrary since d_edgecov ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "edgecov", soname="ergm", 
                                inputs = c(1, 1, 1+length(xm),
                                  NCOL(xm), as.double(xm)),
                                dependence=FALSE)
#                               inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
#                                 NROW(xm), as.double(xm)),
  if(!is.null(attrname))
    cn<-paste("edgecov", as.character(sys.call(0)[[4]][2]), 
              as.character(attrname), sep = ".")
  else
    cn<-paste("edgecov", as.character(sys.call(0)[[4]][2]), sep = ".")
  m$coef.names <- c(m$coef.names, cn)
  m
}

#########################################################
InitErgm.edges<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("edges", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="edges", soname="ergm",
                                inputs=c(0, 1, 0),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,"edges")
  m
}

#########################################################
InitErgm.esp<-function(nw, m, arglist, drop=TRUE, ...) {
# ergm.checkdirected("esp", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("esp", arglist,
    varnames = c("d"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  d<-a$d
  if(drop){
    mesp <- paste("c(",paste(d,collapse=","),")",sep="")
    mesp <- summary(as.formula(paste('nw ~ esp(',mesp,')',sep="")),
                    drop=FALSE)
    if(any(mesp==0)){
      cat(" ")
      cat(paste("Warning: There are no dyads with esp", d[mesp==0],";\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("esp", d[mesp==0],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- d[mesp!=0] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
  m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("esp",d,sep=""))
  m
}

###################################### InitErgm TERMS:  F

###################################### InitErgm TERMS:  G

#########################################################
InitErgm.gwb1degree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwb1degree", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("gwb1degree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb1degree", arglist,
    varnames = c("decay", "fixed", "attrname"),
    vartypes = c("numeric", "logical", "character"),
    defaultvalues = list(0, FALSE, NULL),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  nb1 <- get.network.attribute(nw,"bipartite")
  d <- 1:(network.size(nw) - nb1)
  termnumber<-1+length(m$terms)
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwb1degree term is not yet able to handle a ",
           "non-fixed decay term.", call.=FALSE) # with an attribute.")
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
    m$terms[[termnumber]] <- list(name="b1degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwb1degree=NULL,
                                    gwb1degree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwb1degree#",d,sep=""))
  } else {
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb1degree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb1degree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(m)}
      #  No covariates here, so input component 1 is arbitrary
      m$terms[[termnumber]] <- list(name="gwb1degree_by_attr", soname="ergm",
                                    inputs=c(0, lu, 
                                             1+length(nodecov), 
                                             decay, nodecov),
                                    dependence=TRUE)
      # See comment in d_gwb1degree_by_attr function
      m$coef.names<-c(m$coef.names, paste("gwb1deg", decay, ".", 
                                          attrname, u, sep=""))
    }else{
      m$terms[[termnumber]] <- list(name="gwb1degree", soname="ergm",
                                    inputs=c(0, 1, 1, decay),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names,paste("gwb1deg",decay,sep=""))
    }
  }
  m
}
#########################################################
InitErgm.gwb2degree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwb2degree", is.directed(nw), requirement=FALSE)
  ergm.checkbipartite("gwb2degree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwb2degree", arglist,
    varnames = c("decay", "fixed", "attrname"),
    vartypes = c("numeric", "logical", "character"),
    defaultvalues = list(0, FALSE, NULL),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; fixed<-a$fixed; attrname<-a$attrname
  nb1 <- get.network.attribute(nw,"bipartite")
  d <- 1:nb1
  termnumber<-1+length(m$terms)
  if (!initialfit && !fixed) { # This is a curved exp fam
#    if (!is.null(attrname)) {
      stop("The gwb2degree term is not yet able to handle a ",
           "non-fixed decay term.", call.=FALSE) # with an attribute.")
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
    m$terms[[termnumber]] <- list(name="b2degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwb2degree=NULL,
                                    gwb2degree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwb2degree#",d,sep=""))
  } else { 
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwb2degree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwb2degree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(m)}
      #  No covariates here, so input component 1 is arbitrary
      m$terms[[termnumber]] <- list(name="gwb2degree_by_attr", soname="ergm",
                                    inputs=c(0, lu,
                                             1+length(nodecov), 
                                             decay, nodecov),
                                    dependence=TRUE)
      # See comment in d_gwb2degree_by_attr function
      m$coef.names<-c(m$coef.names, paste("gwb2deg", decay, ".", 
                                          attrname, u, sep=""))
    }else{
      m$terms[[termnumber]] <- list(name="gwb2degree", soname="ergm",
                                    inputs=c(0, 1, 1, decay),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names,paste("gwb2deg",decay,sep=""))
    }
  }
  m
}

#########################################################
InitErgm.gwdegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
 ergm.checkdirected("gwdegree", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwdegree", arglist,
    varnames = c("decay", "fixed", "attrname"),
    vartypes = c("numeric", "logical", "character"),
    defaultvalues = list(0, FALSE, NULL),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  termnumber<-1+length(m$terms)
  d <- 1:(network.size(nw)-1)
  if(!initialfit && !fixed){ # This is a curved exponential family model
    if (!is.null(attrname)) {
      stop("The gwdegree term is not yet able to handle a ",
           "nonfixed decay term with an attribute.", call.=FALSE)
    }
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegree=NULL,
                                    gwdegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegree#",d,sep=""))
  } else if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "gwdegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to gwdegree() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(nrow(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="gwdegree_by_attr", soname="ergm",
                                  inputs=c(0, lu, 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("gwdeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwdegree", soname="ergm",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwdegree")
  }
  m
}

#########################################################
InitErgm.gwdsp<-function(nw, m, arglist, initialfit=FALSE, ...) {
# ergm.checkdirected("gwdsp", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwdsp", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  termnumber<-1+length(m$terms)
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
    if(is.directed(nw)){dname <- "tdsp"}else{dname <- "dsp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdsp=NULL,gwdsp.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdsp#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    if(is.directed(nw)){dname <- "gwtdsp"}else{dname <- "gwdsp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwdsp") # must match params$gwdsp above
  }else{ # fixed == TRUE
    if(is.directed(nw)){dname <- "gwtdsp"}else{dname <- "gwdsp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,paste("gwdsp.fixed.",alpha,sep=""))
  }
  m
}

#########################################################
InitErgm.gwesp<-function(nw, m, arglist, initialfit=FALSE, ...) {
# ergm.checkdirected("gwesp", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwesp", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  termnumber<-1+length(m$terms)
  alpha=alpha[1] # Not sure why anyone would enter a vector here, but...
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...){
      i <- 1:n
      x[1]*exp(x[2])*(1-(1-exp(-x[2]))^i)
    }
    gradient <- function(x,n,...){
      i <- 1:n
      a <- 1-exp(-x[2])
      exp(x[2]) * rbind(1-a^i, x[1] * (1 - a^i - i*a^(i-1) ) )
    }
    if(is.directed(nw)){dname <- "tesp"}else{dname <- "esp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwesp=NULL,gwesp.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("esp#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    if(is.directed(nw)){dname <- "gwtesp"}else{dname <- "gwesp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, 1, alpha))
    m$coef.names<-c(m$coef.names,"gwesp") # Must match params$gwesp above
  }else{ # fixed == TRUE
    if(is.directed(nw)){dname <- "gwtesp"}else{dname <- "gwesp"}
    m$terms[[termnumber]] <- list(name=dname, soname="ergm",
                                  inputs=c(0, 1, 1, alpha))
    m$coef.names<-c(m$coef.names,paste("gwesp.fixed.",alpha,sep=""))
  }
  m
}

#########################################################
InitErgm.gwidegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwidegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("gwidegree", arglist,
                      varnames = c("decay", "fixed", "attrname"),
                      vartypes = c("numeric", "logical", "character"),
                      defaultvalues = list(0, FALSE, NULL),
                      required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  termnumber<-1+length(m$terms)
  d <- 1:(network.size(nw)-1)
  if(!initialfit && !fixed){ # This is a curved exponential family model
    if (!is.null(attrname)) {
      stop("The gwidegree term is not yet able to handle a ",
           "nonfixed decay term with an attribute.", call.=FALSE)
    }
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    m$terms[[termnumber]] <- list(name="idegree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwidegree=NULL,
                                    gwidegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwidegree#",d,sep=""))
  } else { 
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwidegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwidegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(m)}
      #  No covariates here, so input component 1 is arbitrary
      m$terms[[termnumber]] <- list(name="gwidegree_by_attr", soname="ergm",
                                    inputs=c(0, lu, 
                                             1+length(nodecov), 
                                             decay, nodecov),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names, paste("gwideg", decay, ".", 
                                          attrname, u, sep=""))
    }else{
      m$terms[[termnumber]] <- list(name="gwidegree", soname="ergm",
                                    inputs=c(0, 1, length(decay), decay))
      m$coef.names<-c(m$coef.names,paste("gwidegree.fixed.",decay,sep=""))
    }
  }
  m
}

#########################################################
InitErgm.gwodegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwodegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("gwodegree", arglist,
                      varnames = c("decay", "fixed", "attrname"),
                      vartypes = c("numeric", "logical", "character"),
                      defaultvalues = list(0, FALSE, NULL),
                      required = c(TRUE, FALSE, FALSE))
  attach(a)
  decay<-a$decay; attrname<-a$attrname; fixed<-a$fixed  
  termnumber<-1+length(m$terms)
  d <- 1:(network.size(nw)-1)
  if(!initialfit && !fixed){ # This is a curved exponential family model
    if (!is.null(attrname)) {
      stop("The gwodegree term is not yet able to handle a ",
           "nonfixed decay term with an attribute.", call.=FALSE)
    }
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    m$terms[[termnumber]] <- list(name="odegree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwodegree=NULL,
                                    gwodegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwodegree#",d,sep=""))
  } else {
    if(!is.null(attrname)) {
      nodecov <- get.node.attr(nw, attrname, "gwodegree")
      u<-sort(unique(nodecov))
      if(any(is.na(nodecov))){u<-c(u,NA)}
      nodecov <- match(nodecov,u) # Recode to numeric
      if (length(u)==1)
        stop ("Attribute given to gwodegree() has only one value", call.=FALSE)
      # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
      lu <- length(u)
      du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
      if(nrow(du)==0) {return(m)}
      #  No covariates here, so input component 1 is arbitrary
      m$terms[[termnumber]] <- list(name="gwodegree_by_attr", soname="ergm",
                                    inputs=c(0, lu, 
                                             1+length(nodecov), 
                                             decay, nodecov),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names, paste("gwodeg", decay, ".", 
                                          attrname, u, sep=""))
    }else{
      m$terms[[termnumber]] <- list(name="gwodegree", soname="ergm",
                                    inputs=c(0, 1, length(decay), decay))
      #   m$coef.names<-c(m$coef.names,paste("gwodegree.fixed.",decay,sep=""))
      m$coef.names<-c(m$coef.names,"gwodegree")
    }
  }
  m
}

###################################### InitErgm TERMS:  H
#########################################################
InitErgm.hamming<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("hamming", arglist=arglist,
    varnames = c("x","cov","attrname"),
    vartypes = c("matrixnetwork","matrixnetwork","character"),
    defaultvalues = list(nw,NULL,NULL),
    required = c(FALSE,FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  x<-a$x
  cov<-a$cov
  termnumber<-1+length(m$terms)

  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist",attrname)
    x2<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else if(is.null(x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")
  }else{
    xm<-as.matrix(x)
    x2<-paste(quote(x))
  }
  if (is.null(xm) || NCOL(xm)!=2){
    stop("hamming() requires an argument that can be coerced into an edgelist.",
         call.=FALSE)
  }
  if (is.null(cov)) { # calculate unweighted hamming distance
    if(drop){ #   Check for degeneracy
      hamm <- summary(nw ~ hamming(x), drop=FALSE) == 0
      if(hamm){
        cat(paste(" Warning: The Hamming distance is zero;\n",
                  " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        return(m)
      }
    }
    # There is 1 input parameter before the covariate vector, so input
    # element 1 is set to 1 although in this case, input element 1
    # is actually arbitrary since d_hamming ignores the value of inp->attrib.
    m$terms[[termnumber]] <- list(name = "hamming",  soname="ergm",
                                  inputs = c(0, 1, NROW(xm)*2+1, NROW(xm), 
                                             as.integer(xm)),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names, paste("hamming",
                                        as.character(sys.call(0)[[4]][2]),
                                        sep="."))
  } else {
    # Extract dyadic covariate
    if(is.network(cov)){
      covm<-as.matrix.network(cov,matrix.type="adjacency",attrname)
      cov<-paste(quote(cov))
    }else if(is.character(cov)){
      covm<-get.network.attribute(nw,cov)
      covm<-as.matrix.network(covm,matrix.type="adjacency")
    }else{
      covm<-as.matrix(cov)
      cov<-paste(quote(cov))
    }
    if (is.null(covm) || !is.matrix(covm) || nrow(covm)!=get.network.attribute(nw,"bipartite")){
      stop("Improper dyadic covariate passed to hamming()", call.=FALSE)
    }
    ##
    ##   Check for degeneracy
    ##
    #    if(drop){
      #     ematch <- mixmat[u]
      #     mu <- ematch==0
      #     mu[is.na(mu)] <- FALSE
      #     if(any(mu)){
        #      dropterms <- paste(paste("hamming.weighted",attrname,sep="."),
                                  #        apply(u,1,paste,collapse="")[mu],sep="")
        #      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        #      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
        #      u <- u[!mu,]
      #     }
    #    }
    # There is 1 input parameter before the covariate vector, so input
    # component 1 is set to 1 (although in this case, input component 1
                               # is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
    m$terms[[termnumber]] <- list(name = "hamming_weighted", soname="ergm",
                                  inputs = c(1, 1,
                                             1+2*nrow(xm)+nrow(covm)*ncol(covm),
                                             nrow(xm), as.integer(xm),
                                             as.double(covm)),
                                  dependence=FALSE)
    if(!is.null(attrname)){
      cn<-paste("hamming", as.character(sys.call(0)[[4]][2]), "wt",
                as.character(attrname), sep = ".")
    }else{
      cn<-paste("hamming", as.character(sys.call(0)[[4]][2]), "wt",
                as.character(sys.call(0)[[4]][3]), sep = ".")
    }
    m$coef.names <- c(m$coef.names, cn)
  }
  m
}

#########################################################
InitErgm.hammingmix.constant<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingconstantmix", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingmix.constant", arglist=arglist,
    varnames = c("attrname","x","base", contrast),
    vartypes = c("character","matrixnetwork","numeric","logical"),
    defaultvalues = list(NULL,nw,0,FALSE),
    required = c(TRUE,FALSE,FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  x<-a$x
  base<-a$base
  contrast<-a$contrast
  drop<-a$drop
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
    stop("hammingmix.constant() requires an edgelist", call.=FALSE)
  }
    nodecov <- get.node.attr(nw, attrname, "hammingmix.constant")
    mixmat <- mixingmatrix(nw,attrname)$mat
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
#
#   Recode to numeric if necessary
#
    namescov <- sort(unique(nodecov))
    nodecov <- match(nodecov,namescov)
    if (length(nodecov)==1)
        stop ("Argument to hammingmix.constant() has only one value", call.=FALSE)
##
##   Check for degeneracy
##
#    if(drop){
#     ematch <- mixmat[u]
#     mu <- ematch==0
#     mu[is.na(mu)] <- FALSE
#     if(any(mu)){
#      dropterms <- paste(paste("hammingconstantmix",attrname,sep="."),
#        apply(u,1,paste,collapse="")[mu],sep="")
#      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
#      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
#      u <- u[!mu,]
#     }
#    }
  if(contrast){
   u <- u[-1,]
  }
  if(all(base!=0)){
   u <- u[-base,]
  }
  termnumber<-1+length(m$terms)
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui), so that's what
  #  input component 1 equals
  m$terms[[termnumber]] <- list(name="hammingmix_constant", soname="ergm",
    inputs=c(1, 1, nrow(xm)*2+length(nodecov)+1,
            nrow(xm),as.integer(xm), nodecov),
            dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("hammingmix.constant",attrname, sep="."))
  m
}

#########################################################
InitErgm.hammingmix<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingmix", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingmix", arglist=arglist,
    varnames = c("attrname","x","base","contrast"),
    vartypes = c("character","matrixnetwork","numeric","logical"),
    defaultvalues = list(NULL,nw,0,FALSE),
    required = c(TRUE,FALSE,FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  x<-a$x
  base<-a$base
  contrast<-a$contrast
  drop<-a$drop
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
    stop("hammingmix() requires an edgelist", call.=FALSE)
  }
    nodecov <- get.node.attr(nw, attrname, "hammingmix")
    mixmat <- mixingmatrix(nw,attrname)$mat
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
  if(all(base!=0)){
   u <- u[-base,]
  }
  termnumber<-1+length(m$terms)
  #  Number of input parameters before covariates equals twice the number
  #  of used matrix cells, namely 2*length(uui), so that's what
  #  input component 1 equals
  m$terms[[termnumber]] <- list(name="hammingmix", soname="ergm",
    inputs=c(nrow(u), nrow(u), nrow(xm)*2+length(nodecov)+length(u)+1,
            nrow(xm),as.integer(xm), u[,1], u[,2],nodecov),
            dependence=FALSE)
  m$coef.names<-c(m$coef.names,
       paste("hammingmix",attrname, apply(matrix(namescov[u],ncol=2),1,paste,collapse="."), sep="."))
  m
}

###################################### InitErgm TERMS:  I
#########################################################
InitErgm.idegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("idegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("idegree", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "idegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to idegree() has only one value", call.=FALSE)
  }
  if(!is.null(attrname) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      idegreeattr <- summary(
       as.formula(paste('nw ~ idegree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(idegreeattr)){
        dropterms <- paste("ideg", du[1,idegreeattr], ".", attrname,
                           u[du[2,idegreeattr]], sep="")
      cat(" ")
        cat("Warning: These idegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        du <- matrix(du[,!idegreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(drop){
      tmp <- paste("c(",paste(d,collapse=","),")",sep="")
      if(!homophily) {
        midegree <- summary(as.formula(paste('nw ~ idegree(',tmp,')',
                                            sep="")), drop=FALSE) == 0
      } else {
        midegree <- summary(as.formula(paste('nw ~ idegree(',tmp,',"',attrname,
                                                         '", TRUE)', sep="")), 
                                             drop = FALSE) == 0
      }
      if(any(midegree)){
      cat(" ")
        cat("Warning: These idegree terms have extreme counts and will be dropped:\n")
        cat(d[midegree], "\n", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        d <- d[!midegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  termnumber<-1+length(m$terms)
  if(is.null(attrname)) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="idegree", soname="ergm",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("idegree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="idegree_w_homophily", soname="ergm",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("ideg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="idegree_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_idegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("ideg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.intransitive<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("intransitive", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("intransitive", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    nintransitive <- summary(as.formula('nw ~ intransitive'), drop=FALSE)
    if(nintransitive==0){
      cat(" ")
      cat(paste("Warning: There are no intransitive triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
 #    cat(paste("To avoid degeneracy the 'intransitive' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="intransitive", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"intransitive")
  m
}

#########################################################
InitErgm.isolates<-function(nw, m, arglist, drop=TRUE, ...) {
# ergm.checkdirected("isolates", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("isolates", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    mdsp <- summary(as.formula('nw ~ isolates'), drop=FALSE)
    if(mdsp==0){
      cat(" ")
      cat(paste("Warning: There are no isolates;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="isolates", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"isolates")
  m
}

#########################################################
InitErgm.istar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("istar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("istar", arglist,
    varnames = c("k", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  k <- a$k
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "istar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
#     Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to istar() has only one value", call.=FALSE)
    if(drop){
      istarattr <- paste("c(",paste(k,collapse=","),")",sep="")
      istarattr <- summary(as.formula(paste('nw ~ istar(',istarattr,',"',
                                            attrname,'")',sep="")),
                           drop=FALSE) == 0
      if(any(istarattr)){
        dropterms <- paste(paste("istar",attrname,sep="."),k[istarattr],sep="")
      cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        k <- k[!istarattr] 
      }
    }
  }else{
    if(drop){
      mistar <- paste("c(",paste(k,collapse=","),")",sep="")
      mistar <- summary(as.formula(paste('nw ~ istar(',mistar,')',sep="")),
                        drop=FALSE) == 0
      if(any(mistar)){
      cat(" ")
        cat(paste("Warning: There are no order", k[mistar],"stars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("istar", k[mistar],sep="")
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        k <- k[!mistar] 
      }
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="istar", soname="ergm",
                                  inputs=c(lk, lk, lk+length(nodecov),
                                    k, nodecov))
    m$coef.names<-c(m$coef.names,paste("istar",k,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="istar", soname="ergm",
                                  inputs=c(0, lk, lk, k))
    m$coef.names<-c(m$coef.names,paste("istar",k,sep=""))
  }
  m
}

###################################### InitErgm TERMS:  K
#########################################################
InitErgm.kstar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("kstar", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("kstar", arglist,
    varnames = c("k", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  k<-a$k;attrname<-a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "kstar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
#    Recode to numeric if necessary
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to kstar() has only one value", call.=FALSE)
    if(drop){
      kstarattr <- paste("c(",paste(k,collapse=","),")",sep="")
      kstarattr <- summary(as.formula(paste('nw ~ kstar(',kstarattr,
                                            ',"',attrname,'")',sep="")),
                           drop=FALSE) == 0
      if(any(kstarattr)){
        dropterms <- paste(paste("kstar",attrname,sep="."),k[kstarattr],sep="")
      cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        k <- k[!kstarattr] 
      }
    }
  }else{
    if(drop){
      mkstar <- paste("c(",paste(k,collapse=","),")",sep="")
      mkstar <- summary(as.formula(paste('nw ~ kstar(',mkstar,')',sep="")),
                        drop=FALSE) == 0
      if(any(mkstar)){
      cat(" ")
        cat(paste("Warning: There are no order", k[mkstar],"stars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("kstar", k[mkstar],sep="")
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        k <- k[!mkstar] 
      }
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="kstar", soname="ergm",
                                  inputs=c(lk, lk, lk+length(nodecov),
                                    k, nodecov))
    m$coef.names<-c(m$coef.names,paste("kstar",k,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="kstar", soname="ergm",
                                  inputs=c(0, lk, lk, k))
    m$coef.names<-c(m$coef.names,paste("kstar",k,sep=""))
  }
  m
}

###################################### InitErgm TERMS:  L
#########################################################
InitErgm.localtriangle<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("localtriangle", arglist,
    varnames = c("x", "attrname"),
    vartypes = c("matrixnetwork", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  x<-a$x;attrname<-a$attrname
  if(is.network(x))
    xm<-as.matrix.network(x, matrix.type="adjacency", attrname)
  else if(is.character(x))
    xm<-as.matrix.network(nw, matrix.type="adjacency", x)
  else
    xm<-as.matrix(x)
  termnumber <- 1 + length(m$terms)
  m$terms[[termnumber]] <- list(name = "localtriangle", soname="ergm", 
                                inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
                                  NROW(xm), as.double(xm)))
  if(!is.null(attrname))
    cn<-paste("localtriangle", as.character(sys.call(0)[[4]][2]), 
              as.character(sys.call(0)[[5]]), sep = ".")
  else
    cn<-paste("localtriangle", as.character(sys.call(0)[[4]][2]), sep = ".")
  m$coef.names <- c(m$coef.names, cn)
  m
}

###################################### InitErgm TERMS:  M
#########################################################
InitErgm.m2star<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("m2star", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("m2star", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
   degrees <- as.matrix.network.edgelist(nw)
   if(all(is.na(match(degrees[,1],degrees[,2])))){
      cat(" ")
    cat(paste("Warning: The are no mixed 2-stars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#   cat(paste("To avoid degeneracy the 'm2star' term has been dropped.\n"))
    return(m)
   }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(
       name="m2star", 
       soname="ergm",
       inputs=c(0,1,0),
       dependence=TRUE)
  m$coef.names<-c(m$coef.names,"m2star")
  m
}

#########################################################
InitErgm.meandeg<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("meandeg", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="meandeg", soname="ergm",
                                inputs=c(0, 1, 0),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,"meandeg")
  m
}

#########################################################
InitErgm.mutual<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("mutual", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("mutual", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    nmutual <- summary(as.formula('nw ~ mutual'), drop=FALSE)
    if(nmutual==0){
      cat(" ")
      cat(paste("Warning: There are no mutual ties;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'mutual' term has been dropped.\n"))
      return(m)
    }
    if(nmutual==network.dyadcount(nw)){
      cat(" ")
      cat(paste("Warning: All dyads have mutual ties!\n",
                 " the corresponding coefficient has been fixed at its MLE of infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'mutual' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="mutual", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"mutual")
  m
}

###################################### InitErgm TERMS:  N
#########################################################
InitErgm.nearsimmelian<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("nearsimmelian", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("nearsimmelian", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    nsimmelian <- summary(as.formula('nw ~ nearsimmelian'), drop=FALSE)
    if(nsimmelian==0){
      cat(" ")
      cat(paste("Warning: There are no nearsimmelian triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'nearsimmelian' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelian==network.dyadcount(nw)*network.size(nw)*0.5){
      cat(" ")
      cat(paste("Warning: All dyads have nearsimmelian triads!\n",
                 " the corresponding coefficient has been fixed at its MLE of infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'nearsimmelian' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nearsimmelian", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"nearsimmelian")
  m
}

#########################################################
InitErgm.nodecov<-InitErgm.nodemain<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("nodecov", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("nodecov",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "nodecov", numeric=TRUE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodecov", soname="ergm",
                                inputs=c(0,1,length(nodecov),nodecov),
                                dependence=FALSE)
  m
}

#########################################################
InitErgm.nodefactor<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("nodefactor", arglist,
    varnames = c("attrname", "base"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, 1),
    required = c(TRUE, FALSE))
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nodecov <- get.node.attr(nw, attrname, "nodefactor")
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
      dropterms <- paste(paste("nodefactor",attrname,sep="."),u[nfc==0],sep="")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodefactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  
  m$terms[[termnumber]] <- list(name="nodefactor", soname="ergm",
                                inputs=c(lu-length(base), 
                                         lu-length(base), 
                                         lu-length(base)+length(nodecov),
                                         ui[-base], nodecov), dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("nodefactor",
                                      attrname, paste(u[-base]), sep="."))
  m
}

#########################################################
InitErgm.nodeicov<-function (nw, m, arglist, ...) {
  ergm.checkdirected("nodeicov", is.directed(nw), requirement=TRUE,
                     extramessage="See 'nodecov'.")
  a <- ergm.checkargs("nodeicov", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("nodeicov",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "nodeicov", numeric=TRUE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodeicov", soname="ergm",
                                inputs=c(0,1,length(nodecov),nodecov))
  m
}

#########################################################
InitErgm.nodeifactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("nodeifactor", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("nodeifactor", arglist,
    varnames = c("attrname", "base"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, 1),
    required = c(TRUE, FALSE))
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nodecov <- get.node.attr(nw, attrname, "nodeifactor")
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
      dropterms <- paste(paste("nodeifactor",attrname,sep="."),u[nfc==0],sep="")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodeifactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)  
  m$terms[[termnumber]] <- list(name="nodeifactor", soname="ergm",
                                inputs=c(lu-length(base), 
                                         lu-length(base), 
                                         lu-length(base)+length(nodecov),
                                         ui[-base], nodecov), dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("nodeifactor",
                                      attrname, paste(u[-base]), sep="."))
  m
}

#########################################################
InitErgm.nodematch<-InitErgm.match<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("nodematch", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(TRUE, FALSE))
  attach(a)
  attrname<-a$attrname
  nodecov <- get.node.attr(nw, attrname, "nodematch")
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
#   Recode to numeric if necessary
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if (length(u)==1)
    stop ("Argument to nodematch() has only one value", call.=FALSE)
  if(drop){
    mixmat <- mixingmatrix(nw,attrname)$mat
    ematch  <- diag(mixmat)
    if(diff){
      offematch <- apply(mixmat,1,sum)+apply(mixmat,2,sum)-2*ematch
#     Diagonals or off-diagonals are zero
      mu <- ematch==0 | offematch==0
      mu[is.na(mu)] <- FALSE
      if(any(mu)){
        dropterms <- paste(paste("nodematch",attrname,sep="."),u[mu],sep="")
        cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        u <- u[!mu] 
        ui <- ui[!mu] 
      }
    }else{
      offematch <- sum(mixmat)-sum(ematch)
#     Diagonals or off-diagonals are zero
      mu <- sum(ematch)==0 | offematch==0
      mu[is.na(mu)] <- FALSE
      if(mu){
      cat(" ")
        cat(paste("Warning: The number of matching dyads is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("nodematch",attrname,sep=".")
#       cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                 "has been dropped.\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)  
  if (!diff) {
    m$terms[[termnumber]] <- list(name="nodematch", soname="ergm",
                                  inputs=c(0,1,length(nodecov),nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("nodematch",attrname,sep="."))
  } else {
        #  Number of input parameters before covariates equals number of
        #  unique elements in nodecov, namely length(u), so that's what
        #  input element 1 equals
    m$terms[[termnumber]] <- list(name="nodematch", soname="ergm",
                                  inputs=c(length(ui), length(ui),
                                    length(ui)+length(nodecov),
                                    ui, nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("nodematch",
                                       attrname, u, sep="."))
  }
  m
}

#########################################################
InitErgm.nodemix<-InitErgm.mix<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("nodemix", arglist,
    varnames = c("attrname","contrast"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,FALSE),
    required = c(TRUE,FALSE))
  attach(a)
  attrname<-a$attrname
  contrast<-a$contrast
  if(is.bipartite(nw)){
  # So two-mode
    if (is.directed(nw)){ 
      cat(" ")
      cat("Warning!  Bipartite networks are currently\n",
          "automatically treated as undirected\n")
    }
    #  So undirected network storage but directed mixing
    nodecov <- get.node.attr(nw, attrname, "mix")
    nb1 <- get.network.attribute(nw, "bipartite")       
    #  Recode nodecov to numeric (but retain original sorted names in "namescov")
    b1namescov <- sort(unique(nodecov[1:nb1]))
    b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
    namescov <- c(b1namescov, b2namescov)
    b1nodecov <- match(nodecov[1:nb1],b1namescov)
    mixmat <- mixingmatrix(nw,attrname)$mat
    nodecov <- c(b1nodecov, 
     match(nodecov[(1+nb1):network.size(nw)],b2namescov)+nrow(mixmat))
    if (length(nodecov)==1)
        stop ("Argument to mix() has only one value", call.=FALSE)
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)+nrow(mixmat)))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
    #  Check for degeneracy
    if(drop){
     mu <- as.vector(mixmat)==0
     mu[is.na(mu)] <- FALSE
     if(any(mu)){
      dropterms <- paste(paste("mix",attrname,sep="."),
        apply(u,1,paste,collapse="")[mu],sep="")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#         "have been dropped.\n"))
      u <- u[!mu,]
     }
    }
    if(contrast){
     u <- u[-1,]
    }
    termnumber<-1+length(m$terms)
    #  Number of input parameters before covariates equals twice the number
    #  of used matrix cells, namely 2*length(uui), so that's what
    #  input element 1 equals
    m$terms[[termnumber]] <- list(name="mix", soname="ergm",
      inputs=c(NROW(u), NROW(u), length(nodecov)+length(u), u[,1], u[,2],nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,
       paste("mix",attrname, apply(matrix(namescov[u],ncol=2),1,paste,collapse="."), sep="."))
  }else{
#
# So one mode, but could be directed or undirected
#
    nodecov <- get.node.attr(nw, attrname, "nodemix")
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
    if (length(u)==1)   
      stop ("Argument to nodemix() has only one value", call.=FALSE)
    if(drop){
      mixmat <- mixingmatrix(nw,attrname)$mat
      if(is.directed(nw))       #If directed, accumulate in upper triangle
        mixmat[upper.tri(mixmat)] <- t(mixmat)[upper.tri(mixmat)]
      maxmat <- ucount %o% ucount - diag(length(ucount))
      c1mat <- diag(ucount,length(ui),length(ui))==1
      mixvec <- mixmat[upper.tri(mixmat,diag=TRUE)]
  #   Check for extreme cells
      mu <- (mixvec==0) | (mixvec==(maxmat[upper.tri(maxmat,diag=TRUE)])) |  (c1mat[upper.tri(c1mat,diag=TRUE)])
      mu[is.na(mu)] <- FALSE
      if(any(mu)){
        dropterms <- paste(paste("nodemix",attrname,sep="."),uun[mu],sep=".")
      cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#       cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                 "have been dropped.\n"))
        if (sum(!mu)<=1){
          stop ("The attribute to nodemix() must have more than one value", call.=FALSE)
        }
        uun <- uun[!mu]
        uui <- uui[!mu]
        urm <- urm[!mu]
        ucm <- ucm[!mu]
      }
    }
    if(contrast){u <- u[-1]}
    termnumber<-1+length(m$terms)
    #  Number of input parameters before covariates equals twice the number
    #  of used matrix cells, namely 2*length(uui), so that's what
    #  input element 1 equals
    m$terms[[termnumber]] <- list(name="nodemix", soname="ergm",
                                  inputs=c(2*length(uui), length(uui),
                                    2*length(uui)+length(nodecov),
                                    urm, ucm, nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("mix",attrname, uun, sep="."))
  }
  m
}

#########################################################
InitErgm.nodeocov<-function (nw, m, arglist, ...) {
  ergm.checkdirected("nodeocov", is.directed(nw), requirement=TRUE,
                     extramessage="See 'nodecov'.")
  a <- ergm.checkargs("nodeocov", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("nodeocov",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "nodeocov", numeric=TRUE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodeocov", soname="ergm",
                                inputs=c(0,1,length(nodecov),nodecov))
  m
}

#########################################################
InitErgm.nodeofactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("nodeofactor", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("nodeofactor", arglist,
    varnames = c("attrname", "base"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, 1),
    required = c(TRUE, FALSE))
  attach(a)
  attrname<-a$attrname
  base <- a$base
  nodecov <- get.node.attr(nw, attrname, "nodeofactor")
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
      dropterms <- paste(paste("nodeofactor",attrname,sep="."),u[nfc==0],sep="")
      cat(" ")
      cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodeofactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)  
  m$terms[[termnumber]] <- list(name="nodeofactor", soname="ergm",
                                inputs=c(lu-length(base), 
                                         lu-length(base), 
                                         lu-length(base)+length(nodecov),
                                         ui[-base], nodecov), dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("nodeofactor",
                                      attrname, paste(u[-base]), sep="."))
  m
}

###################################### InitErgm TERMS:  O
#########################################################
InitErgm.odegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("odegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("odegree", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "odegree")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to odegree() has only one value", call.=FALSE)
  }
  if(!is.null(attrname) && !homophily) {
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      odegreeattr <- summary(
       as.formula(paste('nw ~ odegree(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(odegreeattr)){
        dropterms <- paste("odeg", du[1,odegreeattr], ".", attrname,
                           u[du[2,odegreeattr]], sep="")
      cat(" ")
        cat("Warning: These odegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        du <- matrix(du[,!odegreeattr], nrow=2)
      }      
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- sum(nodecov==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(drop){
      tmp <- paste("c(",paste(d,collapse=","),")",sep="")
      if(!homophily) {
        modegree <- summary(as.formula(paste('nw ~ odegree(',tmp,')',
                                            sep="")), drop=FALSE) == 0
      } else {
        modegree <- summary(as.formula(paste('nw ~ odegree(',tmp,',"',attrname,
                                                         '", TRUE)', sep="")), 
                                             drop = FALSE) == 0
      }
      if(any(modegree)){
      cat(" ")
        cat("Warning: These odegree terms have extreme counts and will be dropped:\n")
        cat(d[modegree], "\n", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        d <- d[!modegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- network.size(nw)
    }
  }
  termnumber<-1+length(m$terms)
  if(is.null(attrname)) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="odegree", soname="ergm",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("odegree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="odegree_w_homophily", soname="ergm",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    # See comment in d_odegree_w_homophily function
    m$coef.names<-c(m$coef.names,paste("odeg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="odegree_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_odegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("odeg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.ostar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ostar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ostar", arglist,
    varnames = c("k", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  k<-a$k
  attrname <- a$attrname
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "ostar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    if (length(u)==1)
      stop ("Attribute given to ostar() has only one value", call.=FALSE)
    if(drop){
      ostarattr <- paste("c(",paste(k,collapse=","),")",sep="")
      ostarattr <- summary(as.formula(paste('nw ~ ostar(',ostarattr,',"',
                                            attrname,'")',sep="")),
                           drop=FALSE) == 0
      if(any(ostarattr)){
        dropterms <- paste(paste("ostar",attrname,sep="."),k[ostarattr],sep="")
        cat(" ")
        cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        k <- k[!ostarattr]
      }
    }
  }else{
    if(drop){
      mostar <- paste("c(",paste(k,collapse=","),")",sep="")
      mostar <- summary(as.formula(paste('nw ~ ostar(',mostar,')',sep="")),
                        drop=FALSE) == 0
      if(any(mostar)){
      cat(" ")
        cat(paste("Warning: There are no order", k[mostar],"stars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("ostar", k[mostar],sep="")
        k <- k[!mostar] 
      }
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="ostar", soname="ergm",
                                  inputs=c(lk, lk, lk+length(nodecov),
                                           k, nodecov))
    m$coef.names<-c(m$coef.names,paste("ostar",k,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="ostar", soname="ergm",
                                  inputs=c(0, lk, lk, k))
    m$coef.names<-c(m$coef.names,paste("ostar",k,sep=""))
  }
  m
}

###################################### InitErgm TERMS:  R
#########################################################
InitErgm.receiver<-function(nw, m, arglist, drop=FALSE, ...) {
  ergm.checkdirected("receiver", is.directed(nw), requirement=TRUE,
                     extramessage="See 'sociality'.")
  a <- ergm.checkargs("receiver", arglist,
    varnames = "drop",
    vartypes = "logical",
    defaultvalues = list(FALSE),
    required = FALSE)
  drop<-a$drop
  d <- 2:network.size(nw)
  if(drop){
    degrees <- as.numeric(names(table(as.matrix.network.edgelist(nw)[,2])))
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no in ties for the vertex", 
                d[is.na(mdegrees)],";\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("receiver", d[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="receiver", soname="ergm",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("receiver",d,sep=""))
  m
}

###################################### InitErgm TERMS:  S
#########################################################
InitErgm.sender<-function(nw, m, arglist, drop=FALSE, ...) {
  ergm.checkdirected("sender", is.directed(nw), requirement=TRUE,
                     extramessage="See 'sociality'.")
  a <- ergm.checkargs("sender", arglist,
    varnames = "drop",
    vartypes = "logical",
    defaultvalues = list(FALSE),
    required = FALSE)
  drop<-a$drop
  d <- 2:network.size(nw)
  if(drop){
    degrees <- as.numeric(names(table(as.matrix.network.edgelist(nw)[,1])))
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no out ties for the vertex", 
                d[is.na(mdegrees)],";\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("sender", d[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="sender", soname="ergm",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("sender",d,sep=""))
  m
}

#########################################################
InitErgm.simmelian<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("simmelian", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("simmelian", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    nsimmelian <- summary(as.formula('nw ~ simmelian'), drop=FALSE)
    if(nsimmelian==0){
      cat(" ")
      cat(paste("Warning: There are no simmelian triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'simmelian' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelian==network.edgecount(nw)*network.size*0.5){
      cat(" ")
      cat(paste("Warning: All triads are simmelian!\n",
                 " The corresponding coefficient has been fixed at its MLE of infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'simmelian' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="simmelian", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"simmelian")
  m
}

#########################################################
InitErgm.simmelianties<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("simmelianties", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("simmelianties", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  if(drop){
    nsimmelianties <- summary(as.formula('nw ~ simmelianties'), drop=FALSE)
    if(nsimmelianties==0){
      cat(" ")
      cat(paste("Warning: There are no simmelianties ties;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'simmelianties' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelianties==network.edgecount(nw)){
      cat(" ")
      cat(paste("Warning: All ties have simmelianties ties!\n",
                 " the corresponding coefficient has been fixed at its MLE of infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'simmelianties' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="simmelianties", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"simmelianties")
  m
}

#########################################################
InitErgm.smalldiff<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("smalldiff", arglist,
    varnames = c("attrname", "cutoff"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, TRUE))
  attach(a)
  cutoff <- a$cutoff
  attrname <- a$attrname
  if (length(cutoff)>1)
    stop("cutoff for smalldiff() must be a scalar.", call.=FALSE)
  m$coef.names<-c(m$coef.names,paste("smalldiff.",
                                     attrname, cutoff, sep=""))
  nodecov <- get.node.attr(nw, attrname, "smalldiff", numeric=TRUE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="smalldiff", soname="ergm",
                                inputs=c(1, 1, 1+length(nodecov),
                                  cutoff, nodecov), dependence=FALSE)
  m
}

#########################################################
InitErgm.sociality<-function(nw, m, arglist, drop=FALSE, ...) {
  ergm.checkdirected("sociality", is.directed(nw), requirement=FALSE,
                     extramessage = "See 'sender' and 'receiver'.")
  a <- ergm.checkargs("sociality", arglist,
    varnames = c("attrname","drop"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,FALSE),
    required = c(FALSE,FALSE))
  attach(a)
  attrname<-a$attrname
  drop<-a$drop
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "sociality")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to sociality() has only one value", call.=FALSE)
  }
  d <- 1:network.size(nw)
  if(drop){
    if(is.null(attrname)){
      centattr <- summary(nw ~ sociality, drop=FALSE) == 0
    }else{
      centattr <- summary(as.formula(paste('nw ~ sociality(','"',attrname,
                                           '")',sep="")),
                          drop=FALSE) == 0
    }
    if(any(centattr)){
      cat(" ")
      cat(paste("Warning: There are no",attrname," ties for the vertex", 
                d[centattr],";\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("sociality", d[centattr],sep="")
#     cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- d[!centattr] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="sociality", soname="ergm",
                                  inputs=c(0, ld, ld+length(nodecov),
                                    d, nodecov))
    m$coef.names<-c(m$coef.names,paste("sociality",d,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="sociality", soname="ergm",
                                  inputs=c(0, ld, ld, d))
    m$coef.names<-c(m$coef.names,paste("sociality",d,sep=""))
  }
  m
}

###################################### InitErgm TERMS:  T
#########################################################
InitErgm.transitive<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("transitive", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("transitive", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    ntransitive <- summary(as.formula('nw ~ transitive'), drop=FALSE)
    if(ntransitive==0){
      cat(" ")
      cat(paste("Warning: There are no transitive triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'transitive' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="transitive", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"transitive")
  m
}

#########################################################
InitErgm.triadcensus<-function (nw, m, arglist, drop=FALSE, ...) {
  a=ergm.checkargs("triadcensus", arglist,
    varnames = c("d","drop"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  d<-a$d
# drop<-a$drop
  detach(a)

  if(is.directed(nw)){
   tcn <- c("003","012", "102", "021D", "021U", "021C", "111D",
            "111U", "030T", "030C", "201", "120D", "120U", "120C", "210", "300")
   if(is.null(d)){d <- 1:15}
   if(is.character(d)){d <- match(d, tcn)-1}
  }else{
#  Undirected
   tcn <- c("0", "1", "2", "3")
   if(is.null(d)){d <- 1:3}
   if(is.character(d)){d <- match(d, tcn)-1}
  }
  if(drop){
    mdegree <- paste("c(",paste(d,collapse=","),")",sep="")
    mdegree <- summary(
     as.formula(paste('nw ~ triadcensus(',mdegree,')',sep="")),
     drop=FALSE) == 0
    if(any(mdegree)){
     dropterms <- tcn[d[mdegree]]
     cat(" ")
     cat(paste("Warning: There are no triads of type", tcn[d[mdegree]],".\n"))
     cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
     d <- d[!mdegree]
    }
  }
  d <- d + 1
  lengthd<-length(d)
  if(lengthd==0){return(m)}
  termnumber<-1+length(m$terms)
# No covariates here, so input component 1 is arbitrary
  m$terms[[termnumber]] <- list(name="triadcensus", soname="ergm",
                                      inputs=c(0, lengthd, lengthd, d),
                                      dependence=TRUE)
  m$coef.names<-c(m$coef.names, paste("triadcensus",tcn,sep=".")[d])
  m
}

#########################################################
InitErgm.triangle<-InitErgm.triangles<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("triangle", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  attrname <- a$attrname
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "triangle")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to triangle() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ triangle(','"', attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("triangle",attrname,sep="."),
                             u[triattr],sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("triangle",attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="triangle", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("triangle",attrname,sep="."))
    } else {
      m$terms[[termnumber]] <- list(name="triangle", soname="ergm",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("triangle",
                                         attrname, u, sep="."))
    }
  }else{
    m$terms[[termnumber]] <- list(name="triangle", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"triangle")
  }
  m
}

#########################################################
InitErgm.tripercent<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("tripercent", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("tripercent", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  attrname <- a$attrname
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "tripercent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to tripercent() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ tripercent(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(any(triattr)){
        if(diff){
          dropterms <- paste(paste("tripercent",attrname,sep="."),
                             u[triattr],sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }else{
          dropterms <- paste(paste("tripercent",attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="tripercent", soname="ergm",
                                    inputs=c(1, 1, 1+length(nodecov),
                                      1, nodecov))
      m$coef.names<-c(m$coef.names,paste("tripercent",attrname,sep="."))
    } else {
      m$terms[[termnumber]] <- list(name="tripercent", soname="ergm",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("tripercent",
                                         attrname, u, sep="."))
    }
  }else{
    m$terms[[termnumber]] <- list(name="tripercent", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"tripercent")
  }
  m
}

#########################################################
InitErgm.ttriple<-InitErgm.ttriad<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ttriple", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ttriple", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  attrname <- a$attrname
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "ttriple")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ttriple() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ ttriple(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("ttriple",attrname,sep="."),
                             u[triattr],sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("ttriple",attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="ttriple", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("ttriple",attrname,sep="."))
     } else {
       m$terms[[termnumber]] <- list(name="ttriple", soname="ergm",
                                     inputs=c(length(ui), length(ui),
                                       length(ui)+length(nodecov),
                                       ui, nodecov))
       m$coef.names<-c(m$coef.names,paste("ttriple",
                                          attrname, u, sep="."))
     }
  }else{
    m$terms[[termnumber]] <- list(name="ttriple", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"ttriple")
  }
  m
}

#########################################################
InitErgm.twopath<-function(nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("twopath", arglist,
     varnames = NULL,
     vartypes = NULL,
     defaultvalues = list(),
     required = NULL)
  if(is.directed(nw)){
   if(drop){
    degrees <- as.matrix.network.edgelist(nw)
    if(all(is.na(match(degrees[,1],degrees[,2])))){
      cat(" ")
     cat(paste("Warning: The are no two-paths;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#    cat(paste("To avoid degeneracy the 'twopath' term has been dropped.\n"))
     return(m)
    }
   }
   termnumber<-1+length(m$terms)
   m$terms[[termnumber]] <- list(
        name="m2star", 
        soname="ergm",
        inputs=c(0,1,0),
        dependence=TRUE)
  }else{
   k<-2
   if(drop){
    mkstar <- paste("c(",paste(k,collapse=","),")",sep="")
    mkstar <- summary(as.formula(paste('nw ~ kstar(',mkstar,')',sep="")),
                      drop=FALSE) == 0
    if(any(mkstar)){
      cat(" ")
      cat(paste("Warning: There are no two paths;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the twopath term has been dropped.\n"))
      return(m)
    }
   }
   lk<-length(k)
   if(lk==0){return(m)}
   termnumber<-1+length(m$terms)
   m$terms[[termnumber]] <- list(name="kstar", soname="ergm",
                                 inputs=c(0, lk, lk, k))
  }
  m$coef.names<-c(m$coef.names,"twopath")
  m
}






