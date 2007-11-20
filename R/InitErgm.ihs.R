#  See InitErgm.R for a general explanation 
#  of InitErgm functions

###################################### InitErgm TERMS:  A
#########################################################
InitErgm.altistar<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("altistar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("altistar", arglist,
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
    m$terms[[termnumber]] <- list(name="idegree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altistar=NULL,
                                    altistar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altistar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altistar", soname="ergm",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altistar")
  }
  m
}

#########################################################
InitErgm.altostar<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("altostar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("altostar", arglist,
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
    m$terms[[termnumber]] <- list(name="odegree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altostar=NULL,
                                    altostar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altostar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altostar", soname="ergm",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altostar")
  }
  m
}

###################################### InitErgm TERMS:  B
#########################################################
InitErgm.balance<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("balance", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  attrname <- a$attrname
  detach(a)
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
# ergm.checkdirected("concurrent", is.directed(nw), requirement=TRUE)
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
        cat(dropterms, "\n", fill=T)
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
        cat(paste("Warning: There are no concurrent actors;\n",
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

###################################### InitErgm TERMS:  D
#########################################################
InitErgm.degreep<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("degreep", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("degreep", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname; homophily <- a$homophily
  emptynwstats<-NULL
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "degreep")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
         stop ("Attribute given to degreep() has only one value", call.=FALSE)
  }
  if(!is.null(attrname) && !homophily) {
    # Combine degreep and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      degreepattr <- summary(
       as.formula(paste('nw ~ degreep(',tmp,',"',attrname,'")',sep="")),
       drop=FALSE) == 0
      if(any(degreepattr)){
        dropterms <- paste("pdeg", du[1,degreepattr], ".", attrname,
                           u[du[2,degreepattr]], sep="")
      cat(" ")
        cat("Warning: These degreep terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        du <- matrix(du[,!degreepattr], nrow=2)
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
        mdegreep <- summary(as.formula(paste('nw ~ degreep(',tmp,')',
                                            sep="")), drop=FALSE) == 0
      } else {
        mdegreep <- summary(as.formula(paste('nw ~ degreep(',tmp,',"',attrname,
                                                         '", TRUE)', sep="")), 
                                             drop = FALSE) == 0
      }
      if(any(mdegreep)){
      cat(" ")
        cat("Warning: These degreep terms have extreme counts and will be dropped:\n")
        cat(d[mdegreep], "\n", fill=T)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        d <- d[!mdegreep] 
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
    m$terms[[termnumber]] <- list(name="degreep", soname="ergm",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("degreep",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="degreep_w_homophily", soname="ergm",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    # See comment in d_degreep_w_homophily function
    m$coef.names<-c(m$coef.names,paste("pdeg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="degreep_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_degreep_by_attr function
    m$coef.names<-c(m$coef.names, paste("pdeg", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }
  if (!is.null(emptynwstats))
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.duration<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("duration", arglist,
    varnames = c("form", "dissolve", "x"),
    vartypes = c("matrixnetwork", "matrixnetwork", "matrixnetwork"),
    defaultvalues = list(NULL, NULL, NULL),
    required = c(TRUE, TRUE, FALSE))
  attach(a)
  x<-a$x;form<-a$form;dissolve<-a$dissolve
  m$coef.names<-c(m$coef.names, paste("duration.",x,sep=""))
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
  if (is.null(xm) || NCOL(xm)!=2){
    stop("duration requires the edgelist of the base network")
  }
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
  if (is.null(formm) || dim(formm)!=c(network.size(nw),network.size(nw))){
    stop("duration requires a matrix of formation rates")
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
  if (is.null(dissolvem) || dim(dissolvem)!=
      c(network.size(nw),network.size(nw))){
    stop("duration requires a matrix of dissolution rates")
  }
  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1 (although in this case, input element 1
# is actually arbitrary since d_duration ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "duration", soname="ergm",
                                inputs = c(1, 1, 
                                  NROW(xm)*2+2*NROW(formm)^2, NROW(xm),
                                  as.double(c(xm, formm, dissolvem))))
  m
}

###################################### InitErgm TERMS:  E

###################################### InitErgm TERMS:  F
#########################################################
InitErgm.factor<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("factor", arglist,
    varnames = c("factorname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  factorname<-a$factorname
  x <- get(factorname)
  #Coerce x to an adjacency matrix
  if(!is.factor(x)){
    stop ("Attribute given to factor() should be have class 'factor'", call.=FALSE)
  }
  xm <- model.matrix(as.formula("~ x"))[,-1]
  if(drop){
    mfactor <- summary(nw ~ factor(x), drop=FALSE)
    if(all(mfactor==0)){return(m)}
    if(any(mfactor==0)){
      cat(" ")
      cat(paste("Warning: There are no dyads with factor level",
          names(mfactor)[mfactor==0],";\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the terms",names(mfactor)[mfactor==0],
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      xm <- xm[,mfactor!=0] 
    }
  }

  #Update the term number
  termnumber <- 1 + length(m$terms)
  #Update the terms list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1 (although in this case, input element 1
# is actually arbitrary since d_factor ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "factor",  soname="ergm",
                                inputs = c(1, 3, 1+NROW(xm)*NCOL(xm),
                                  NCOL(xm), as.double(as.numeric(xm))),
                                dependence=FALSE)
  cn <- paste(factorname,substring(colnames(xm),first=2),sep=".")
  m$coef.names<-c(m$coef.names, cn)
  m
}

###################################### InitErgm TERMS:  G
#########################################################
InitErgm.geodegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("geodegree", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("geodegree", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map<- function(x,n,...){
      i <- 1:n
      x[1]*(exp(-(x[2])*i)-1)
    }
    gradient <- function(x,n,...) {
      i=1:n
      rbind(exp(-(x[2])*i)-1, -x[1]*i*exp(-x[2]*i))
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(geodegree=NULL,
                                    geodegree.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("geodegree#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="geodegree", soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"geodegree")
  }
  m
}

#########################################################
InitErgm.gwdegree706<-function(nw, m, arglist, initialfit=FALSE, ...) {
# Slight modification to the parameterization in gwdegree.
# ergm.checkdirected("gwdegree706", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwdegree706", arglist,
    varnames = c("decay","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  decay<-a$decay;fixed<-a$fixed
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      0.5 + x[1]*(exp(x[2])*(1-(1-exp(-x[2]))^i))
      # Note small change to add number of edges to usual gwdegree.
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(exp(x[2])*(1-(1-exp(-x[2]))^i),
            x[1]*(exp(x[2])-(1-exp(-x[2]))^{i-1}*(1+i-exp(-x[2])))
           )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegree=NULL,
                                    gwdegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegree#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegree706", soname="ergm",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwdegree706")
  }
  m
}

#########################################################
InitErgm.gwdegreealpha<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwdegreealpha", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwdegreealpha", arglist,
    varnames = c("alpha","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(0, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  alpha<-a$alpha;fixed<-a$fixed
  if(!initialfit && !fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    ld<-length(d)
    if(ld==0){return(m)}
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(i*exp(x[2]) + exp(2*x[2])*(((1-exp(-x[2]))^i)-1))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind((i*exp(x[2]) + exp(2*x[2])*(((1-exp(-x[2]))^i)-1)),
            x[1]*exp(x[2])*(i - 2*exp(x[2]) + 2*exp(x[2])*(1-exp(-x[2]))^i +
                            i*(1-exp(-x[2]))^(i-1)))
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegreealpha=NULL,
                                    gwdegree.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegreealpha#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegreealpha", soname="ergm",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwdegreealpha")
  }
  m
}

#########################################################
InitErgm.gwdegreelambda<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwdegreelambda", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("gwdegreelambda", arglist,
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
      x[1]*(x[2]*(1-(1-1/x[2])^i))
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(x[2]*(1-(1-1/x[2])^i),
            x[1]*(1 - (x[2]*x[2]-x[2]+i)*((1-1/x[2])^(i-1))/(x[2]*x[2]) )
           )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="degree", soname="ergm",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegreelambda=NULL,
                                    gwdegree.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegreelambda#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegreelambda", soname="ergm",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"gwdegree")
  }
  m
}

###################################### InitErgm TERMS:  H
#########################################################
InitErgm.heideriandynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("heideriandynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("heideriandynamic", arglist,
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
    xm<-get.network.attribute(nw,x)
  else
    xm<-as.matrix(x)

  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1
  m$terms[[termnumber]] <- list(name = "heideriandynamic", soname="ergm", 
                                inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
                                  NROW(xm), as.double(xm)),
                                dependence=FALSE)
  if(!is.null(attrname))
    cn<-paste("heideriandynamic", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
  else
    cn<-paste("heideriandynamic", as.character(sys.call(0)[[4]][2]), sep = ".")

  m$coef.names <- c(m$coef.names, cn)
  m
}

#########################################################
InitErgm.hiertriad<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("hiertriad", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("hiertriad", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    nhiertriad <- summary(as.formula('nw ~ hiertriad'), drop=FALSE)
    if(nhiertriad==0){
      cat(" ")
      cat(paste("Warning: There are no hiertriad ties;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'hiertriad' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="hiertriad", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"hiertriad")
  m
}

#########################################################
InitErgm.hiertriaddegree<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("hiertriaddegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("hiertriaddegree", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    nhiertriad <- summary(as.formula('nw ~ hiertriaddegree'), drop=FALSE)
    if(nhiertriad==0){
      cat(" ")
      cat(paste("Warning: There are no hiertriaddegree ties;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'hiertriaddegree' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="hiertriaddegree", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"hiertriaddegree")
  m
}

###################################### InitErgm TERMS:  I
#########################################################
InitErgm.icvar<-function(nw, m, arglist, ...) {
  ergm.checkdirected("icvar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("icvar", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
    termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="icvar", soname="ergm",
                                inputs=c(0, 1, 0),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,"icvar")
  m
}

#########################################################
InitErgm.idc<-function(nw, m, arglist, ...) {
  ergm.checkdirected("idc", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("idc", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="idc", soname="ergm",
                                inputs=c(0, 1, 0),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,"idc")
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
InitErgm.intransitivedynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("intransitivedynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("intransitivedynamic", arglist,
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

   termnumber <- 1 + length(m$terms)
#  There is 1 input parameter before the covariate vector, so input
#  element 1 is set to 1 (although in this case, input element 1
#  is actually arbitrary since d_intransitivedynamic ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "intransitivedynamic", soname="ergm", 
                                 inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
                                   NROW(xm), as.double(xm)),
                                 dependence=FALSE)
   if(!is.null(attrname))
     cn<-paste("intransitivedynamic", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
   else
     cn<-paste("intransitivedynamic", as.character(sys.call(0)[[4]][2]), sep = ".")

   m$coef.names <- c(m$coef.names, cn)
  m
}

#########################################################
InitErgm.intransitivity<-function (nw, m, arglist, drop=TRUE, ...) {
# ergm.checkdirected("intransitivity", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("intransitivity", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    nintransitive <- summary(as.formula('nw ~ intransitivity'), drop=FALSE)
    if(nintransitive==0){
      cat(" ")
      cat(paste("Warning: There are no intransitive triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'intransitivity' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="intransitivity", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"intransitivity")
  m
}

###################################### InitErgm TERMS:  K
#########################################################
InitErgm.kappa<-function(nw, m, arglist, ...) {
  ergm.checkdirected("kappa", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("kappa", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  if(is.bipartite(nw)){
   m$terms[[termnumber]] <- list(name="bkappa", soname="ergm",
                                inputs=c(0, 1, 0))
  }else{
   m$terms[[termnumber]] <- list(name="kappa", soname="ergm",
                                inputs=c(0, 1, 0))
  }
  m$coef.names<-c(m$coef.names,"kappa")
  m
}

###################################### InitErgm TERMS:  L

###################################### InitErgm TERMS:  M

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

###################################### InitErgm TERMS:  O

###################################### InitErgm TERMS:  R

###################################### InitErgm TERMS:  S
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
InitErgm.simmeliandynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("simmeliandynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("simmeliandynamic", arglist,
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

   termnumber <- 1 + length(m$terms)
#  There is 1 input parameter before the covariate vector, so input
#  element 1 is set to 1 (although in this case, input element 1
#  is actually arbitrary since d_simmeliandynamic ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "simmeliandynamic", soname="ergm", 
                                 inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
                                   NROW(xm), as.double(xm)),
                                 dependence=FALSE)
   if(!is.null(attrname))
     cn<-paste("simmeliandynamic", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
   else
     cn<-paste("simmeliandynamic", as.character(sys.call(0)[[4]][2]), sep = ".")

   m$coef.names <- c(m$coef.names, cn)
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
InitErgm.spatial<-function(nw, m, arglist, initialfit=FALSE, ...) {
  a <- ergm.checkargs("spatial", arglist,
    varnames = c("pb","alpha","gamma","d"),
    vartypes = c("numeric","numeric","numeric","numeric"),
    defaultvalues = list(10, -10, 0, NA),
    required = c(FALSE, FALSE, FALSE, TRUE))
  attach(a)
  if(is.directed(nw)){
    d<-as.vector(d)
    nstat<-network.size(nw)*(network.size(nw)-1)
  }else{
    d<-d[lower.tri(d)]
    nstat<-choose(network.size(nw),2)
  }
  if(!initialfit){
    map<- function(x,n,cv){
      -log(((1+exp(x[1]))*(1+exp(x[2])*cv)^exp(x[3]))/exp(x[1])-1)
    }
    gradient <- function(x,n,cv) {
      rbind(
        1/(1+exp(x[1])*(1-(1+exp(x[2])*cv)^(-exp(x[3])))),
        -cv*exp(x[2]+x[3])*(1+exp(x[1])) / ((1+exp(x[2])*cv)*(1+exp(x[1])*(1-(1+exp(x[2])*cv)^(-exp(x[3]))))),
        -(exp(x[3])*(1+exp(x[1]))*(1+exp(x[2])*cv)^exp(x[3])*log(1+exp(x[2])*cv)) / ((1+exp(x[1]))*(1+exp(x[2])*cv)^exp(x[3])-exp(x[1]))
      )
    }
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="berninhom", soname="ergm",
                                inputs=c(0, nstat, 0, 0),
                                dependence=TRUE,
                                params=list(spatial.pb=pb,
                                spatial.alpha=alpha,spatial.gamma=gamma),
                                map=map, gradient=gradient, 
                                eta.cov=d)
    m$coef.names<-c(m$coef.names,paste("spatial#",1:nstat,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="spatial", soname="ergm",
                                inputs=c(0, 1, 3+nstat, c(pb, alpha, gamma, d)),
                                dependence=TRUE)
    m$coef.names<-c(m$coef.names,"spatial.pb")
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
InitErgm.transitivedynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("transitivedynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("transitivedynamic", arglist,
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

   termnumber <- 1 + length(m$terms)
#  There is 1 input parameter before the covariate vector, so input
#  element 1 is set to 1 (although in this case, input element 1
#  is actually arbitrary since d_transitivedynamic ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "transitivedynamic", soname="ergm", 
                                 inputs = c(1, 1, 1+NROW(xm)*NROW(xm),
                                   NROW(xm), as.double(xm)),
                                 dependence=FALSE)
   if(!is.null(attrname))
     cn<-paste("transitivedynamic", as.character(sys.call(0)[[4]][2]), 
               as.character(attrname), sep = ".")
   else
     cn<-paste("transitivedynamic", as.character(sys.call(0)[[4]][2]), sep = ".")

   m$coef.names <- c(m$coef.names, cn)
  m
}

#########################################################
InitErgm.transitivity<-function (nw, m, arglist, drop=TRUE, ...) {
# ergm.checkdirected("transitivity", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("transitivity", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  if(drop){
    ntransitive <- summary(as.formula('nw ~ transitivity'), drop=FALSE)
    if(ntransitive==0){
      cat(" ")
      cat(paste("Warning: There are no transitive triads;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#     cat(paste("To avoid degeneracy the 'transitivity' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="transitivity", soname="ergm",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"transitivity")
  m
}




