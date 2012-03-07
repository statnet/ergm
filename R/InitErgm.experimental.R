#  File ergm/R/InitErgm.experimental.R
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
#########################################################
InitErgm.altistar<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("altistar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("altistar", arglist,
    varnames = c("lambda","fixed"),
    vartypes = c("numeric","logical"),
    defaultvalues = list(1, FALSE),
    required = c(FALSE, FALSE))
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
InitErgm.bounded.degree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.degree", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("bounded.degree", arglist,
    varnames = c("d","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  bound <- a$bound; d <- a$d
  if(drop){
    degrees <- as.numeric(names(table(table(as.edgelist(nw)))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no degree", d[is.na(mdegrees)],
                "vertices;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.degree", d[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  m$coef.names<-c(m$coef.names,paste("bounded.degree",d,sep=""))
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.degree() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  bound <- rep(bound, length=ld) 
  m$terms[[termnumber]] <- list(name="boundeddegree", soname="ergm",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  if (any(d==0)) {
    emptynwstats <- rep(0, ld)
    emptynwstats[d==0] <- min(bound[d==0],network.size(nw))
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  }
  m
}

#########################################################
InitErgm.bounded.idegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.idegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.idegree", arglist,
    varnames = c("d", "bound"),
    vartypes = c("numeric", "numeric"),
    defaultvalues = list(NULL, 5),
    required = c(TRUE, TRUE))
  d<-a$d
  bound<-a$bound
  if(drop){
    degrees <-
      as.numeric(names(table(table(as.edgelist(nw)[,2]))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no indegree", d[is.na(mdegrees)],
                "vertices;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.idegree", d[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.idegree() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedidegree",
                                soname="ergm",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  m$coef.names<-c(m$coef.names,paste("bounded.idegree",d,sep=""))
  m
}

#########################################################
InitErgm.bounded.istar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.istar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.istar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  k<-a$k;bound<-a$bound
  if(drop){
    mistar <- paste("c(",paste(k,collapse=","),")",sep="")
    mistar <- summary(as.formula(paste('nw ~ bounded.istar(',mistar,')',sep="")),
                        drop=FALSE) == 0
    if(any(mistar)){
      cat(" ")
      cat(paste("Warning: There are no order", k[mistar],"bounded.istars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.istar", k[mistar],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      k <- k[!mistar] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.istar() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedistar", soname="ergm",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("istar",k,".bound",bound,sep=""))
  m
}

#########################################################
InitErgm.bounded.kstar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.kstar", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("bounded.kstar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  k<-a$k;bound<-a$bound
  if(drop){
    degrees <- as.numeric(names(table(table(as.edgelist(nw)))))
    mdegrees <- match(k, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no degree", k[is.na(mdegrees)],
                "vertices;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.kstar", k[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      k <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.kstar() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedkstar", soname="ergm",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("kstar",k,".bound",bound,sep=""))
  m
}

#########################################################
InitErgm.bounded.odegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.odegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.odegree", arglist,
    varnames = c("d", "bound"),
    vartypes = c("numeric", "numeric"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, TRUE))
  d<-a$d
  bound<-a$bound
  if(drop){
    degrees <-
      as.numeric(names(table(table(as.edgelist(nw)[,1]))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(" ")
      cat(paste("Warning: There are no outdegree", d[is.na(mdegrees)],
                "vertices;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.odegree", d[is.na(mdegrees)],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.odegree() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedodegree",
                                soname="ergm",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  m$coef.names<-c(m$coef.names,paste("bounded.odegree",d,sep=""))
  m
}

#########################################################
InitErgm.bounded.ostar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.ostar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.ostar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  k<-a$k;bound<-a$bound
  if(drop){
    mostar <- paste("c(",paste(k,collapse=","),")",sep="")
    mostar <- summary(as.formula(paste('nw ~ bounded.ostar(',mostar,')',sep="")),
                        drop=FALSE) == 0
    if(any(mostar)){
      cat(" ")
      cat(paste("Warning: There are no order", k[mostar],"bounded.ostars;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
      dropterms <- paste("bounded.ostar", k[mostar],sep="")
#     cat(paste("To avoid degeneracy the terms",
#               paste(dropterms,collapse=" and, "),
#               "have been dropped.\n"))
      k <- k[!mostar] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.ostar() expects its 2 arglist to be of the ",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedostar", soname="ergm",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("ostar",k,".bound",bound,sep=""))
  m
}

#########################################################
InitErgm.bounded.triangle<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("bounded.triangle", arglist,
    varnames = c("bound"),
    vartypes = c("numeric"),
    defaultvalues = list(5),
    required = c(TRUE))
  bound<-a$bound
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedtriangle", soname="ergm",
                                inputs = c(0, 1, 1, bound))
  m$coef.names<-c(m$coef.names,paste("triangle.bound",bound,sep=""))
  m
}

###################################### InitErgm TERMS:  C

###################################### InitErgm TERMS:  D
#########################################################
InitErgm.degreep<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("degreep", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("degreep", arglist,
    varnames = c("d", "attrname", "homophily"),
    vartypes = c("numeric", "character", "logical"),
    defaultvalues = list(NULL, NULL, FALSE),
    required = c(TRUE, FALSE, FALSE))
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


###################################### InitErgm TERMS:  E

###################################### InitErgm TERMS:  F
#########################################################
InitErgm.factor<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("factor", arglist,
    varnames = c("factorname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
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
  m$terms[[termnumber]]$emptynwstats <- network.size(nw)
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
InitErgm.hamming.weighted<-function (nw, m, arglist, ...) {
# Note:  This term is now obsolete; all its functionality is in 
# InitErgm.hamming
# ergm.checkdirected("hamming.weighted", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hamming.weighted", arglist=arglist,
    varnames = c("cov","x","attrname"),
    vartypes = c("matrix,network","matrix,network","character"),
    defaultvalues = list(NULL,nw,NULL),
    required = c(TRUE,FALSE,FALSE))
  attrname<-a$attrname
  x<-a$x
  cov<-a$cov
#
# Extract dyadic covariate
#
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
    stop("hamming.weighted() requires a proper dyadic covariate", call.=FALSE)
  }
#
# Extract reference network as an edgelist
#
  if(is.network(x)){
    xm<-as.edgelist(x)
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.edgelist(xm)
  }else if(is.null(x)){
    xm<-as.edgelist(nw)
  }else{
    xm<-as.matrix(x)
    x<-paste(quote(x))
  }
  if (is.null(xm) || ncol(xm)!=2){
    stop("hamming.weighted() requires a proper network as its reference", 
         call.=FALSE)
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
  termnumber<-1+length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# component 1 is set to 1 (although in this case, input component 1
# is actually arbitrary since d_dyadcov ignores the value of inp->attrib).
   m$terms[[termnumber]] <- list(name = "hamming_weighted", soname="ergm",
                                 inputs = c(1, 1,
                                   1+2*nrow(xm)+nrow(covm)*ncol(covm),
                                   nrow(xm), as.integer(xm),
                                   as.double(covm)),
                                 dependence=TRUE)
   if(!is.null(attrname)){
     cn<-paste("hamming.weighted", as.character(sys.call(0)[[4]][2]),
               as.character(attrname), sep = ".")
   }else{
     cn<-paste("hamming.weighted", as.character(sys.call(0)[[4]][2]), sep = ".")
   }
   m$coef.names <- c(m$coef.names, cn)
   m
}

#########################################################
InitErgm.hammingmix.constant<-function (nw, m, arglist, ...) {
# ergm.checkdirected("hammingconstantmix", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("hammingmix.constant", arglist=arglist,
    varnames = c("attrname","x","base", "contrast"),
    vartypes = c("character","matrix,network","numeric","logical"),
    defaultvalues = list(NULL,nw,0,FALSE),
    required = c(TRUE,FALSE,FALSE,FALSE))
  attrname<-a$attrname
  x<-a$x
  base<-a$base
  contrast<-a$contrast
  drop<-a$drop
  drop<-TRUE
  if(is.network(x)){
    xm<-as.edgelist(x,attrname)
    x<-paste(quote(x))
  }else if(is.character(x)){
    xm<-get.network.attribute(nw,x)
    xm<-as.edgelist(xm)
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
InitErgm.heideriandynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("heideriandynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("heideriandynamic", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrix,network","character"),
    defaultvalues = list(NULL,NULL),
    required = c(TRUE,FALSE))
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
  m$terms[[termnumber]]$emptynwstats <- summary(nw~asymmetric)
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
InitErgm.intransitivedynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("intransitivedynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("intransitivedynamic", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrix,network","character"),
    defaultvalues = list(NULL,NULL),
    required = c(TRUE,FALSE))
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
   m$terms[[termnumber]]$emptynwstats <- summary(nw~intransitive)
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
  attrname <- a$attrname
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

###################################### InitErgm TERMS:  O

###################################### InitErgm TERMS:  R

###################################### InitErgm TERMS:  S
#########################################################
InitErgm.simmeliandynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("simmeliandynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("simmeliandynamic", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrix,network","character"),
    defaultvalues = list(NULL,NULL),
    required = c(TRUE,FALSE))
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
InitErgm.spatial<-function(nw, m, arglist, initialfit=FALSE, ...) {
  a <- ergm.checkargs("spatial", arglist,
    varnames = c("pb","alpha","gamma","d"),
    vartypes = c("numeric","numeric","numeric","numeric"),
    defaultvalues = list(10, -10, 0, NA),
    required = c(FALSE, FALSE, FALSE, TRUE))
  pb <- a$pb; alpha <- a$alpha; gamma <- a$gamma; d <- a$d
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
InitErgm.transitivedynamic<-function (nw, m, arglist, ...) {
  ergm.checkdirected("transitivedynamic", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("transitivedynamic", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrix,network","character"),
    defaultvalues = list(NULL,NULL),
    required = c(TRUE,FALSE))
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




