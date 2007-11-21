#  See InitErgm.R for a general explanation 
#  of InitErgm functions

###################################### InitErgm TERMS:  A
#########################################################
InitErgm.aconcurrent<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("aconcurrent", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("aconcurrent", arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  attach(a)
  attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "aconcurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to aconcurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    if(drop){ #   Check for degeneracy
      aconcurrentattr <- paste('nw ~ aconcurrent("',attrname,'")',sep="")
      aconcurrentattr <- summary(as.formula(aconcurrentattr),
                                 drop=FALSE) == 0
      if(any(aconcurrentattr)){
        cat(" ")
        cat(paste("Warning: There are no aconcurrent", ".", attrname,
                           u[aconcurrentattr],
                  "actors;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        dropterms <- paste("aconcurrent", ".", attrname,
                           u[aconcurrentattr], sep="")
        u <- u[-aconcurrentattr]
      }
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      maconcurrent <- summary(
                          as.formula(paste('nw ~ aconcurrent',sep="")),
                          drop=FALSE) == 0
      if(any(maconcurrent)){
        cat(paste("Warning: There are no concurrent actors.\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="aconcurrent_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           u, nodecov),
                                  dependence=TRUE)
    # See comment in d_aconcurrent_by_attr function
    m$coef.names<-c(m$coef.names, paste("aconcurrent",".", attrname,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="aconcurrent", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("aconcurrent",sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
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
      cat(" ")
      cat(paste("Warning: There are no actorfactor", ".", attrname,
                         u[nfc==0],
                "actors;\n",
               " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
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
  m$terms[[termnumber]] <- list(name="actorfactor", soname="ergm",
                                inputs=c(lu, lu, lu+length(nodecov),
                                         ui, nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("actorfactor",
                                      attrname, paste(u), sep="."))
  m
}

#########################################################
InitErgm.adegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("adegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("adegree", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
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
      adegreeattr <- summary(as.formula
                             (paste('nw ~ adegree(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(adegreeattr)){
        cat(" ")
        cat(paste("Warning: There are no degree", du[1,adegreeattr], ".",
                   attrname, u[du[2,adegreeattr]],
                  "actors;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        du <- matrix(du[,!adegreeattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nactors]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      madegree <- paste("c(",paste(d,collapse=","),")",sep="")
      madegree <- summary(
                          as.formula(paste('nw ~ adegree(',madegree,')',sep="")),
                          drop=FALSE) == 0
      if(any(madegree)){
        cat(" ")
        cat(paste("Warning: There are no degree", d[madegree],
                  "actors;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        d <- d[!madegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nactors
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="adegree_by_attr", soname="ergm",
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
    m$terms[[termnumber]] <- list(name="adegree", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("adegree",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.akappa<-function(nw, m, arglist, ...) {
  ergm.checkdirected("akappa", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("akappa", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="akappa", soname="ergm",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"akappa")
  m
}

#########################################################
InitErgm.astar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("astar", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("astar", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "astar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to astar() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      astarattr <- summary(as.formula
                             (paste('nw ~ astar(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(astarattr)){
        cat(" ")
        cat(paste("Warning: There are no actor stars", du[1,astarattr], ".",
                   attrname, u[du[2,astarattr]],
                  ";\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        du <- matrix(du[,!astarattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nactors]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mastar <- paste("c(",paste(d,collapse=","),")",sep="")
      mastar <- summary(
                          as.formula(paste('nw ~ astar(',mastar,')',sep="")),
                          drop=FALSE) == 0
      if(any(mastar)){
        cat(" ")
        cat(paste("Warning: There are no actor stars", d[mastar],
                  "actors;\n",
                  " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        d <- d[!mastar] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nactors
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ostar_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_astar_by_attr function
    m$coef.names<-c(m$coef.names, paste("astar", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ostar", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("astar",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

###################################### InitErgm TERMS:  E
#########################################################
InitErgm.econcurrent<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("econcurrent", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("econcurrent", arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  attach(a)
  attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "econcurrent")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to econcurrent() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    if(drop){ #   Check for degeneracy
      econcurrentattr <- paste('nw ~ econcurrent("',attrname,'")',sep="")
      econcurrentattr <- summary(as.formula(econcurrentattr),
                                 drop=FALSE) == 0
      if(any(econcurrentattr)){
        cat(" ")
        cat(paste("Warning: There are no econcurrent", ".", attrname,
                           u[econcurrentattr],
                  "events;\n",
                 " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
        u <- u[-econcurrentattr]
      }
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      meconcurrent <- summary(
                          as.formula(paste('nw ~ econcurrent',sep="")),
                          drop=FALSE) == 0
      if(any(meconcurrent)){
        cat(paste("Warning: There are no concurrent events\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="econcurrent_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           u, nodecov),
                                  dependence=TRUE)
    # See comment in d_econcurrent_by_attr function
    m$coef.names<-c(m$coef.names, paste("econcurrent",".", attrname,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="econcurrent", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("econcurrent",sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.edegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("edegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("edegree", arglist,
    varnames = c("d", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
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
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[(1+nactors):n]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      medegree <- paste("c(",paste(d,collapse=","),")",sep="")
      medegree <- summary(
       as.formula(paste('nw ~ edegree(',medegree,')',sep="")),
       drop=FALSE) == 0
      if(any(medegree)){
        cat(" ")
        cat(paste("Warning: There are no degree", d[mdegree],
                  "events;\n",
                  " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
       dropterms <- paste("edegree", d[medegree],sep="")
       d <- d[!medegree] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- n-nactors
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="edegree_by_attr", soname="ergm",
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
    m$terms[[termnumber]] <- list(name="edegree", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("edegree",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
InitErgm.ekappa<-function(nw, m, arglist, ...) {
  ergm.checkdirected("ekappa", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("ekappa", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="ekappa", soname="ergm",
                                inputs=c(0, 1, 0))
  m$coef.names<-c(m$coef.names,"ekappa")
  m
}

#########################################################
InitErgm.estar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkbipartite("estar", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("estar", arglist,
                      varnames = c("d", "attrname"),
                      vartypes = c("numeric", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  attach(a)
  d<-a$d; attrname <- a$attrname
  emptynwstats<-NULL
  nactors <- get.network.attribute(nw, "bipartite")
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "estar")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to estar() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    if(drop){ #   Check for degeneracy
      tmp <- paste("c(",paste(d,collapse=","),")")
      estarattr <- summary(as.formula
                             (paste('nw ~ estar(', tmp,',"',attrname,'")',sep="")),
                             drop=FALSE) == 0
      if(any(astarattr)){
        dropterms <- paste("estar", du[1,estarattr], ".", attrname,
                           u[du[2,estarattr]], sep="")
        cat("Warning: These estar terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
        du <- matrix(du[,!estarattr], nrow=2)
      }
    }
    if (any(du[1,]==0)) {
      emptynwstats <- rep(0, ncol(du))
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) tmp[i] <- 
        sum(nodecov[1:nactors]==tmp[i])
        emptynwstats[du[1,]==0] <- tmp
    }
  } else {
    if(is.logical(attrname)){drop <- attrname}
    if(drop){
      mestar <- paste("c(",paste(d,collapse=","),")",sep="")
      mestar <- summary(
                          as.formula(paste('nw ~ estar(',mestar,')',sep="")),
                          drop=FALSE) == 0
      if(any(mestar)){
        cat(paste("Warning: There are no order", d[mestar],"estars.\n"))
        dropterms <- paste("estar", d[mestar],sep="")
        cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
        d <- d[!mestar] 
      }
    }
    if (any(d==0)) {
      emptynwstats <- rep(0, length(d))
      emptynwstats[d==0] <- nactors
    }
  }
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="istar_by_attr", soname="ergm",
                                  inputs=c(0, ncol(du), 
                                           length(du)+length(nodecov), 
                                           as.vector(du), nodecov),
                                  dependence=TRUE)
    # See comment in d_estar_by_attr function
    m$coef.names<-c(m$coef.names, paste("estar", du[1,], ".", attrname,
                                        u[du[2,]], sep=""))
  }else{
    lengthd<-length(d)
    if(lengthd==0){return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="istar", soname="ergm",
                                       inputs=c(0, lengthd, lengthd, d),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("estar",d,sep=""))
  }
  if (!is.null(emptynwstats)) 
    m$terms[[termnumber]]$emptynwstats <- emptynwstats
  m
}

#########################################################
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
      cat(" ")
      cat(paste("Warning: There are no eventfactor", ".", attrname,
                         u[nfc==0],
                "events;\n",
               " the corresponding coefficient has been fixed at its MLE of nenegative infinity.\n",sep=" "))
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
  m$terms[[termnumber]] <- list(name="eventfactor", soname="ergm",
                                inputs=c(lu, lu, lu+length(nodecov),
                                         ui, nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("eventfactor",
                                      attrname, paste(u), sep="."))
  m
}

###################################### InitErgm TERMS:  G
#########################################################
InitErgm.gwadegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkbipartite("gwadegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwadegree", arglist,
    varnames = c("decay", "fixed", "attrname"),
    vartypes = c("numeric", "logical", "character"),
    defaultvalues = list(0, FALSE, NULL),
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
    m$terms[[termnumber]] <- list(name="adegree", soname="ergm",
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
    m$terms[[termnumber]] <- list(name="gwadegree_by_attr", soname="ergm",
                                  inputs=c(0, lu, 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    # See comment in d_gwadegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("gwadeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwadegree", soname="ergm",
                                       inputs=c(0, 1, 1, decay),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("gwadeg",decay,sep=""))
  }
  m
}

#########################################################
InitErgm.gwedegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkbipartite("gwedegree", is.bipartite(nw), requirement=TRUE)
  a <- ergm.checkargs("gwedegree", arglist,
    varnames = c("decay", "fixed", "attrname"),
    vartypes = c("numeric", "logical", "character"),
    defaultvalues = list(0, FALSE, NULL),
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
    m$terms[[termnumber]] <- list(name="adegree", soname="ergm",
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
    m$terms[[termnumber]] <- list(name="gwedegree_by_attr", soname="ergm",
                                  inputs=c(0, lu,
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    # See comment in d_gwedegree_by_attr function
    m$coef.names<-c(m$coef.names, paste("gwedeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwedegree", soname="ergm",
                                       inputs=c(0, 1, 1, decay),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("gwedeg",decay,sep=""))
  }
  m
}


###################################### InitErgm TERMS:  H
#########################################################
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
    stop("hammingdyadcov() requires a proper dyadic covariate")
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
    stop("hammingdyadcov() requires a proper network as its reference")
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
   m$terms[[termnumber]] <- list(name = "hammingdyadcov", soname="ergm",
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

#########################################################
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
  m$terms[[termnumber]] <- list(name="hammingfixmix", soname="ergm",
    inputs=c(1, 1, nrow(xm)*2+length(nodecov)+1,
            nrow(xm),as.integer(xm), nodecov),
            dependence=FALSE)
  m$coef.names<-c(m$coef.names, paste("hammingfixmix",attrname, sep="."))
  m
}

#########################################################
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
  if(all(omitwhich!=0)){
   u <- u[-omitwhich,]
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




