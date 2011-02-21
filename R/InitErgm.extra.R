# These are InitErgm functions that were never included in the public release.

InitErgm.concurrentties<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("concurrentties", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("concurrentties", arglist,
                      varnames = c("byarg"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  byarg <- a$byarg
  if(!is.null(byarg)) {
    nodecov <- get.node.attr(nw, byarg, "concurrentties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrentties() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
    if(drop){ #   Check for degeneracy
      concurrentattr <- summary(as.formula
                             (paste('nw ~ concurrentties(',byarg,'")',sep="")),
                             drop=FALSE) == 0
      if(any(concurrentattr)){
        dropterms <- paste("concurrentties", ".", byarg,
                           u[concurrentattr], sep="")
      cat(" ")
        cat("Warning: These concurrentties terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        u <- u[-concurrentattr]
        ui <- ui[-concurrentattr]
      }
    }
  } else {
    if(is.logical(byarg)){drop <- byarg}
    if(drop){
      mconcurrent <- summary(
                          as.formula(paste('nw ~ concurrentties',sep="")),
                          drop=FALSE) == 0
      if(any(mconcurrent)){
      cat(" ")
        cat(paste("Warning: There are no concurrentties b1s;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        return(m)
      }
    }                                                         
  }
  termnumber<-1+length(m$terms)
  if(!is.null(byarg)) {
    if(length(u)==0) {return(m)}
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="concurrent_ties_by_attr", soname="ergm",
                                  inputs=c(0, length(u), 
                                           length(u)+length(nodecov), 
                                           ui, nodecov),
                                  dependence=TRUE)
    # See comment in d_concurrent_ties_by_attr function
    m$coef.names<-c(m$coef.names, paste("concurrentties",".", byarg,
                                        u, sep=""))
  }else{
    #  No covariates here, so input component 1 is arbitrary
    m$terms[[termnumber]] <- list(name="concurrent_ties", soname="ergm",
                                       inputs=c(0, 1, 0),
                                       dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("concurrentties",sep=""))
  }
  m
}

#########################################################
InitErgm.transitiveties2<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("transitiveties2", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("transitiveties2", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "transitiveties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to transitiveties2() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ transitiveties2(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("transitiveties2",attrname,sep="."),
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
          dropterms <- paste(paste("transitiveties2",attrname,sep="."),sep="")
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
      m$terms[[termnumber]] <- list(name="transitiveties2", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("transitiveties2",attrname,sep="."))
     } else {
       m$terms[[termnumber]] <- list(name="transitiveties2", soname="ergm",
                                     inputs=c(length(ui), length(ui),
                                       length(ui)+length(nodecov),
                                       ui, nodecov))
       m$coef.names<-c(m$coef.names,paste("transitiveties2",
                                          attrname, u, sep="."))
     }
  }else{
    m$terms[[termnumber]] <- list(name="transitiveties2", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"transitiveties2")
  }
  m
}

#########################################################
InitErgm.cyclicalties2<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("cyclicalties2", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("cyclicalties2", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "cyclicalties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties2() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ cyclicalties2(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("cyclicalties2",attrname,sep="."),
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
          dropterms <- paste(paste("cyclicalties2",attrname,sep="."),sep="")
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
      m$terms[[termnumber]] <- list(name="cyclicalties2", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("cyclicalties2",attrname,sep="."))
     } else {
       m$terms[[termnumber]] <- list(name="cyclicalties2", soname="ergm",
                                     inputs=c(length(ui), length(ui),
                                       length(ui)+length(nodecov),
                                       ui, nodecov))
       m$coef.names<-c(m$coef.names,paste("cyclicalties2",
                                          attrname, u, sep="."))
     }
  }else{
    m$terms[[termnumber]] <- list(name="cyclicalties2", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"cyclicalties2")
  }
  m
}

#########################################################
InitErgm.cyclicalties<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("cyclicalties", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("cyclicalties", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attrname <- a$attrname
  diff <- a$diff
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "cyclicalties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ cyclicalties(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("cyclicalties",attrname,sep="."),
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
          dropterms <- paste(paste("cyclicalties",attrname,sep="."),sep="")
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
      m$terms[[termnumber]] <- list(name="cyclicalties", soname="ergm",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("cyclicalties",attrname,sep="."))
     } else {
       m$terms[[termnumber]] <- list(name="cyclicalties", soname="ergm",
                                     inputs=c(length(ui), length(ui),
                                       length(ui)+length(nodecov),
                                       ui, nodecov))
       m$coef.names<-c(m$coef.names,paste("cyclicalties",
                                          attrname, u, sep="."))
     }
  }else{
    m$terms[[termnumber]] <- list(name="cyclicalties", soname="ergm",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"cyclicalties")
  }
  m
}

