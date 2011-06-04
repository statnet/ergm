# These are InitErgm functions that were never included in the public release.

#########################################################
InitErgmTerm.cyclicalties<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
 ### Process the arguments
  ### Construct the list to return
  out <- list(name="cyclicalties",                      #name: required
              coef.names = "cyclicalties"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "cyclicalties")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties() has only one value", call.=FALSE)
    if(drop){ # Check for zero statistics, print -Inf messages if applicable
      triattr <- summary(as.formula(paste('nw ~ cyclicalties(','"',a$attrname,
                                          '",diff=',a$diff,')',sep="")),
                         drop=FALSE) == 0
      if(a$diff){
        if(any(triattr)){
          dropterms <- paste(paste("cyclicalties",a$attrname,sep="."),
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
          dropterms <- paste(paste("cyclicalties",a$attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!a$diff) {
      out$coef.names <- paste("cyclicalties", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("cyclicalties", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "cyclicalties"
    out$inputs <- c(0,1,0)
  }
  out
}

#########################################################
InitErgmTerm.cyclicalties2<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
 ### Process the arguments
  ### Construct the list to return
  out <- list(name="cyclicalties2",                      #name: required
              coef.names = "cyclicalties2"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "cyclicalties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to cyclicalties2() has only one value", call.=FALSE)
    if(drop){ # Check for zero statistics, print -Inf messages if applicable
      triattr <- summary(as.formula(paste('nw ~ cyclicalties2(','"',a$attrname,
                                          '",diff=',a$diff,')',sep="")),
                         drop=FALSE) == 0
      if(a$diff){
        if(any(triattr)){
          dropterms <- paste(paste("cyclicalties2",a$attrname,sep="."),
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
          dropterms <- paste(paste("cyclicalties2",a$attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!a$diff) {
      out$coef.names <- paste("cyclicalties2", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("cyclicalties2", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "cyclicalties2"
    out$inputs <- c(0,1,0)
  }
  out
}

InitErgmTerm.concurrentties<-function(nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
                      varnames = c("byarg"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  if(!is.null(a$byarg)) {
    nodecov <- get.node.attr(nw, a$byarg, "concurrentties")
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
                             (paste('nw ~ concurrentties(',a$byarg,'")',sep="")),
                             drop=FALSE) == 0
      if(any(concurrentattr)){
        dropterms <- paste("concurrentties", ".", a$byarg,
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
    if(is.logical(a$byarg)){drop <- a$byarg}
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
  out <- list(name="concurrentties",                      #name: required
              coef.names = "concurrentties"               #coef.names: required
              ) 
  if(!is.null(a$byarg)) {
    #  No covariates here, so input component 1 is arbitrary
    out$coef.names <- paste("concurrentties",".", a$byarg, u, sep="")
    out$inputs <- c(0, length(u), length(u)+length(nodecov), ui, nodecov)
    out$name="concurrent_ties_by_attr"
    # See comment in d_concurrent_ties_by_attr function
  }else{
    #  No covariates here, so input component 1 is arbitrary
    out$inputs <- c(0, 1, 0)
    out$coef.names <- "concurrentties"
  }
  out
}

#########################################################
InitErgmTerm.transitiveties2<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  out <- list(name="transitiveties2",                      #name: required
              coef.names = "transitiveties2"               #coef.names: required
              ) 
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "transitiveties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to transitiveties2() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ transitiveties2(','"',a$attrname,
                                          '",diff=',a$diff,')',sep="")),
                         drop=FALSE) == 0
      if(a$diff){
        if(any(triattr)){
          dropterms <- paste(paste("transitiveties2",a$attrname,sep="."),
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
          dropterms <- paste(paste("transitiveties2",a$attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!a$diff) {
      out$coef.names <- paste("transitiveties2", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("transitiveties2", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "transitiveties2"
    out$inputs <- c(0,1,0)
  }
  out
}

InitErgmTerm.concurrentties2<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("concurrentties2", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("concurrentties2", arglist,
                      varnames = c("byarg"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(FALSE))
  out <- list(name="cyclicalties",                      #name: required
              coef.names = "cyclicalties"               #coef.names: required
              ) 
  if(!is.null(a$byarg)) {
    nodecov <- get.node.attr(nw, a$byarg, "concurrentties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    if (length(u)==1)
      stop ("Attribute given to concurrentties2() has only one value", call.=FALSE)
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    ui <- seq(along=u)
    if(drop){ #   Check for degeneracy
      concurrentattr <- summary(as.formula
                             (paste('nw ~ concurrentties2(',a$byarg,'")',sep="")),
                             drop=FALSE) == 0
      if(any(concurrentattr)){
        dropterms <- paste("concurrentties2", ".", a$byarg,
                           u[concurrentattr], sep="")
      cat(" ")
        cat("Warning: These concurrentties2 terms have extreme counts and will be dropped:\n")
        cat(dropterms, "", fill=TRUE)
        cat("  The corresponding coefficients have been fixed at their MLE of negative infinity.\n")
        u <- u[-concurrentattr]
        ui <- ui[-concurrentattr]
      }
    }
  } else {
    if(is.logical(a$byarg)){drop <- a$byarg}
    if(drop){
      mconcurrent <- summary(
                          as.formula(paste('nw ~ concurrentties2',sep="")),
                          drop=FALSE) == 0
      if(any(mconcurrent)){
      cat(" ")
        cat(paste("Warning: There are no concurrentties2 b1s;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
        return(m)
      }
    }                                                         
  }
  out <- list(name="concurrentties2",                      #name: required
              coef.names = "concurrentties2"               #coef.names: required
              ) 
  if(!is.null(a$byarg)) {
    #  No covariates here, so input component 1 is arbitrary
    out$coef.names <- paste("concurrentties2",".", a$byarg, u, sep="")
    out$inputs <- c(0, length(u), length(u)+length(nodecov), ui, nodecov)
    out$name="concurrent_ties2_by_attr"
    # See comment in d_concurrent_ties_by_attr function
  }else{
    #  No covariates here, so input component 1 is arbitrary
    out$inputs <- c(0, 1, 0)
    out$coef.names <- "concurrentties2"
  }
  out
}

#########################################################
InitErgmTerm.transitiveties2<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,bipartite=NULL,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attrname <- a$attrname
  if(!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname, "transitiveties2")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to transitiveties2() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ transitiveties2(','"',a$attrname,
                                          '",diff=',a$diff,')',sep="")),
                         drop=FALSE) == 0
      if(a$diff){
        if(any(triattr)){
          dropterms <- paste(paste("transitiveties2",a$attrname,sep="."),
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
          dropterms <- paste(paste("transitiveties2",a$attrname,sep="."),sep="")
      cat(" ")
          cat(paste("Warning: The count of", dropterms, "is extreme;\n",
                 " the corresponding coefficient has been fixed at its MLE of negative infinity.\n",sep=" "))
#         cat(paste("To avoid degeneracy the term",
#               paste(dropterms,collapse=" and, "),
#                   "have been dropped.\n"))
        }
      }
    }
    if (!a$diff) {
      out$coef.names <- paste("transitiveties2", a$attrname, sep=".")
      out$inputs <- c(0,1,length(nodecov),nodecov)
    } else {
      out$coef.names <- paste("transitiveties2", a$attrname, u, sep=".")
      out$inputs <- c(length(ui), length(ui),
                      length(ui)+length(nodecov),
                      ui, nodecov)
    }
  }else{
    out$coef.names <- "transitiveties2"
    out$inputs <- c(0,1,0)
  }
  out
}
