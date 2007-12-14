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
# Required arguments of output list:
# -----------------------------------
#       name: This is the (text) name of the term.  It is expected that there
#             is a C function called d_[name].
# changestatnames: Names of the changestats.  This is a vector of character
#                  strings whose length must be equal to the number of change
#                  statistics calculated by this ERGM term.
# inputs: vector of numeric inputs to the ERGM term.  This can include parameters
#         as well as nodal or dyadic covariates.

#
# Optional arguments of output list:
# ----------------------------------
#        soname: This is the (text) name of the package containing the C function
#                called d_[name].
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


#################################################################
# InitERGMterm functions (new version:  June 2007)
# The old InitErgm functions should still work.
#
# INPUT:
# Each InitERGMterm function takes two arguments, nw and arglist,
# which are automatically supplied by ergm.getmodel.  There may be
# other arguments passed by ergm.getmodel, so each InitERGMterm
# function should also include the ... argument in its list.
#
# OUTPUT:
# Each InitERGMterm function should return a list.  
#    REQUIRED LIST ITEMS:
#          name: Name of the C function that produces the change
#                statistics.  (Note:  The name will have "d_" 
#                prepended.  For example, the C function for the
#                absdiff change statistics is called "d_absdiff"
#                even though InitERGMterms.absdiff only returns
#                names = "absdiff")
#    coef.names: Vector of names for the coefficients (parameters)
#                as they will be reported in the output.
#
#    OPTIONAL LIST ITEMS:
#        inputs: Vector of inputs (of type double) that the
#                d_xxx function will require.  Default is NULL.
#        soname: This is the (text) name of the package containing the C function
#                called d_[name].  Default is "statnet"
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

InitERGMterm.absdiff <- function(nw, arglist, ...) {
  # First, make sure this ERGM term is right for this network (optional!)
  check.InitERGM.situation(nw, directed=NULL, bipartite=NULL)
  # Second, check the arguments to make sure they are appropriate.
  a <- get.InitERGMterm.args(arglist,
                             varnames = c("attrname"),
                             vartypes = c("character"),
                             defaultvalues = list(NULL),
                             required = c(TRUE))  
  # Process the arguments
  nodecov <- get.node.attr(nw, a$attrname, "absdiff")
  # Construct the output list
  list(name="absdiff",                                    #name: required
       coef.names = paste("absdiff", a$attrname, sep=""), #coef.names: required
       inputs = nodecov,
       soname = "statnet",
       dependence = FALSE
       )
}
     


#InitErgm.absdiffcat<-function (nw, m, arglist, ...) {
#  a <- ergm.checkargs("absdiffcat", arglist,
#                   varnames = c("attrname","omitwhich"),
#                   vartypes = c("character","numeric"),
#                   defaultvalues = list(NULL,NULL),
#                   required = c(TRUE,FALSE))
#  # omitwhich:  If not NULL, indices of non-zero absdiff categories to delete
#  # (ordered from 1=smallest to largest).
#  attach(a)
#  attrname<-a$attrname
#  omitwhich <- a$omitwhich
#  nodecov <- get.node.attr(nw, attrname, "absdiffcat")  
#  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
#  u <- u[u>0]
#  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=T)) # Arbitrary unused (and nonzero) value
#  napositions <- is.na(nodecov)
#  nodecov[napositions] <- NAsubstitute
#  if(any(napositions)){u<-c(u,NA)}
#  if(!is.null(omitwhich)) u <- u[-omitwhich]
#  if (length(u)==0)
#    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
#  termnumber<-1+length(m$terms)  
#  u2 <- u[!is.na(u)]
#  m$terms[[termnumber]] <- list(name="absdiffcat", soname="statnet",
#                                inputs=c(length(u2)+1, length(u),
#                                         length(u2)+1+length(nodecov),
#                                         u2, NAsubstitute, nodecov),
#                                dependence=FALSE)
#  m$coef.names<-c(m$coef.names,paste("absdiff", attrname, u, sep="."))
#  m
#}
#
#InitErgm.bounded.degree<-function(nw, m, arglist, drop=TRUE, ...) {
#  ergm.checkdirected("bounded.degree", is.directed(nw), requirement=FALSE)
#  a <- ergm.checkargs("bounded.degree", arglist,
#    varnames = c("d","bound"),
#    vartypes = c("numeric","numeric"),
#    defaultvalues = list(NULL,5),
#    required = c(TRUE,TRUE))
#  attach(a)
#  if(drop){
#    degrees <- as.numeric(names(table(table(as.matrix.network.edgelist(nw)))))
#    degrees[degrees > max(d)] <- max(d)
#    mdegrees <- match(d, degrees)  
#    if(any(is.na(mdegrees))){
#      cat(paste("Warning: There are no degree", d[is.na(mdegrees)],
#                "vertices.\n"))
#      dropterms <- paste("bounded.degree", d[is.na(mdegrees)],sep="")
#      cat(paste("To avoid degeneracy the terms",dropterms,
#                "have been dropped.\n"))
#      d <- degrees[mdegrees[!is.na(mdegrees)]] 
#    }
#  }
#  m$coef.names<-c(m$coef.names,paste("bounded.degree",d,sep=""))
#  ld<-length(d)
#  if(ld==0){return(m)}
#  if (length(bound)!=ld)
#    stop(paste("bounded.degree() expects its 2 arglist to be of the",
#               "same length"), call.=FALSE)
#  termnumber<-1+length(m$terms)
#  m$terms[[termnumber]] <- list(name="boundeddegree", soname="statnet",
#                                inputs = c(0, ld, ld+ld, c(d,bound)))
#  m
#}


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
      cat(paste("Warning: There are no order", k[mcycle],"cycles.\n"))
      dropterms <- paste("cycle", k[mcycle],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
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
  m$terms[[termnumber]] <- list(name="cycle", soname="statnet",
                                inputs=c(0, lk, mk+1, direct, mk, usestats))
  m$coef.names<-c(m$coef.names,paste("cycle",k,sep=""))
  m
}

InitErgm.density<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("density", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="density", soname="statnet",
                                inputs=c(0, 1, 0),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,"density")
  m
}

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
        cat("Warning: These degree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
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
        cat("Warning: These degree terms have extreme counts and will be dropped:\n")
        cat(d[mdegree], "\n", fill=T)
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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("degree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="degree_w_homophily", soname="statnet",
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
    m$terms[[termnumber]] <- list(name="degree_by_attr", soname="statnet",
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




