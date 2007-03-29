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

ergm.checkargs <- function(fname, arglist, varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL) {
# ergm.checkargs will check to make sure that the inputs to each model
# term are of the correct type, assign default values if applicable,
# and then return a list with elements (var1=var1value, var2=var2value, etc.)
# Note that the list returned will contain the maximum possible number of
# arguments; any arguments without values are returned as NULL.
#   Inputs:  fname is a character giving the name of the model term.
#            arglist is the list of arguments passed from ergm.getmodel
#            varnames is a vector of variable names
#            vartypes is a vector of corresponding variable types
#            defaultvalues is a list of default values (NULL means no default)
#            required is a vector of logicals:  Is this var required or not?
  sr=sum(required)
  lv=length(varnames)
  la=length(arglist)
  if(la < sr || la > lv) {
    if (sr < lv)
      expected = paste("from",sr,"to",lv,"arguments,")
    else if(sr==1)
      expected = "1 argument,"
    else
      expected = paste(sr,"arguments,")
    stop(paste(fname,"model term expected", expected, "got", la), call.=FALSE)
  }
# The correctness of what the user typed is checked, but it is assumed
# that each InitErgm function faithfully passes in what the user typed;
# thus, the correctness of input from the InitErgm function isn't checked.
  out = defaultvalues
  names(out)=varnames
  m=NULL
  if (la>0) {
    for(i in 1:la) { # check each arglist entry
      if (!is.null(names(arglist)) && (name <- names(arglist)[i]) != "") {
        m = pmatch(name, varnames)# try to match user-typed name if applicable
        if(is.na(m)) { # User typed an unrecognizable name
          stop(paste(fname,"model term does not recognize",
                     name, "argument"), call.=FALSE)
        }
        # valid name match with mth variable if we got to here
        if (!eval(call(paste("is.",vartypes[m],sep=""),arglist[[i]]))) {
          # Wrong type
          stop(paste(name, "argument to", fname, "model term is not of",
                     "the expected", vartypes[m], "type"), call.=FALSE)
        }
        # correct type if we got to here
        out[[m]]=arglist[[i]]
      } else { # no user-typed name for this argument
        if (!is.null(m)) {
          stop(paste("unnamed argument follows named argument in",
                     fname,"model term"), call.=FALSE)
        }
        if (!eval(call(paste("is.",vartypes[i],sep=""),arglist[[i]]))) {
          # Wrong type
          stop(paste("argument number", i, "to", fname, "model term is not",
                     "of the expected", vartypes[i], "type"), call.=FALSE)
        }
        # correct type if we got to here
        out[[i]]=arglist[[i]]
      }
    }
  }
  c(.conflicts.OK=TRUE,out)
}

ergm.checkdirected <- function(fname, nw.directedflag, requirement,
                               extramessage="") {
  if (!nw.directedflag && requirement)
    stop(paste(fname, "model term may not be used with an undirected network.",
               extramessage), call.=FALSE)
  if (nw.directedflag && !requirement)
    stop(paste(fname, "model term may not be used with a directed network.",
               extramessage), call.=FALSE)
}

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
  nodecov <- get.node.attr(nw, attrname, "absdiff")
  m$terms[[termnumber]] <- list(name="absdiff", soname="statnet",
                                inputs=c(0,1,length(nodecov),nodecov),
                                dependence=FALSE)
  m
}

InitErgm.absdiffcat<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("absdiffcat", arglist,
                   varnames = c("attrname","omitwhich"),
                   vartypes = c("character","numeric"),
                   defaultvalues = list(NULL,NULL),
                   required = c(TRUE,FALSE))
  # omitwhich:  If not NULL, indices of non-zero absdiff categories to delete
  # (ordered from 1=smallest to largest).
  attach(a)
  attrname<-a$attrname
  omitwhich <- a$omitwhich
  nodecov <- get.node.attr(nw, attrname, "absdiffcat")  
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=T)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(omitwhich)) u <- u[-omitwhich]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  termnumber<-1+length(m$terms)  
  u2 <- u[!is.na(u)]
  m$terms[[termnumber]] <- list(name="absdiffcat", soname="statnet",
                                inputs=c(length(u2)+1, length(u),
                                         length(u2)+1+length(nodecov),
                                         u2, NAsubstitute, nodecov),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,paste("absdiff", attrname, u, sep="."))
  m
}

InitErgm.bounded.degree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.degree", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("bounded.degree", arglist,
    varnames = c("d","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  attach(a)
  if(drop){
    degrees <- as.numeric(names(table(table(as.matrix.network.edgelist(nw)))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(paste("Warning: There are no degree", d[is.na(mdegrees)],
                "vertices.\n"))
      dropterms <- paste("bounded.degree", d[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  m$coef.names<-c(m$coef.names,paste("bounded.degree",d,sep=""))
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.degree() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundeddegree", soname="statnet",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  m
}

InitErgm.bounded.idegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.idegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.idegree", arglist,
    varnames = c("d", "bound"),
    vartypes = c("numeric", "numeric"),
    defaultvalues = list(NULL, 5),
    required = c(TRUE, TRUE))
  attach(a)
  d<-a$d
  if(drop){
    degrees <-
      as.numeric(names(table(table(as.matrix.network.edgelist(nw)[,2]))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(paste("Warning: There are no indegree", d[is.na(mdegrees)],
                "vertices.\n"))
      dropterms <- paste("bounded.idegree", d[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.idegree() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedidegree",
                                soname="statnet",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  m$coef.names<-c(m$coef.names,paste("bounded.idegree",d,sep=""))
  m
}

InitErgm.bounded.istar<-function(nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("bounded.istar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  attach(a)
  k<-a$k;bound<-a$bound
  if(drop){
    mistar <- paste("c(",paste(k,collapse=","),")",sep="")
    mistar <- summary(as.formula(paste('nw ~ bounded.istar(',mistar,')',sep="")),
                        drop=FALSE) == 0
    if(any(mistar)){
      cat(paste("Warning: There are no order", k[mistar],"bounded.istars.\n"))
      dropterms <- paste("bounded.istar", k[mistar],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      k <- k[!mistar] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.istar() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedistar", soname="statnet",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("istar",k,".bound",bound,sep=""))
  m
}

InitErgm.bounded.kstar<-function(nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("bounded.kstar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  attach(a)
  k<-a$k;bound<-a$bound
  if(drop){
    degrees <- as.numeric(names(table(table(as.matrix.network.edgelist(nw)))))
    mdegrees <- match(k, degrees)  
    if(any(is.na(mdegrees))){
      cat(paste("Warning: There are no degree", k[is.na(mdegrees)],
                "vertices.\n"))
      dropterms <- paste("bounded.kstar", k[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      k <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.kstar() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedkstar", soname="statnet",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("kstar",k,".bound",bound,sep=""))
  m
}

InitErgm.bounded.odegree<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("bounded.odegree", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("bounded.odegree", arglist,
    varnames = c("d", "bound"),
    vartypes = c("numeric", "numeric"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, TRUE))
  attach(a)
  d<-a$d
  if(drop){
    degrees <-
      as.numeric(names(table(table(as.matrix.network.edgelist(nw)[,1]))))
    degrees[degrees > max(d)] <- max(d)
    mdegrees <- match(d, degrees)  
    if(any(is.na(mdegrees))){
      cat(paste("Warning: There are no outdegree", d[is.na(mdegrees)],
                "vertices.\n"))
      dropterms <- paste("bounded.odegree", d[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  if (length(bound)!=ld)
    stop(paste("bounded.odegree() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedodegree",
                                soname="statnet",
                                inputs = c(0, ld, ld+ld, c(d,bound)))
  m$coef.names<-c(m$coef.names,paste("bounded.odegree",d,sep=""))
  m
}

InitErgm.bounded.ostar<-function(nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("bounded.ostar", arglist,
    varnames = c("k","bound"),
    vartypes = c("numeric","numeric"),
    defaultvalues = list(NULL,5),
    required = c(TRUE,TRUE))
  attach(a)
  k<-a$k;bound<-a$bound
  if(drop){
    mostar <- paste("c(",paste(k,collapse=","),")",sep="")
    mostar <- summary(as.formula(paste('nw ~ bounded.ostar(',mostar,')',sep="")),
                        drop=FALSE) == 0
    if(any(mostar)){
      cat(paste("Warning: There are no order", k[mostar],"bounded.ostars.\n"))
      dropterms <- paste("bounded.ostar", k[mostar],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      k <- k[!mostar] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  if (length(bound)!=lk)
    stop(paste("bounded.ostar() expects its 2 arglist to be of the",
               "same length"), call.=FALSE)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedostar", soname="statnet",
                                inputs = c(0, lk, lk+lk, c(k,bound)))
  m$coef.names<-c(m$coef.names,paste("ostar",k,".bound",bound,sep=""))
  m
}

InitErgm.bounded.triangle<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("bounded.triangle", arglist,
    varnames = c("bound"),
    vartypes = c("numeric"),
    defaultvalues = list(5),
    required = c(TRUE))
  attach(a)
  bound<-a$bound
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="boundedtriangle", soname="statnet",
                                inputs = c(0, 1, 1, bound))
  m$coef.names<-c(m$coef.names,paste("triangle.bound",bound,sep=""))
  m
}

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
      cat(paste("Warning: There are no",attrname," ties for the vertex", 
                d[centattr],"\n"))
      dropterms <- paste("sociality", d[centattr],sep="")
      cat(paste("To avoid degeneracy the term",dropterms,
                "have been dropped.\n"))
      d <- d[!centattr] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="sociality", soname="statnet",
                                  inputs=c(0, ld, ld+length(nodecov),
                                    d, nodecov))
    m$coef.names<-c(m$coef.names,paste("sociality",d,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="sociality", soname="statnet",
                                  inputs=c(0, ld, ld, d))
    m$coef.names<-c(m$coef.names,paste("sociality",d,sep=""))
  }
  m
}

InitErgm.ctriad<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ctriad", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ctriad", arglist,
    varnames = c("attrname","diff"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,FALSE),
    required = c(FALSE,FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    nodecov <- get.node.attr(nw, attrname, "ctriad")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ctriad() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ ctriad(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("ctriad",attrname,sep="."),
                             u[triattr],sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the terms",dropterms,
                    "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("ctriad",attrname,sep="."),sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the term",dropterms,
                    "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="ctriad", soname="statnet",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("ctriad",attrname,sep="."))
    } else {
      #  Number of input parameters before covariates equals number of
      #  unique elements in nodecov, namely length(u), so that's what
      #  input element 1 equals
      m$terms[[termnumber]] <- list(name="ctriad", soname="statnet",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("ctriad", attrname, u, sep="."))
    }
  }else{
#    No attributes (or diff)
#    No covariates, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="ctriad", soname="statnet",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"ctriad")
  }
  m
}

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

InitErgm.dsp<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("dsp", is.directed(nw), requirement=FALSE)
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
      cat(paste("Warning: There are no dsp", d[mdsp==0],"dyads.\n"))
      dropterms <- paste("dsp", d[mdsp==0],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      d <- d[mdsp!=0] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="dsp", soname="statnet",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("dsp",d,sep=""))
  m
}

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
  m$terms[[termnumber]] <- list(name = "duration", soname="statnet",
                                inputs = c(1, 1, 
                                  NROW(xm)*2+2*NROW(formm)^2, NROW(xm),
                                  as.double(c(xm, formm, dissolvem))))
  m
}

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
   m$terms[[termnumber]] <- list(name = "dyadcov",  soname="statnet",
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
   m$terms[[termnumber]] <- list(name = "dyadcov", soname="statnet", 
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
   m$terms[[termnumber]] <- list(name = "simmeliandynamic", soname="statnet", 
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
   m$terms[[termnumber]] <- list(name = "intransitivedynamic", soname="statnet", 
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
   m$terms[[termnumber]] <- list(name = "transitivedynamic", soname="statnet", 
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
  m$terms[[termnumber]] <- list(name = "heideriandynamic", soname="statnet", 
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
    mesp <- summary(nw ~ factor(x), drop=FALSE)
    if(all(mesp==0)){return(m)}
    if(any(mesp==0)){
      cat(paste("Warning: There are no dyads with factor level",
          names(mesp)[mesp==0],".\n"))
      cat(paste("To avoid degeneracy the terms",names(mesp)[mesp==0],
                "have been dropped.\n"))
      xm <- xm[,mesp!=0] 
    }
  }

  #Update the term number
  termnumber <- 1 + length(m$terms)
  #Update the terms list, adding the vectorized adjacency matrix

# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1 (although in this case, input element 1
# is actually arbitrary since d_factor ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "factor",  soname="statnet",
                                inputs = c(1, 3, 1+NROW(xm)*NCOL(xm),
                                  NCOL(xm), as.double(as.numeric(xm))),
                                dependence=FALSE)
  cn <- paste(factorname,substring(colnames(xm),first=2),sep=".")
  m$coef.names<-c(m$coef.names, cn)
  m
}

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
  m$terms[[termnumber]] <- list(name = "edgecov", soname="statnet", 
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

InitErgm.edges<-function(nw, m, arglist, ...) {
  a <- ergm.checkargs("edges", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="edges", soname="statnet",
                                inputs=c(0, 1, 0),
                                dependence=FALSE)
  m$coef.names<-c(m$coef.names,"edges")
  m
}

InitErgm.esp<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("esp", is.directed(nw), requirement=FALSE)
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
      cat(paste("Warning: There are no dyads with esp", d[mesp==0],".\n"))
      dropterms <- paste("esp", d[mesp==0],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      d <- d[mesp!=0] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="esp", soname="statnet",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("esp",d,sep=""))
  m
}

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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegreealpha=NULL,
                                    gwdegree.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegreealpha#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegreealpha", soname="statnet",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwdegreealpha")
  }
  m
}

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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegreelambda=NULL,
                                    gwdegree.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegreelambda#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegreelambda", soname="statnet",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"gwdegree")
  }
  m
}

InitErgm.gwdegree<-function(nw, m, arglist, initialfit=FALSE, ...) {
# ergm.checkdirected("gwdegree", is.directed(nw), requirement=FALSE)
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
      stop("The gwadegree term is not yet able to handle a",
           "nonfixed decay term with an attribute.")
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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
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
    m$terms[[termnumber]] <- list(name="gwdegree_by_attr", soname="statnet",
                                  inputs=c(0, lu, 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("gwdeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwdegree", soname="statnet",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwdegree")
  }
  m
}

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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdegree=NULL,
                                    gwdegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdegree#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="gwdegree706", soname="statnet",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwdegree706")
  }
  m
}

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
      stop("The gwidegree term is not yet able to handle a",
           "nonfixed decay term with an attribute.")
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
    m$terms[[termnumber]] <- list(name="idegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwidegree=NULL,
                                    gwidegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwidegree#",d,sep=""))
  } else if(!is.null(attrname)) {
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
    m$terms[[termnumber]] <- list(name="gwidegree_by_attr", soname="statnet",
                                  inputs=c(0, lu, 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("gwideg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwidegree", soname="statnet",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwidegree")
  }
  m
}

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
      stop("The gwodegree term is not yet able to handle a",
           "nonfixed decay term with an attribute.")
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
    m$terms[[termnumber]] <- list(name="odegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwodegree=NULL,
                                    gwodegree.decay=decay),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwodegree#",d,sep=""))
  } else if(!is.null(attrname)) {
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
    m$terms[[termnumber]] <- list(name="gwodegree_by_attr", soname="statnet",
                                  inputs=c(0, lu, 
                                           1+length(nodecov), 
                                           decay, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names, paste("gwodeg", decay, ".", 
                                        attrname, u, sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="gwodegree", soname="statnet",
                                  inputs=c(0, 1, length(decay), decay))
    m$coef.names<-c(m$coef.names,"gwodegree")
  }
  m
}

InitErgm.altkstar<-function(nw, m, arglist, initialfit=FALSE, ...) {
# ergm.checkdirected("altkstar", is.directed(nw), requirement=FALSE)
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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altkstar=NULL,
                                    altkstar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altkstar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altkstar", soname="statnet",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altkstar")
  }
  m
}

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
    m$terms[[termnumber]] <- list(name="idegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altistar=NULL,
                                    altistar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altistar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altistar", soname="statnet",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altistar")
  }
  m
}

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
    m$terms[[termnumber]] <- list(name="odegree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(altostar=NULL,
                                    altostar.lambda=lambda),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("altostar#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="altostar", soname="statnet",
                                  inputs=c(0, 1, length(lambda), lambda))
    m$coef.names<-c(m$coef.names,"altostar")
  }
  m
}

InitErgm.gwdsp<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwdsp", is.directed(nw), requirement=FALSE)
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
    m$terms[[termnumber]] <- list(name="dsp", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwdsp=NULL,gwdsp.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("gwdsp#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    m$terms[[termnumber]] <- list(name="gwdsp", soname="statnet",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"gwdsp") # must match params$gwdsp above
  }else{ # fixed == TRUE
    m$terms[[termnumber]] <- list(name="gwdsp", soname="statnet",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,paste("gwdsp.fixed.",alpha,sep=""))
  }
  m
}

InitErgm.gwesp<-function(nw, m, arglist, initialfit=FALSE, ...) {
  ergm.checkdirected("gwesp", is.directed(nw), requirement=FALSE)
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
    m$terms[[termnumber]] <- list(name="esp", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(gwesp=NULL,gwesp.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("esp#",d,sep=""))
  }else if (initialfit && !fixed) { # First pass to get MPLE coefficient
    m$terms[[termnumber]] <- list(name="gwesp", soname="statnet",
                                  inputs=c(0, 1, 1, alpha))
    m$coef.names<-c(m$coef.names,"gwesp") # Must match params$gwesp above
  }else{ # fixed == TRUE
    m$terms[[termnumber]] <- list(name="gwesp", soname="statnet",
                                  inputs=c(0, 1, 1, alpha))
    m$coef.names<-c(m$coef.names,paste("gwesp.fixed.",alpha,sep=""))
  }
  m
}

InitErgm.r0<-function(nw, m, arglist, ...) {
  ergm.checkdirected("r0", is.directed(nw), requirement=FALSE)
  a <- ergm.checkargs("r0", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  if(is.bipartite(nw)){
   m$terms[[termnumber]] <- list(name="birnought", soname="statnet",
                                inputs=c(0, 1, 0))
  }else{
   m$terms[[termnumber]] <- list(name="rnought", soname="statnet",
                                inputs=c(0, 1, 0))
  }
  m$coef.names<-c(m$coef.names,"r0")
  m
}

InitErgm.hamming<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("hamming", arglist,
    varnames = c("x","attrname"),
    vartypes = c("matrixnetwork", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(FALSE, FALSE))
  attach(a)
  x<-a$x;attrname<-a$attrname
  if(is.network(x)){
    xm<-as.matrix.network(x,matrix.type="edgelist",attrname)
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
  if (is.null(xm) || NCOL(xm)!=2){
    stop("hamming() requires an edgelist")
  }
  termnumber <- 1 + length(m$terms)
# There is 1 input parameter before the covariate vector, so input
# element 1 is set to 1 (although in this case, input element 1
# is actually arbitrary since d_hamming ignores the value of inp->attrib).
  m$terms[[termnumber]] <- list(name = "hamming",  soname="statnet",
                                inputs = c(0, 1, NROW(xm)*2+1, NROW(xm), as.integer(xm)))
  m$coef.names<-c(m$coef.names, paste("hamming",x,sep="."))
  m
}

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
        cat("Warning: These idegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
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
        cat("Warning: These idegree terms have extreme counts and will be dropped:\n")
        cat(d[midegree], "\n", fill=T)
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
    m$terms[[termnumber]] <- list(name="idegree", soname="statnet",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("idegree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="idegree_w_homophily", soname="statnet",
                                  inputs=c(0, length(d), 
                                           length(d) + length(nodecov), 
                                           d, nodecov),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("ideg", d, ".homophily.",
                                       attrname, sep=""))
  } else {
    if(ncol(du)==0) {return(m)}
    #  No covariates here, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="idegree_by_attr", soname="statnet",
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


InitErgm.istar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("istar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("istar", arglist,
    varnames = c("k", "attrname"),
    vartypes = c("numeric", "character"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, FALSE))
  attach(a)
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
        cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        k <- k[!istarattr] 
      }
    }
  }else{
    if(drop){
      mistar <- paste("c(",paste(k,collapse=","),")",sep="")
      mistar <- summary(as.formula(paste('nw ~ istar(',mistar,')',sep="")),
                        drop=FALSE) == 0
      if(any(mistar)){
        cat(paste("Warning: There are no order", k[mistar],"stars.\n"))
        dropterms <- paste("istar", k[mistar],sep="")
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        k <- k[!mistar] 
      }
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="istar", soname="statnet",
                                  inputs=c(lk, lk, lk+length(nodecov),
                                    k, nodecov))
    m$coef.names<-c(m$coef.names,paste("istar",k,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="istar", soname="statnet",
                                  inputs=c(0, lk, lk, k))
    m$coef.names<-c(m$coef.names,paste("istar",k,sep=""))
  }
  m
}

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
      cat(paste("Warning: There are no isolates.\n"))
      cat(paste("To avoid degeneracy the term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="isolates", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"isolates")
  m
}

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
        cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        k <- k[!kstarattr] 
      }
    }
  }else{
    if(drop){
      mkstar <- paste("c(",paste(k,collapse=","),")",sep="")
      mkstar <- summary(as.formula(paste('nw ~ kstar(',mkstar,')',sep="")),
                        drop=FALSE) == 0
      if(any(mkstar)){
        cat(paste("Warning: There are no order", k[mkstar],"stars.\n"))
        dropterms <- paste("kstar", k[mkstar],sep="")
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        k <- k[!mkstar] 
      }
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)){
    m$terms[[termnumber]] <- list(name="kstar", soname="statnet",
                                  inputs=c(lk, lk, lk+length(nodecov),
                                    k, nodecov))
    m$coef.names<-c(m$coef.names,paste("kstar",k,".",attrname,sep=""))
  }else{
    m$terms[[termnumber]] <- list(name="kstar", soname="statnet",
                                  inputs=c(0, lk, lk, k))
    m$coef.names<-c(m$coef.names,paste("kstar",k,sep=""))
  }
  m
}

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
  m$terms[[termnumber]] <- list(name = "localtriangle", soname="statnet", 
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
      cat(paste("Warning: There are no mutual ties.\n"))
      cat(paste("To avoid degeneracy the 'mutual' term has been dropped.\n"))
      return(m)
    }
    if(nmutual==network.dyadcount(nw)){
      cat(paste("Warning: All dyads have mutual ties!\n"))
      cat(paste("To avoid degeneracy the 'mutual' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="mutual", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"mutual")
  m
}

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
      cat(paste("Warning: There are no asymmetric ties.\n"))
      cat(paste("To avoid degeneracy the 'asymmetric' term has been dropped.\n"))
      return(m)
    }
    if(nasymmetric==network.dyadcount(nw)){
      cat(paste("Warning: All dyads have asymmetric ties!\n"))
      cat(paste("To avoid degeneracy the 'asymmetric' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="asymmetric", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"asymmetric")
  m
}

InitErgm.nodematch<-function (nw, m, arglist, drop=TRUE, ...) {
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
    mixmat <- mixingmatrix(nw,attrname)
    mixmat <- mixmat[-NROW(mixmat),-NROW(mixmat)]
    ematch  <- diag(mixmat)
    if(diff){
      offematch <- apply(mixmat,1,sum)+apply(mixmat,2,sum)-2*ematch
#     Diagonals or off-diagonals are zero
      mu <- ematch==0 | offematch==0
      mu[is.na(mu)] <- FALSE
      if(any(mu)){
        dropterms <- paste(paste("nodematch",attrname,sep="."),u[mu],sep="")
        cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        u <- u[!mu] 
        ui <- ui[!mu] 
      }
    }else{
      offematch <- sum(mixmat)-sum(ematch)
#     Diagonals or off-diagonals are zero
      mu <- sum(ematch)==0 | offematch==0
      mu[is.na(mu)] <- FALSE
      if(mu){
        cat(paste("Warning: The number of matching dyads is extreme.\n"))
        dropterms <- paste("nodematch",attrname,sep=".")
        cat(paste("To avoid degeneracy the term",dropterms,
                  "has been dropped.\n"))
        return(m)
      }
    }
  }
  termnumber<-1+length(m$terms)  
  if (!diff) {
    m$terms[[termnumber]] <- list(name="nodematch", soname="statnet",
                                  inputs=c(0,1,length(nodecov),nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("nodematch",attrname,sep="."))
  } else {
        #  Number of input parameters before covariates equals number of
        #  unique elements in nodecov, namely length(u), so that's what
        #  input element 1 equals
    m$terms[[termnumber]] <- list(name="nodematch", soname="statnet",
                                  inputs=c(length(ui), length(ui),
                                    length(ui)+length(nodecov),
                                    ui, nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("nodematch",
                                       attrname, u, sep="."))
  }
  m
}

InitErgm.nodemain<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("nodemain", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("nodemain",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "nodemain")
  if(!is.numeric(nodecov)){
    stop("nodemain() attribute must be numeric", call.=FALSE) 
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodemain", soname="statnet",
                                inputs=c(0,1,length(nodecov),nodecov),
                                dependence=FALSE)
  m
}

InitErgm.nodemix<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("nodemix", arglist,
    varnames = c("attrname","contrast"),
    vartypes = c("character","logical"),
    defaultvalues = list(NULL,FALSE),
    required = c(TRUE,FALSE))
  attach(a)
  attrname<-a$attrname
  contrast<-a$contrast
  if(is.bipartite(nw)){
#
#   So undirected network storage but directed mixing
#
    nodecov <- get.node.attr(nw, attrname, "mix")
    mixmat <- mixingmatrix(nw,attrname)
    mixmat <- mixmat[-NROW(mixmat),-NROW(mixmat)]
    u <- cbind(as.vector(row(mixmat)), 
               as.vector(col(mixmat)))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}
  #
  #   Recode to numeric if necessary
  #
    namescov <- sort(unique(nodecov))
    nodecov <- match(nodecov,namescov)
    if (length(nodecov)==1)
        stop ("Argument to mix() has only one value", call.=FALSE)
  #
  # Check for degeneracy
  #
    if(drop){
     ematch <- mixmat[u]
     mu <- ematch==0
     mu[is.na(mu)] <- FALSE
     if(any(mu)){
      dropterms <- paste(paste("mix",attrname,sep="."),
        apply(u,1,paste,collapse="")[mu],sep="")
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
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
    m$terms[[termnumber]] <- list(name="mix", soname="statnet",
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
    uui <- uui[upper.tri(uui,diag=TRUE)]
    urm <- t(sapply(ui,rep,length(ui)))   #This is the reverse of what you'd
    urm <- urm[upper.tri(urm,diag=TRUE)]  #expect for r/c, but it's correct
    ucm <- sapply(ui,rep,length(ui))
    ucm <- ucm[upper.tri(ucm,diag=TRUE)]
    uun <- outer(u,u,paste,sep=".")
    uun <- uun[upper.tri(uun,diag=TRUE)]
    if (length(u)==1)
      stop ("Argument to nodemix() has only one value", call.=FALSE)
    if(drop){
      mixmat <- mixingmatrix(nw,attrname)
      mixmat <- mixmat[-NROW(mixmat),-NROW(mixmat)]
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
        cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        cat(paste("To avoid degeneracy the terms",dropterms,
                  "have been dropped.\n"))
        if (sum(!mu)<=1){
          stop ("One or less values of the attribute to nodemix()", call.=FALSE)
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
    m$terms[[termnumber]] <- list(name="nodemix", soname="statnet",
                                  inputs=c(2*length(uui), length(uui),
                                    2*length(uui)+length(nodecov),
                                    urm, ucm, nodecov),
                                  dependence=FALSE)
    m$coef.names<-c(m$coef.names,paste("mix",attrname, uun, sep="."))
  }
  m
}

InitErgm.mix<-InitErgm.nodemix 

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
      cat(paste("Warning: There are no in ties for the vertex", 
                d[is.na(mdegrees)],"\n"))
      dropterms <- paste("receiver", d[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the term",dropterms,
                "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="receiver", soname="statnet",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("receiver",d,sep=""))
  m
}
InitErgm.receivercov<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("receivercov", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("receivercov",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "receivercov")
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="receivercov", soname="statnet",
                                inputs=c(0,1,length(nodecov),nodecov))
  m
}

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
      cat(paste("Warning: There are no out ties for the vertex", 
                d[is.na(mdegrees)],"\n"))
      dropterms <- paste("sender", d[is.na(mdegrees)],sep="")
      cat(paste("To avoid degeneracy the term",dropterms,
                "have been dropped.\n"))
      d <- degrees[mdegrees[!is.na(mdegrees)]] 
    }
  }
  ld<-length(d)
  if(ld==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="sender", soname="statnet",
                                inputs=c(0, ld, ld, d))
  m$coef.names<-c(m$coef.names,paste("sender",d,sep=""))
  m
}

InitErgm.sendercov<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("sendercov", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
  m$coef.names<-c(m$coef.names, paste("sendercov",attrname,sep="."))
  nodecov <- get.node.attr(nw, attrname, "sendercov")
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="sendercov", soname="statnet",
                                inputs=c(0,1,length(nodecov),nodecov))
  m
}


InitErgm.nodefactor<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("nodefactor", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
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
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodefactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodefactor", soname="statnet",
                                inputs=c(lu-1, lu-1, lu-1+length(nodecov),
                                         ui[-1], nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("nodefactor",
                                      attrname, paste(u[-1]), sep="."))
  m
}


InitErgm.nodeifactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("nodeifactor", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("nodeifactor", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
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
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodeifactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodeifactor", soname="statnet",
                                inputs=c(lu-1, lu-1, lu-1+length(nodecov),
                                         ui[-1], nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("nodeifactor",
                                      attrname, paste(u[-1]), sep="."))
  m
}

InitErgm.nodeofactor<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("nodeofactor", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("nodeofactor", arglist,
    varnames = c("attrname"),
    vartypes = c("character"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  attrname<-a$attrname
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
      cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      u<-u[nfc>0]
      ui<-ui[nfc>0]
    }
  }
  lu <- length(ui)
  if (lu==1){
    stop ("Argument to nodeofactor() has only one value", call.=FALSE)
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nodeofactor", soname="statnet",
                                inputs=c(lu-1, lu-1, lu-1+length(nodecov),
                                         ui[-1], nodecov), dependence=FALSE)
  # smallest value of u is "control group"
  m$coef.names<-c(m$coef.names, paste("nodeofactor",
                                      attrname, paste(u[-1]), sep="."))
  m
}

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
        cat("Warning: These odegree terms have extreme counts and will be dropped:\n")
        cat(dropterms, "\n", fill=T)
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
        cat("Warning: These odegree terms have extreme counts and will be dropped:\n")
        cat(d[modegree], "\n", fill=T)
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
    m$terms[[termnumber]] <- list(name="odegree", soname="statnet",
                                  inputs=c(0, length(d), length(d), d),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,paste("odegree",d,sep=""))
  } else if (homophily) {
    if(length(d)==0){return(m)}
    m$terms[[termnumber]] <- list(name="odegree_w_homophily", soname="statnet",
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
    m$terms[[termnumber]] <- list(name="odegree_by_attr", soname="statnet",
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

InitErgm.ostar<-function(nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ostar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ostar", arglist,
    varnames = c("k"),
    vartypes = c("numeric"),
    defaultvalues = list(NULL),
    required = c(TRUE))
  attach(a)
  k<-a$k
  if(drop){
    ks <- tabulate(tabulate(as.matrix.network.edgelist(nw)[,1]))
    mistar <- NULL
    for(cdeg in k){
      aaa <- (cdeg:length(ks))
      aaa <- aaa[ks[aaa] > 0]
      mistar <- c(mistar, sum(choose(aaa,cdeg),na.rm=TRUE) == 0)
    }
    if(any(mistar)){
      cat(paste("Warning: There are no order", k[mistar],"stars.\n"))
      dropterms <- paste("istar", k[mistar],sep="")
      cat(paste("To avoid degeneracy the terms",dropterms,
                "have been dropped.\n"))
      k <- ks[mistar[!mistar]] 
    }
  }
  lk<-length(k)
  if(lk==0){return(m)}
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="ostar", soname="statnet",
                                inputs=c(0, lk, lk, k))
  m$coef.names<-c(m$coef.names,paste("ostar",k,sep=""))
  m
}

InitErgm.smalldiff<-function (nw, m, arglist, ...) {
  a <- ergm.checkargs("smalldiff", arglist,
    varnames = c("attrname", "cutoff"),
    vartypes = c("character", "numeric"),
    defaultvalues = list(NULL, NULL),
    required = c(TRUE, TRUE))
  attach(a)
  if (length(cutoff)>1)
    stop("cutoff for smalldiff() must be a scalar.", call.=FALSE)
  m$coef.names<-c(m$coef.names,paste("smalldiff.",
                                     attrname, cutoff, sep=""))
  nodecov <- get.node.attr(nw, attrname, "smalldiff")
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="smalldiff", soname="statnet",
                                inputs=c(1, 1, 1+length(nodecov),
                                  cutoff, nodecov), dependence=FALSE)
  m
}

InitErgm.triangle<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("triangle", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
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
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the terms",dropterms,
                    "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("triangle",attrname,sep="."),sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the term",dropterms,
                    "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="triangle", soname="statnet",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("triangle",attrname,sep="."))
    } else {
      m$terms[[termnumber]] <- list(name="triangle", soname="statnet",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("triangle",
                                         attrname, u, sep="."))
    }
  }else{
    m$terms[[termnumber]] <- list(name="triangle", soname="statnet",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"triangle")
  }
  m
}
InitErgm.triangles<-InitErgm.triangle

InitErgm.tripercent<-function (nw, m, arglist, drop=TRUE, ...) {
  a <- ergm.checkargs("tripercent", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
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
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the terms",dropterms,
                    "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }else{
          dropterms <- paste(paste("tripercent",attrname,sep="."),sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the term",dropterms,
                    "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="tripercent", soname="statnet",
                                    inputs=c(1, 1, 1+length(nodecov),
                                      1, nodecov))
      m$coef.names<-c(m$coef.names,paste("tripercent",attrname,sep="."))
    } else {
      m$terms[[termnumber]] <- list(name="tripercent", soname="statnet",
                                    inputs=c(length(ui), length(ui),
                                      length(ui)+length(nodecov),
                                      ui, nodecov))
      m$coef.names<-c(m$coef.names,paste("tripercent",
                                         attrname, u, sep="."))
    }
  }else{
    m$terms[[termnumber]] <- list(name="tripercent", soname="statnet",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"tripercent")
  }
  m
}

InitErgm.ttriad<-function (nw, m, arglist, drop=TRUE, ...) {
  ergm.checkdirected("ttriad", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("ttriad", arglist,
    varnames = c("attrname", "diff"),
    vartypes = c("character", "logical"),
    defaultvalues = list(NULL, FALSE),
    required = c(FALSE, FALSE))
  attach(a)
  termnumber<-1+length(m$terms)
  if(!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname, "ttriad")
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    ui <- seq(along=u)
    if (length(u)==1)
      stop ("Attribute given to ttriad() has only one value", call.=FALSE)
    if(drop){
      triattr <- summary(as.formula(paste('nw ~ ttriad(','"',attrname,
                                          '",diff=',diff,')',sep="")),
                         drop=FALSE) == 0
      if(diff){
        if(any(triattr)){
          dropterms <- paste(paste("ttriad",attrname,sep="."),
                             u[triattr],sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the terms",dropterms,
                    "have been dropped.\n"))
          u <- u[!triattr] 
          ui <- ui[!triattr] 
        }
      }else{
        if(triattr){
          dropterms <- paste(paste("ttriad",attrname,sep="."),sep="")
          cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
          cat(paste("To avoid degeneracy the term",dropterms,
                    "have been dropped.\n"))
        }
      }
    }
    if (!diff) {
      m$terms[[termnumber]] <- list(name="ttriad", soname="statnet",
                                    inputs=c(0,1,length(nodecov),nodecov))
      m$coef.names<-c(m$coef.names,paste("ttriad",attrname,sep="."))
     } else {
       m$terms[[termnumber]] <- list(name="ttriad", soname="statnet",
                                     inputs=c(length(ui), length(ui),
                                       length(ui)+length(nodecov),
                                       ui, nodecov))
       m$coef.names<-c(m$coef.names,paste("ttriad",
                                          attrname, u, sep="."))
     }
  }else{
    m$terms[[termnumber]] <- list(name="ttriad", soname="statnet",
                                  inputs=c(0,1,0))
    m$coef.names<-c(m$coef.names,"ttriad")
  }
  m
}

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
      cat(paste("Warning: There are no hiertriad ties.\n"))
      cat(paste("To avoid degeneracy the 'hiertriad' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="hiertriad", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"hiertriad")
  m
}

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
      cat(paste("Warning: There are no intransitive triads\n"))
      cat(paste("To avoid degeneracy the 'intransitive' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="intransitive", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"intransitive")
  m
}

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
      cat(paste("Warning: There are no transitive triads\n"))
      cat(paste("To avoid degeneracy the 'transitive' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="transitive", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"transitive")
  m
}

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
      cat(paste("Warning: There are no intransitive triads\n"))
      cat(paste("To avoid degeneracy the 'intransitivity' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="intransitivity", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"intransitivity")
  m
}

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
      cat(paste("Warning: There are no transitive triads\n"))
      cat(paste("To avoid degeneracy the 'transitivity' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="transitivity", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"transitivity")
  m
}

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
      cat(paste("Warning: There are no hiertriaddegree ties.\n"))
      cat(paste("To avoid degeneracy the 'hiertriaddegree' term has been dropped.\n"))
      return(m)
    }
  }
  m$terms[[termnumber]] <- list(name="hiertriaddegree", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"hiertriaddegree")
  m
}

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
    m$terms[[termnumber]] <- list(name="berninhom", soname="statnet",
                                inputs=c(0, nstat, 0, 0),
                                dependence=TRUE,
                                params=list(spatial.pb=pb,
                                spatial.alpha=alpha,spatial.gamma=gamma),
                                map=map, gradient=gradient, 
                                eta.cov=d)
    m$coef.names<-c(m$coef.names,paste("spatial#",1:nstat,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="spatial", soname="statnet",
                                inputs=c(0, 1, 3+nstat, c(pb, alpha, gamma, d)),
                                dependence=TRUE)
    m$coef.names<-c(m$coef.names,"spatial.pb")
  }
  m
}


#
# Depreciated terms
#

InitErgm.match<-InitErgm.nodematch
InitErgm.nodecov <- InitErgm.nodemain

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
    m$terms[[termnumber]] <- list(name="degree", soname="statnet",
                                  inputs=c(0, ld, ld, d),
                                  params=list(geodegree=NULL,
                                    geodegree.alpha=alpha),
                                  map=map, gradient=gradient)
    m$coef.names<-c(m$coef.names,paste("geodegree#",d,sep=""))
  }else{
    termnumber<-1+length(m$terms)
    m$terms[[termnumber]] <- list(name="geodegree", soname="statnet",
                                  inputs=c(0, 1, length(alpha), alpha))
    m$coef.names<-c(m$coef.names,"geodegree")
  }
  m
}

#InitErgm.balance<-function (nw, m, arglist, ...) {
#  warning(paste("balance() model term is only available in the ",
#                "'statnet' package\n", 
#                "Ingoring this term.\n", sep = ""), call. = FALSE)
#  m
#}

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
        cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
        cat(paste("To avoid degeneracy the terms",dropterms,"have been dropped.\n"))
        u <- u[!triattr] 
        ui <- ui[!triattr] 
       }
      }else{
       if(triattr){
         dropterms <- paste(paste("balance",attrname,sep="."),sep="")
         cat(paste("Warning: The count of", dropterms, "is extreme.\n"))
         cat(paste("To avoid degeneracy the term",dropterms,"have been dropped.\n"))
       }
      }
     }
     if (!diff) {
#     No parameters before covariates here, so input element 1 equals 0
      m$terms[[termnumber]] <- list(name="balance", soname="statnet",
                                    inputs=c(0,1,length(nodecov),nodecov),
                                    dependence=TRUE)
      m$coef.names<-c(m$coef.names,paste("balance",attrname,sep="."))
     } else {
      #  Number of input parameters before covariates equals number of
      #  unique elements in nodecov, namely length(u), so that's what
      #  input element 1 equals
      m$terms[[termnumber]] <- list(name="balance", soname="statnet",
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
#       cat(paste("Warning: There are no balanced triads\n"))
#       cat(paste("To avoid degeneracy the balance term has been dropped.\n"))
#    }
#   }
#   No covariates, so input element 1 is arbitrary
    m$terms[[termnumber]] <- list(name="balance", soname="statnet",
                                  inputs=c(0,1,0),
                                  dependence=TRUE)
    m$coef.names<-c(m$coef.names,"balance")
   }
   m
}


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
      cat(paste("Warning: There are no simmelian triads\n"))
      cat(paste("To avoid degeneracy the 'simmelian' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelian==network.edgecount(nw)*network.size*0.5){
      cat(paste("Warning: All triads are simmelian!\n"))
      cat(paste("To avoid degeneracy the 'simmelian' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="simmelian", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"simmelian")
  m
}

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
      cat(paste("Warning: There are no simmelianties ties\n"))
      cat(paste("To avoid degeneracy the 'simmelianties' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelianties==network.edgecount(nw)){
      cat(paste("Warning: All ties have simmelianties ties!\n"))
      cat(paste("To avoid degeneracy the 'simmelianties' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="simmelianties", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"simmelianties")
  m
}

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
      cat(paste("Warning: There are no nearsimmelian triads\n"))
      cat(paste("To avoid degeneracy the 'nearsimmelian' term has been dropped.\n"))
      return(m)
    }
    if(nsimmelian==network.dyadcount(nw)*network.size(nw)*0.5){
      cat(paste("Warning: All dyads have nearsimmelian triads!\n"))
      cat(paste("To avoid degeneracy the 'nearsimmelian' term has been dropped.\n"))
      return(m)
    }
  }
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="nearsimmelian", soname="statnet",
                                inputs=c(0,1,0))
  m$coef.names<-c(m$coef.names,"nearsimmelian")
  m
}

InitErgm.icvar<-function(nw, m, arglist, ...) {
  ergm.checkdirected("icvar", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("icvar", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
    termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="icvar", soname="statnet",
                                inputs=c(0, 1, 0),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,"icvar")
  m
}

InitErgm.idc<-function(nw, m, arglist, ...) {
  ergm.checkdirected("idc", is.directed(nw), requirement=TRUE)
  a <- ergm.checkargs("idc", arglist,
    varnames = NULL,
    vartypes = NULL,
    defaultvalues = list(),
    required = NULL)
  termnumber<-1+length(m$terms)
  m$terms[[termnumber]] <- list(name="idc", soname="statnet",
                                inputs=c(0, 1, 0),
                                dependence=TRUE)
  m$coef.names<-c(m$coef.names,"idc")
  m
}
