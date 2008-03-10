#################################################################
# InitErgmTerm functions (new, easier-to-write version)
# The old InitErgm functions should still work.
#
# INPUT:
# Each InitErgmTerm function takes two arguments, nw and arglist,
# which are automatically supplied by ergm.getmodel.  There may be
# other arguments passed by ergm.getmodel, so each InitErgmTerm
# function must also include the ... argument in its list.
#
# OUTPUT:
# Each InitErgmTerm function should return a list.  
#    REQUIRED LIST ITEMS:
#          name: Name of the C function that produces the change
#                statistics.  (Note:  The name will have "d_" 
#                prepended.  For example, the C function for the
#                absdiff change statistics is called "d_absdiff"
#                even though InitErgmTerm.absdiff only returns
#                names = "absdiff")
#    coef.names: Vector of names for the coefficients (parameters)
#                as they will be reported in the output.
#
#    OPTIONAL LIST ITEMS:
#        inputs: Vector of inputs (of type double) that the
#                d_xxx function will require.  Default is NULL.  This vector
#                may have an attribute named "ParamsBeforeCov", which is the
#                number that used to be the old Element 1 (number of input
#                parameters BEFORE the beginning of the covariate vector)                                                         
#                when using the old InitErgm specification; see the comment
#                at the top of the InitErgm.R file for details.  This 
#                ParamsBeforeCov value is only necessary for compatibility 
#                with some of the existing d_xxx changestatistic functions.  
#        soname: This is the (text) name of the package containing the C function
#                called d_[name].  Default is "ergm"
#    dependence: Logical variable telling whether addition of this term to
#                the model makes the model into a dyadic dependence model.
#                If none of the terms sets dependence==TRUE, then the model
#                is assumed to be a dyadic independence model, which means
#                that the pseudolikelihood estimate coincides with the
#                maximum likelihood estimate.  Default value:  TRUE
#  emptynwstats: Vector of values (if nonzero) for the statistics evaluated
#                on the empty network.  If all are zero for this term, this
#                argument may be omitted.  Example:  If the degree0 term is
#                among the statistics, this argument is necessary because
#                degree0 = number of nodes for the empty network.
#        params: For curved exponential family model terms only: This argument 
#                is a list:  Each item in the list should be named with the
#                corresponding parameter name (one or more of these will
#                probably coincide with the coef.names used when
#                initialfit=TRUE; the initial values of such parameter values
#                will be set by MPLE, so their values in params are ignored.)
#                Any parameter not having its initial value set by MPLE
#                should be given its initial value in this params list.
#           map: For curved exponential family model terms only: A function 
#                that gives the map from theta (the canonical
#                parameters associated with the statistics for this term)
#                to eta (the corresponding curved parameters).  The length
#                of eta is the same as the length of the params list above.
#                This function takes two args:  theta and length(eta).
#      gradient: For curved exponential family model terms only: A function 
#                that gives the gradient of the eta map above.
#                If theta has length p and eta has length q, then gradient
#                should return a p by q matrix.
#                This function takes two args:  theta and length(eta).
#

# Prototype InitErgmTerm functions
#########################################################
InitErgmTerm.absdiff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  ### Construct the output list
  list(name="absdiff",                                     #name: required
       coef.names = paste("absdiff", a$attrname, sep="."), #coef.names: required
       inputs = nodecov,  # We need to include the nodal covariate for this term
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

########################################################
InitErgmTerm.absdiffcat <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","base"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list(NULL,NULL),
                      required = c(TRUE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=TRUE)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(a$base)) u <- u[-(a$base)]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  u2 <- u[!is.na(u)]
  ### Construct the output list
  inputs <- c(u2, NAsubstitute, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(u2)+1 # See comment at top of file
  list(name="absdiffcat",                                  #name: required
       coef.names = paste("absdiff", a$attrname, u, sep="."), #coef.names: required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

#########################################################
InitErgmTerm.altkstar <- function(nw, arglist, initialfit=FALSE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=FALSE, bipartite=NULL,
                      varnames = c("lambda","fixed"),
                      vartypes = c("numeric","logical"),
                      defaultvalues = list(1,FALSE),
                      required = c(FALSE,FALSE))
  ### Process the arguments
  if(!initialfit && !a$fixed){ # This is a curved exponential family model
    d <- 1:(network.size(nw)-1)
    map <- function(x,n,...) {
      i <- 1:n
      x[1]*(x[2]*((1-1/x[2])^i + i) - 1)
    }
    gradient <- function(x,n,...) {
      i <- 1:n
      rbind(x[2]*((1-1/x[2])^i + i) - 1,
            x[1]*(i - 1 + (x[2]*x[2]-x[2]+i)*((1-1/x[2])^(i-1))/(x[2]*x[2])))
    }
    ### Construct the output list
    outlist <- list(name="degree",                 #name: required
       coef.names = paste("altkstar#", d, sep=""), #coef.names: required
       inputs = d, map=map, gradient=gradient,
       params=list(altkstar=NULL, altkstar.lambda=a$lambda)
       )
  } else {
    outlist <- list (name="altkstar",                      #name: required
       coef.names = paste("altkstar", a$lambda, sep="."),  #coef.names: required
       inputs=a$lambda
       )
  }
  outlist
}

#########################################################
InitErgmTerm.asymmetric <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attrname", "diff", "keep"),
                      vartypes = c("character", "logical", "numeric"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE))
  ### Process the arguments
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u <- sort(unique(nodecov))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    if (is.null(a$attrname)) {
      n <- network.size(nw)
      ndc <- n * (n-1) / 2 # temporary until network.dyadcount is fixed
      if (extremewarnings(obsstats, maxval=ndc)) {
        return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
      }
    } else {
      ew <- extremewarnings(obsstats)
      u <- u[!ew]
      ui <- ui[!ew]
    }
  }
  ### Construct the output list
  out <- list(name="asymmetric",                      #name: required
              coef.names = "asymmetric"               #coef.names: required
              ) 
  if (!is.null(a$attrname)) {
    if (a$diff) {
      out$coef.names <- paste("asymmetric", a$attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("asymmetric", a$attrname, sep=".")
      out$inputs <- nodecov
    }
  }
  out
}


#########################################################
InitErgmTerm.isolates <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                     varnames = NULL,
                     vartypes = NULL,
                     defaultvalues = list(),
                     required = NULL)
  ### Process the arguments
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    if (extremewarnings(obsstats, maxval=network.size(nw))) {
      return (NULL)  # Do not add this term at all if isolates==0 or n
    }
  }
  ### Construct the output list
  list(name="isolates",                               #name: required
       coef.names = "isolates",                       #coef.names: required
       emptynwstats = network.size(nw) # When nw is empty, isolates=n, not 0
       )                                                               
}

#########################################################
InitErgmTerm.mutual<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attrname", "diff", "keep"),
                      vartypes = c("character", "logical", "numeric"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(FALSE, FALSE, FALSE))
  ### Process the arguments
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u <- sort(unique(nodecov))
    if (!is.null(a$keep)) {
      u <- u[a$keep]
    }
    #   Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    # All of the "nomatch" should be given unique IDs so they never match:
    dontmatch <- nodecov==(length(u)+1)
    nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
    ui <- seq(along=u)
  }
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    if (is.null(a$attrname)) {
      n <- network.size(nw)
      ndc <- n * (n-1) / 2 # temporary until network.dyadcount is fixed
      if (extremewarnings(obsstats, maxval=ndc)) {
        return(NULL) # In this case the obs nw has 0 or n(n-1)/2 asymmetric dyads
      }
    } else {
      ew <- extremewarnings(obsstats)
      u <- u[!ew]
      ui <- ui[!ew]
    }
  }
  ### Construct the output list
  out <- list(name="mutual",                      #name: required
              coef.names = "mutual"               #coef.names: required
              ) 
  if (!is.null(a$attrname)) {
    if (a$diff) {
      out$coef.names <- paste("mutual", a$attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("mutual", a$attrname, sep=".")
      out$inputs <- nodecov
    }
  }
  out
}

#########################################################
InitErgmTerm.nodefactor<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    u <- u[!ew]
    ui <- ui[!ew]
  }
  ### Construct the output list
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodefactor",                                        #required
       coef.names = paste("nodefactor", a$attrname, u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  

#########################################################
InitErgmTerm.nodeifactor<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = c("attrname", "base"),          
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    u <- u[!ew]
    ui <- ui[!ew]
  }
  ### Construct the output list
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeifactor",                                        #required
       coef.names = paste("nodeifactor", a$attrname, u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  

#########################################################
InitErgmTerm.nodematch<-InitErgmTerm.match<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("attrname", "diff", "keep"),
                      vartypes = c("character", "logical", "numeric"),
                      defaultvalues = list(NULL, FALSE, NULL),
                      required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(nodecov))
  if (!is.null(a$keep)) {
    u <- u[a$keep]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  # All of the "nomatch" should be given unique IDs so they never match:
  dontmatch <- nodecov==(length(u)+1)
  nodecov[dontmatch] <- length(u) + (1:sum(dontmatch))
  ui <- seq(along=u)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    u <- u[!ew]
    ui <- ui[!ew]
  }
  ### Construct the output list
  if (a$diff) {
    coef.names <- paste("nodematch", a$attrname, u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", a$attrname, sep=".")
    inputs <- nodecov
  }
  list(name="nodematch",                                 #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

#########################################################
InitErgmTerm.nodeofactor<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, 
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, 1),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
  }
  #   Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  ui <- seq(along=u)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    u <- u[!ew]
    ui <- ui[!ew]
  }
  ### Construct the output list
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeofactor",                                        #required
       coef.names = paste("nodeofactor", a$attrname, u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  


