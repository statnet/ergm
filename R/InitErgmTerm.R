#  File ergm/R/InitErgmTerm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
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
#        inputs: Vector of (double-precision numeric) inputs that the 
#                changestat function called d_<name> will require
#                (see WHAT THE C CHANGESTAT FUNCTION RECEIVES below).
#                The default is NULL; no inputs are required.  But it MUST
#                be a vector!  Thus, if some of the inputs are, say, matrices,
#                they must be "flattened" to vectors; if some are categorical
#                character-valued variables, they must be converted to numbers.
#                Optionally, the inputs vector may have an attribute named 
#                "ParamsBeforeCov", which is the
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
# WHAT THE C CHANGESTAT FUNCTION RECEIVES:
#                The changestat function, written in C and called d_<name>,
#                where <name> is the character string passed as the required
#                output item called "name" (see above), will have access to
#                the vector of double-precision values created by the 
#                InitErgmTerm function as the optional output item called
#                "inputs".  This array will be called INPUT_PARAMS in the C
#                code and its entries may accessed as INPUT_PARAMS[0],
#                INPUT_PARAMS[1], and so on.  The size of the INPUT_PARAMS 
#                array is equal to N_INPUT_PARAMS, a value which is 
#                automatically set for you and which is available inside the
#                C function.  Thus INPUT_PARAMS[N_INPUT_PARAMS-1] is the last
#                element in the vector. Note in particular that it is NOT 
#                necessary to add the number of inputs to the "inputs" vector
#                since this is done automatically.
#
# Prototype InitErgmTerm functions
############################## InitErgmTerm functions:  A
#########################################################
InitErgmTerm.absdiff <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attrname","pow"),
                      vartypes = c("character","numeric"),
                      defaultvalues = list(NULL,1),
                      required = c(TRUE,FALSE))
  ### Process the arguments
  nodecov <- get.node.attr(nw, a$attrname)
  ### Construct the list to return
  list(name="absdiff",                                     #name: required
       coef.names = paste(paste("absdiff",if(a$pow!=1) a$pow else "",sep=""), a$attrname, sep="."), #coef.names: required
       inputs = c(a$pow,nodecov),  # We need to include the nodal covariate for this term
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
  ### Construct the list to return
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
    ### Construct the list to return
    outlist <- list(name="degree",                 #name: required
       coef.names = paste("altkstar#", d, sep=""), #coef.names: required
       inputs = d, map=map, gradient=gradient,
       params=list(altkstar=NULL, altkstar.lambda=a$lambda)
       )
  } else {
    if (initialfit) { # coef.names must match "altkstar" from params list above
      coef.names = "altkstar"
    } else { # Not necessary to match; provide more informative label
      coef.names = paste("altkstar", a$lambda, sep=".")
    }
    outlist <- list (name="altkstar",                      #name: required
                     coef.names = coef.names,
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
  ### Construct the list to return
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


############################## InitErgmTerm functions:  B
#########################################################
InitErgmTerm.b1degree <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "byarg"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(a$byarg)) {  # CASE 1:  a$byarg GIVEN
    nodecov <- get.node.attr(nw, a$byarg)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    du <- matrix(du[,!ew], nrow=2) # Drop any zero-obs-value rows
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov[1:nb1]==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b1degree_by_attr"
    coef.names <- paste("b1deg", du[1,], ".", a$byarg, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$byarg NOT GIVEN
    a$d <- a$d[!ew] # Drop any zero-obs-value values
    name <- "b1degree"
    coef.names <- paste("b1deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- nb1
    }
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats)
}

#########################################################
InitErgmTerm.b1star <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "ostar"
    coef.names <- paste("b1star", a$k, ".", a$attrname, sep="")
    inputs <- c(a$k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- length(a$k)
  } 
  else {
    name <- "ostar"
    coef.names <- paste("b1star",a$k,sep="")
    inputs <- a$k
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs)
}

#########################################################
InitErgmTerm.b1starmix <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname", "base", "diff"),
                       vartypes = c("numeric", "character", "numeric", "logical"),
                       defaultvalues = list(NULL, NULL, NULL, TRUE),
                       required = c(TRUE, TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  nodecov <- get.node.attr(nw, a$attrname)
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { stop("Only a single scalar k may be used with each b1starmix term") }
  b1namescov <- sort(unique(nodecov[1:nb1]))
  b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
  b1nodecov <- match(nodecov[1:nb1],b1namescov)
  b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
  namescov <- u[c(b1namescov, b2namescov)]
  nr <- length(b1namescov)
  nc <- length(b2namescov)
  nodecov <- c(b1nodecov, b2nodecov + nr)
  if (a$diff) {
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base,] }
    u <- u[!ew,]
    name <- "b1starmix"
    coef.names <- paste("b1starmix", a$k, a$attrname,
                        apply(matrix(namescov[u],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- 1:nr
    if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base] }
    u <- u[!ew]
    name <- "b1starmixhomophily"
    coef.names <- paste("b1starmix", a$k, a$attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs)
}

#########################################################
InitErgmTerm.b1twostar <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("b1attrname", "b2attrname", "base"),
                       vartypes = c("character", "character", "numeric"),
                       defaultvalues = list(NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  b1nodecov <- get.node.attr(nw, a$b1attrname)[1:nb1]
  b1u<-sort(unique(b1nodecov))
  if(any(is.na(b1nodecov))){ b1u<-c(b1u,NA) }
  if(is.null(a$b2attrname)) { a$b2attrname <- a$b1attrname }
  b2nodecov <- get.node.attr(nw, a$b2attrname)[(1+nb1):n]
  b2u<-sort(unique(b2nodecov))
  if(any(is.na(b2nodecov))){b2u<-c(b2u,NA)}
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  nr <- length(b1u)
  nc <- length(b2u)
  u <- cbind(rep(1:nr, nc*nc), rep(rep(1:nc, each=nr), nc), rep(1:nc, each=nc*nr))
  u <- u[u[,2] <= u[,3],]  
  if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base,] }
  u <- u[!ew,]
  coef.names <- paste("b1twostar", a$b1attrname, b1u[u[,1]],  a$b2attrname,
                      apply(matrix(b2u[u[,2:3]],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b1twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]) )
}

#########################################################
InitErgmTerm.b2degree <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "byarg"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(a$byarg)) {  # CASE 1:  a$byarg GIVEN
    nodecov <- get.node.attr(nw, a$byarg)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(a$d)*length(u)
    lu <- length(u)
    du <- rbind(rep(a$d,lu), rep(1:lu, rep(length(a$d), lu)))
    du <- matrix(du[,!ew], nrow=2) # Drop any zero-obs-value rows
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov[(1+nb1):n]==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b2degree_by_attr"
    coef.names <- paste("b2deg", du[1,], ".", a$byarg, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  a$byarg NOT GIVEN
    a$d <- a$d[!ew] # Drop any zero-obs-value values
    name <- "b2degree"
    coef.names <- paste("b2deg", a$d, sep="")
    inputs <- a$d
    emptynwstats <- rep(0, length(a$d))
    if (any(a$d==0)) { # alter emptynwstats
      emptynwstats[a$d==0] <- n-nb1
    }
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats)
}

#########################################################
InitErgmTerm.b2star <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  ### Process the arguments
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(a$attrname)) {
    nodecov <- get.node.attr(nw, a$attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "istar"
    coef.names <- paste("b2star", a$k, ".", a$attrname, sep="")
    inputs <- c(a$k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- length(a$k)
  } 
  else {
    name <- "istar"
    coef.names <- paste("b2star",a$k,sep="")
    inputs <- a$k
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs)
}

#########################################################
InitErgmTerm.b2starmix <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("k", "attrname", "base", "diff"),
                       vartypes = c("numeric", "character", "numeric", "logical"),
                       defaultvalues = list(NULL, NULL, NULL, TRUE),
                       required = c(TRUE, TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  nodecov <- get.node.attr(nw, a$attrname)
  u<-sort(unique(nodecov))
  if(any(is.na(nodecov))){u<-c(u,NA)}
  # Recode to numeric
  nodecov <- match(nodecov,u,nomatch=length(u)+1)
  if (length(a$k) > 1) 
    { stop("Only a single scalar k may be used with each b2starmix term") }
  b1namescov <- sort(unique(nodecov[1:nb1]))
  b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
  b1nodecov <- match(nodecov[1:nb1],b1namescov)
  b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
  namescov <- u[c(b1namescov, b2namescov)]
  nr <- length(b1namescov)
  nc <- length(b2namescov)
  nodecov <- c(b1nodecov, b2nodecov + nr)
  if (a$diff) {
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base,] }
    u <- u[!ew,]
    name <- "b2starmix"
    coef.names <- paste("b2starmix", a$k, a$attrname,
                        apply(matrix(namescov[u[,2:1]],ncol=2), 1,paste,collapse="."), 
                        sep=".")
    inputs <- c(a$k, nodecov, u[,1], u[,2])
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  else {
    u <- nr+(1:nc)
    if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base] }
    u <- u[!ew]
    name <- "b2starmixhomophily"
    coef.names <- paste("b2starmix", a$k, a$attrname, namescov[u], sep=".")
    inputs <- c(a$k, nodecov, u)
    attr(inputs, "ParamsBeforeCov") <- length(a$k) # should be 1
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs)
}

#########################################################
InitErgmTerm.b2twostar <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("b1attrname", "b2attrname", "base"),
                       vartypes = c("character", "character", "numeric"),
                       defaultvalues = list(NULL, NULL, NULL),
                       required = c(TRUE, FALSE, FALSE))
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  b1nodecov <- get.node.attr(nw, a$b1attrname)[1:nb1]
  b1u<-sort(unique(b1nodecov))
  if(any(is.na(b1nodecov))){ b1u<-c(b1u,NA) }
  if(is.null(a$b2attrname)) { a$b2attrname <- a$b1attrname }
  b2nodecov <- get.node.attr(nw, a$b2attrname)[(1+nb1):n]
  b2u<-sort(unique(b2nodecov))
  if(any(is.na(b2nodecov))){b2u<-c(b2u,NA)}
  # Recode to numeric
  b1nodecov <- match(b1nodecov,b1u,nomatch=length(b1u)+1)
  b2nodecov <- match(b2nodecov,b2u,nomatch=length(b2u)+1)
  nr <- length(b1u)
  nc <- length(b2u)
  u <- cbind(rep(1:nc, nr*nr), rep(rep(1:nr, each=nc), nr), rep(1:nr, each=nc*nr))
  u <- u[u[,2] <= u[,3],]  
  if (!is.null(a$base) && !identical(a$base,0)) { u <- u[-a$base,] }
  u <- u[!ew,]
  coef.names <- paste("b2twostar", a$b2attrname, b2u[u[,1]],  a$b1attrname,
                      apply(matrix(b1u[u[,2:3]],ncol=2), 1, paste, collapse="."),
                      sep=".")
  list(name = "b2twostar", coef.names = coef.names, #name and coef.names: required
       inputs = c(b1nodecov, b2nodecov, u[,1], u[,2], u[,3]) )
}

############################## InitErgmTerm functions:  C
#########################################################
InitErgmTerm.cycle <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                     varnames = c("k"),
                     vartypes = c("numeric"),
                     defaultvalues = list(NULL),
                     required = c(TRUE))
  ### Process the arguments
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    a$k <- a$k[!ew]
  }
  if (length(a$k)==0) return(NULL)
  ### Construct the list to return
  list(name="cycle",                            #name: required
       coef.names = paste("cycle", a$k, sep=""),  #coef.names: required
       inputs = c(max(a$k), (2:max(a$k)) %in% a$k)
       )
}

############################## InitErgmTerm functions:  E
#########################################################
InitErgmTerm.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                     varnames = c("x", "attrname"),
                     vartypes = c("matrixnetwork", "character"),
                     defaultvalues = list(NULL, NULL),
                     required = c(TRUE, FALSE))
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix.network(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), 
              as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov",   #name: required
       coef.names = cn,  #coef.names: required
       inputs = inputs,
       dependence=FALSE
       )
}

############################## InitErgmTerm functions:  H
#########################################################
InitErgmTerm.hamming<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm (nw, arglist,
	    varnames = c("x","cov","attrname","defaultweight"),
	    vartypes = c("matrixnetwork","matrixnetwork","character","numeric"),
	    defaultvalues = list(nw, NULL, NULL, NULL),
	    required = c(FALSE, FALSE, FALSE, FALSE))

  ## Process hamming network ##
  if(is.network(a$x)){													# Arg to hamming is a network
    xm<-as.matrix.network(a$x,matrix.type="edgelist",a$attrname)
  }else if(is.character(a$x)){												# Arg to hamming is the name of an attribute in nw
    xm<-get.network.attribute(nw,a$x)
    xm<-as.matrix.network(xm,matrix.type="edgelist")
  }else if(is.null(a$x)){
    xm<-as.matrix.network(nw,matrix.type="edgelist")								# Arg to hamming does not exist; uses nw
  }else{
    xm<-as.matrix(a$x)													# Arg to hamming is anything else; attempts to coerce
  }
  if (is.vector(xm)) xm <- matrix(xm, ncol=2)

  ## Process case without dyadcov (i.e. unweighted) ##
  sc03 <- sys.call(0)[[3]]
  coef.names <- "hamming"  # This might be modified later
  if (is.null(a$cov)) {
    if(drop){ #   Check for zero statistics (Should this happen when !is.null(a$cov)?)
      obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
      if (extremewarnings(obsstats))
        return(NULL) # In this case the observed Hamming distance is zero.
    }
#    name <- "hamhamming_weighted"
    if (length(sc03)>1) 
      coef.names <- paste("hamming", as.character(sc03[[2]]), sep=".")
    covm <- NULL
    if (is.null(a$defaultweight))
      a$defaultweight <- 1.0
    emptynwstats <- NROW(xm) * a$defaultweight

  ## Process case with dyadcov (i.e. weighted) ##
  } else {
    # Extract dyadic covariate
    if(is.network(a$cov)){
      covm<-as.matrix.network(a$cov,matrix.type="edgelist",a$attrname)
      if(length(covm)==2){covm <- matrix(covm,ncol=2)}
      if(length(covm)==3){covm <- matrix(covm,ncol=3)}
      if (NCOL(covm)==2)
        covm <- cbind(covm,1)
    }else if(is.character(a$cov)){
      covm<-get.network.attribute(nw,a$cov)
      covm<-as.matrix.network(covm,matrix.type="edgelist") # DH:  Not really sure what should happen here
    }else{
      covm<-as.matrix(a$cov)
    }
    if (is.null(covm) || !is.matrix(covm) || NCOL(covm)!=3){
      stop("Improper dyadic covariate passed to hamming()", call.=FALSE)
    }
    #    name = "hamhamming_weighted"
    emptynwstats <- sum(apply(xm, 1, function(a,b) sum(b[(a[1]==b[,1] & a[2]==b[,2]),3]), covm))
    if (is.null(a$defaultweight))
      a$defaultweight <- 0
    if(!is.null(a$attrname) && length(sc03)>1){
      coef.names<-paste("hamming", as.character(sc03[2]), "wt",
                as.character(a$attrname), sep = ".")
    }else if (length(sc03)>1) {
      coef.names<-paste("hamming", as.character(sc03[2]), "wt",
                as.character(sys.call(0)[[3]][3]), sep = ".")
    }
  }
  ## Return ##
  if (!is.null(xm)) {
    if (!is.directed(nw)) {
      tmp <- apply(xm, 1, function(a) a[1]>a[2])
      xm[tmp,] <- xm[tmp,2:1]
    }
    xm <- xm[order(xm[,1], xm[,2]), , drop=FALSE]
  }
  if (!is.null(covm)) {
    if (!is.directed(nw)) {
      tmp <- apply(covm, 1, function(a) a[1]>a[2])
      covm[tmp,] <- covm[tmp,c(2,1,3)]
    }
   covm <- covm[order(covm[,1], covm[,2]), , drop=FALSE]
  }
  inputs <- c(NROW(xm), as.vector(xm), a$defaultweight, NROW(covm), as.vector(covm))
  list(name="hamhamming", coef.names=coef.names, #name and coef.names: required 
       inputs = inputs, emptynwstats = emptynwstats, dependence = FALSE)
}

############################## InitErgmTerm functions:  I
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
  ### Construct the list to return
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
  ### Construct the list to return
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

  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      print("Warning:  nodefactor term deleted because it contributes no statistics")
      return()
    }
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
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodefactor",                                        #required
       coef.names = paste("nodefactor", paste(a$attrname,collapse="."), u, sep="."), #required
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

  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      print("Warning:  nodeifactor term deleted because it contributes no statistics")
      return()
    }
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
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeifactor",                                        #required
       coef.names = paste("nodeifactor", paste(a$attrname,collapse="."), u, sep="."), #required
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
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
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
  ### Construct the list to return
  if (a$diff) {
    coef.names <- paste("nodematch", paste(a$attrname,collapse="."), u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", paste(a$attrname,collapse="."), sep=".")
    inputs <- nodecov
  }
  list(name="nodematch",                                 #name: required
       coef.names = coef.names,                          #coef.names: required
       inputs =  inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}

#########################################################
InitErgmTerm.nodemix<-function (nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname", "base"),
                      vartypes = c("character", "numeric"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  ### Process the arguments
  if (is.bipartite(nw) && is.directed(nw)) {
    stop("Directed bipartite networks are not currently possible")
  }
  
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
    
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # deal with ew variable later
  } else { ew <- FALSE }
  if (is.bipartite(nw)) {
    #  So undirected network storage but directed mixing
    nb1 <- get.network.attribute(nw, "bipartite")       
    #  Recode nodecov to numeric (but retain original sorted names in "namescov")
    b1namescov <- sort(unique(nodecov[1:nb1]))
    b2namescov <- sort(unique(nodecov[(1+nb1):network.size(nw)]))
    namescov <- c(b1namescov, b2namescov)
    b1nodecov <- match(nodecov[1:nb1],b1namescov)
    b2nodecov <- match(nodecov[(1+nb1):network.size(nw)],b2namescov)
    nr <- length(b1namescov)
    nc <- length(b2namescov)
    nodecov <- c(b1nodecov, b2nodecov + nr)
    u <- cbind(rep(1:nr,nc), nr + rep(1:nc, each=nr))
    if(any(is.na(nodecov))){u<-rbind(u,NA)}    
    if (!is.null(a$base) && !identical(a$base,0)) {
      u <- u[-a$base,]
    }
    u <- u[!ew,]
    name <- "mix"
    cn <- paste("mix", paste(a$attrname,collapse="."), apply(matrix(namescov[u],ncol=2),
                                       1,paste,collapse="."), sep=".")
    inputs <- c(u[,1], u[,2], nodecov)
    attr(inputs, "ParamsBeforeCov") <- NROW(u)
  } else {# So one mode, but could be directed or undirected
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
    if (!is.null(a$base) && !identical(a$base,0)) {
      urm <- as.vector(urm)[-a$base]
      ucm <- as.vector(ucm)[-a$base]
      uun <- as.vector(uun)[-a$base]
    }
    urm <- urm[!ew]
    ucm <- ucm[!ew]
    uun <- uun[!ew]
    name <- "nodemix"
    cn <- paste("mix", paste(a$attrname,collapse="."), uun, sep=".")
    inputs <- c(urm, ucm, nodecov)
    #attr(inputs, "ParamsBeforeCov") <- 2*length(uui)
    attr(inputs, "ParamsBeforeCov") <- 2*length(uun)
  }
  ### Construct the list to return
  list(name = name, coef.names = cn, # required
       inputs = inputs, 
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
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }

  u <- sort(unique(nodecov))
  if (!is.null(a$base) && !identical(a$base,0)) {
    u <- u[-a$base]
    if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
      print("Warning:  nodeofactor term deleted because it contributes no statistics")
      return()
    }
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
  ### Construct the list to return
  inputs <- c(ui, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(ui) # See comment at top of file
  list(name="nodeofactor",                                        #required
       coef.names = paste("nodeofactor", paste(a$attrname,collapse="."), u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  

############################## InitErgmTerm functions:  T
#########################################################
InitErgmTerm.threepath <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, 
                       varnames = c("keep"),
                       vartypes = c("numeric"),
                       defaultvalues = list(1:4),
                       required = c(FALSE))
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } 
  else 
    ew <- FALSE
  types <- c("RRR","RRL","LRR","LRL")[a$keep]
  types <- types[!ew]  
  if (is.directed(nw)) {
    return(list(name = "threepath", 
                coef.names = paste("threepath", types, sep="."),
                inputs=a$keep[!ew]))
  }
  else {
    return(list(name = "threepath", coef.names = "threepath"))
  }
}

