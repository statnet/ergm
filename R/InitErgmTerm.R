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
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  ### Construct the list to return
  list(name="absdiff",                                     #name: required
       coef.names = paste("absdiff", attrname, sep="."), #coef.names: required
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  u <- sort(unique(as.vector(abs(outer(nodecov,nodecov,"-")))),na.last=NA)
  u <- u[u>0]
  NAsubstitute <- 2*(1+max(abs(c(nodecov,u)),na.rm=TRUE)) # Arbitrary unused (and nonzero) value
  napositions <- is.na(nodecov)
  nodecov[napositions] <- NAsubstitute
  if(any(napositions)){u<-c(u,NA)}
  if(!is.null(base)) u <- u[-(base)]
  if (length(u)==0)
    stop ("Argument to absdiffcat() has too few distinct differences", call.=FALSE)
  u2 <- u[!is.na(u)]
  ### Construct the list to return
  inputs <- c(u2, NAsubstitute, nodecov)
  attr(inputs, "ParamsBeforeCov") <- length(u2)+1 # See comment at top of file
  list(name="absdiffcat",                                  #name: required
       coef.names = paste("absdiff", attrname, u, sep="."), #coef.names: required
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if(!initialfit && !fixed){ # This is a curved exponential family model
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
       params=list(altkstar=NULL, altkstar.lambda=lambda)
       )
  } else {
    outlist <- list (name="altkstar",                      #name: required
       coef.names = paste("altkstar", lambda, sep="."),  #coef.names: required
       inputs=lambda
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if (!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname)
    u <- sort(unique(nodecov))
    if (!is.null(keep)) {
      u <- u[keep]
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
    if (is.null(attrname)) {
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
  if (!is.null(attrname)) {
    if (diff) {
      out$coef.names <- paste("asymmetric", attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("asymmetric", attrname, sep=".")
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
                       varnames = c("d", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(attrname)) {  # CASE 1:  attrname GIVEN
    nodecov <- get.node.attr(nw, attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    du <- matrix(du[,!ew], nrow=2) # Drop any zero-obs-value rows
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov[1:nb1]==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b1degree_by_attr"
    coef.names <- paste("b1deg", du[1,], ".", attrname, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  attrname NOT GIVEN
    d <- d[!ew] # Drop any zero-obs-value values
    name <- "b1degree"
    coef.names <- paste("b1deg", d, sep="")
    inputs <- d
    emptynwstats <- rep(0, length(d))
    if (any(d==0)) { # alter emptynwstats
      emptynwstats[d==0] <- nb1
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    # Recode to numeric
    nodecov <- match(nodecov,u,nomatch=length(u)+1)
    name <- "ostar"
    coef.names <- paste("b1star", k, ".", attrname, sep="")
    inputs <- c(k, nodecov)
    attr(inputs, "ParamsBeforeCov") <- length(k)
  } 
  else {
    name <- "ostar"
    coef.names <- paste("b1star",k,sep="")
    inputs <- k
  }
  list(name = name, coef.names = coef.names, #name and coef.names: required
       inputs = inputs)
}

#########################################################
InitErgmTerm.b2degree <- function(nw, arglist, drop=TRUE, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d", "attrname"),
                       vartypes = c("numeric", "character"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, FALSE))
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  n <- network.size(nw)
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    # will process the ew variable later
  } else {ew <- FALSE}
  if (!is.null(attrname)) {  # CASE 1:  attrname GIVEN
    nodecov <- get.node.attr(nw, attrname)
    u<-sort(unique(nodecov))
    if(any(is.na(nodecov))){u<-c(u,NA)}
    nodecov <- match(nodecov,u) # Recode to numeric
    # Combine degree and u into 2xk matrix, where k=length(d)*length(u)
    lu <- length(u)
    du <- rbind(rep(d,lu), rep(1:lu, rep(length(d), lu)))
    du <- matrix(du[,!ew], nrow=2) # Drop any zero-obs-value rows
    emptynwstats <- rep(0, ncol(du))
    if (any(du[1,]==0)) { # Alter emptynwstats
      tmp <- du[2,du[1,]==0]
      for(i in 1:length(tmp)) 
        tmp[i] <- sum(nodecov[(1+nb1):n]==tmp[i])
      emptynwstats[du[1,]==0] <- tmp
    }
    name <- "b2degree_by_attr"
    coef.names <- paste("b2deg", du[1,], ".", attrname, u[du[2,]], sep="")
    inputs <- c(as.vector(du), nodecov)
  } else { # CASE 2:  attrname NOT GIVEN
    d <- d[!ew] # Drop any zero-obs-value values
    name <- "b2degree"
    coef.names <- paste("b2deg", d, sep="")
    inputs <- d
    emptynwstats <- rep(0, length(d))
    if (any(d==0)) { # alter emptynwstats
      emptynwstats[d==0] <- n-nb1
    }
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats)
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if(drop) { # Check for zero statistics, print -Inf messages if applicable
    obsstats <- check.ErgmTerm.summarystats(nw, arglist, ...)
    ew <- extremewarnings(obsstats)
    k <- k[!ew]
  }
  if (length(k)==0) return(NULL)
  ### Construct the list to return
  list(name="cycle",                            #name: required
       coef.names = paste("cycle", k, sep=""),  #coef.names: required
       inputs = c(max(k), (2:max(k)) %in% k)
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if(is.network(x))
    xm<-as.matrix.network(x,matrix.type="adjacency",attrname)
  else if(is.character(x))
    xm<-get.network.attribute(nw,x)
  else
    xm<-as.matrix(x)
  ### Construct the list to return
  if(!is.null(attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov", as.character(sys.call(0)[[3]][2]), 
              as.character(attrname), sep = ".")
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if (!is.null(attrname)) {
    nodecov <- get.node.attr(nw, attrname)
    u <- sort(unique(nodecov))
    if (!is.null(keep)) {
      u <- u[keep]
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
    if (is.null(attrname)) {
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
  if (!is.null(attrname)) {
    if (diff) {
      out$coef.names <- paste("mutual", attrname, u, sep=".")
      out$inputs <- c(ui, nodecov)
    } else {
      out$coef.names <- paste("mutual", attrname, sep=".")
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  u <- sort(unique(nodecov))
  if (!is.null(base) && !identical(base,0)) {
    u <- u[-base]
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
       coef.names = paste("nodefactor", attrname, u, sep="."), #required
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  u <- sort(unique(nodecov))
  if (!is.null(base) && !identical(base,0)) {
    u <- u[-base]
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
       coef.names = paste("nodeifactor", attrname, u, sep="."), #required
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  u <- sort(unique(nodecov))
  if (!is.null(keep)) {
    u <- u[keep]
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
  if (diff) {
    coef.names <- paste("nodematch", attrname, u, sep=".")
    inputs <- c(ui, nodecov)
  } else {
    coef.names <- paste("nodematch", attrname, sep=".")
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  if (is.bipartite(nw) && is.directed(nw)) {
    stop("Directed bipartite networks are not currently possible")
  }
  nodecov <- get.node.attr(nw, attrname)
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
    if (!is.null(base) && !identical(base,0)) {
      u <- u[-base,]
    }
    u <- u[!ew,]
    name <- "mix"
    cn <- paste("mix", attrname, apply(matrix(namescov[u],ncol=2),
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
    if (!is.null(base) && !identical(base,0)) {
      urm <- as.vector(urm)[-base]
      ucm <- as.vector(ucm)[-base]
      uun <- as.vector(uun)[-base]
    }
    urm <- urm[!ew]
    ucm <- ucm[!ew]
    uun <- uun[!ew]
    name <- "nodemix"
    cn <- paste("mix", attrname, uun, sep=".")
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
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nodecov <- get.node.attr(nw, attrname)
  u <- sort(unique(nodecov))
  if (!is.null(base) && !identical(base,0)) {
    u <- u[-base]
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
       coef.names = paste("nodeofactor", attrname, u, sep="."), #required
       inputs = inputs,
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}  


