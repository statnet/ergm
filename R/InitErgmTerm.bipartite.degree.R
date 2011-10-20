

###########  InitErgmTerm.b1degree.edgecov ###################


InitErgmTerm.b1degree.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d","edgecov"),
                       vartypes = c("numeric", "network"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, TRUE))
  assignvariables(a) # create local variables with names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1

  d <- sort(unique(d))
  if(sum(dim(edgecov) != c(nb1,nb2))>0)
    stop("Covariate is not of same dimensions as network", call.=FALSE)

  name <- "b1degree_edgecov"
  coef.names <- paste("b1deg.edgecov",as.character(sys.call(0)[[3]][[3]]),d,sep=".")
  edgecov.vector <- as.vector(edgelist.ergm(edgecov))
  inputs <- c(d,edgecov.vector)
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb1
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}




###########  InitErgmTerm.b2degree.edgecov ###################


InitErgmTerm.b2degree.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d","edgecov"),
                       vartypes = c("numeric", "network"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, TRUE))
  assignvariables(a) # create local variables with names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1

  d <- sort(unique(d))
  if(sum(dim(edgecov) != c(nb1,nb2))>0)
    stop("Covariate is not of same dimensions as network", call.=FALSE)

  name <- "b2degree_edgecov"
  coef.names <- paste("b2deg.edgecov",as.character(sys.call(0)[[3]][[3]]),d,sep=".")
  edgecov.vector <- as.vector(edgelist.ergm(edgecov))
  inputs <- c(d,edgecov.vector)
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb2
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}


###########  InitErgmTerm.b1mindegree  ###################

InitErgmTerm.b1mindegree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d"),
                       vartypes = c("numeric"),
                       defaultvalues = list(NULL),
                       required = c(TRUE))
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1
  name <- "b1mindegree"
  coef.names <- paste("b1mindeg", d, sep="")
  inputs <- d
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb1
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}

###########  InitErgmTerm.b2mindegree  ###################

InitErgmTerm.b2mindegree <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d"),
                       vartypes = c("numeric"),
                       defaultvalues = list(NULL),
                       required = c(TRUE))
  assignvariables(a) # create local variables from names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1
  name <- "b2mindegree"
  coef.names <- paste("b2mindeg", d, sep="")
  inputs <- d
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb2
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0)
}

###########  InitErgmTerm.b1mindegree.edgecov ###################

InitErgmTerm.b1mindegree.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d","edgecov"),
                       vartypes = c("numeric", "network"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, TRUE))
  assignvariables(a) # create local variables with names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1

  d <- sort(unique(d))
  if(sum(dim(edgecov) != c(nb1,nb2))>0)
    stop("Covariate is not of same dimensions as network", call.=FALSE)

  name <- "b1mindegree_edgecov"
  coef.names <- paste("b1mindeg.edgecov",as.character(sys.call(0)[[3]][[3]]),d,sep=".")
  edgecov.vector <- as.vector(edgelist.ergm(edgecov))
  inputs <- c(d,edgecov.vector)
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb1
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}




###########  InitErgmTerm.b2mindegree.edgecov ###################


InitErgmTerm.b2mindegree.edgecov <- function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm (nw, arglist, directed=FALSE, bipartite=TRUE,
                       varnames = c("d","edgecov"),
                       vartypes = c("numeric", "network"),
                       defaultvalues = list(NULL, NULL),
                       required = c(TRUE, TRUE))
  assignvariables(a) # create local variables with names in 'varnames'
  ### Process the arguments
  nb1 <- get.network.attribute(nw, "bipartite")
  nb2 <- get.network.attribute(nw, "n") - nb1

  d <- sort(unique(d))
  if(sum(dim(edgecov) != c(nb1,nb2))>0)
    stop("Covariate is not of same dimensions as network", call.=FALSE)

  name <- "b2mindegree_edgecov"
  coef.names <- paste("b2mindeg.edgecov",as.character(sys.call(0)[[3]][[3]]),d,sep=".")
  edgecov.vector <- as.vector(edgelist.ergm(edgecov))
  inputs <- c(d,edgecov.vector)
  emptynwstats <- rep(0, length(d))
  if (any(d==0)) { # alter emptynwstats
    emptynwstats[d==0] <- nb2
  }
  list(name=name, coef.names=coef.names, #name and coef.names: required
       inputs = inputs, emptynwstats=emptynwstats,
       dependence = TRUE, minval = 0
       )
}


