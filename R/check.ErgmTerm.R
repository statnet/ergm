check.ErgmTerm <- function(nw, arglist, directed=NULL, bipartite=NULL,
                           varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL) {
  ## This function (called by InitErgmTerm.xxx functions) will
  ## (a) check to make sure that the network has the correct "directed"
  ##     and "bipartite" attributes, if applicable, i.e., if either directed
  ##     of bipartite is non-NULL.
  ## (b) check that the inputs to the model term are of the correct type
  ## (c) assign default values if applicable
  ## (d) return a list with elements (var1=var1value, var2=var2value, etc.)
  ## Note that the list returned will contain the maximum possible number of
  ## arguments; any arguments without values are returned as NULL.
  ##
  ##   Inputs:  nw (required) is the network
  ##            arglist (required) is the list of arguments passed from ergm.getmodel
  ##            directed is logical if directed=T or F is required, NULL o/w
  ##            bipartite is logical if bipartite=T or F is required, NULL o/w
  ##            varnames is a vector of variable names
  ##            vartypes is a vector of corresponding variable types
  ##            defaultvalues is a list of default values (NULL means no default)
  ##            required is a vector of logicals:  Is this var required or not?

  fname <- get.InitErgm.fname() # From what InitErgm function was this called?
  fname <- sub('.*[.]', '', fname) # truncate up to last '.'
  message <- NULL
  if (!is.null(directed) && directed != (dnw<-is.directed(nw))) {
    #directed != (dnw<-eval(expression(nw$gal$dir),parent.frame()))) {
    message <- paste("networks with directed==", dnw, sep="")
  }
  if (!is.null(bipartite) && bipartite != ((bnw <- nw %n% "bipartite") > 0)) {
    #bipartite != (bnw <- eval(expression(nw %n% "bipartite"),parent.frame()) > 0)) {
    message <- paste("networks with bipartite", 
                     ifelse(bnw>0, " > 0", "==FALSE"), sep="")
  }
  if (!is.null(message)) {
    stop(paste("The ERGM term",fname,"may not be used with",message))
  }

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
# that each InitErgmTerm function faithfully passes in what the user typed;
# thus, the correctness of input from the InitErgmTerm function isn't checked.
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

check.ErgmTerm.summarystats <- function(nw, args, ...) {
  fname <- get.InitErgm.fname() # From what InitErgm function was this called?
  Initfn <- get(fname,envir=.GlobalEnv)
  outlist <- Initfn(nw, args, drop=FALSE, ...)
  m <- updatemodel.ErgmTerm(list(), outlist)
  gs <- ergm.getglobalstats(nw, m)
}

# Search back in time through sys.calls() to find the name of the last
# function whose name begins with "InitErgm"
get.InitErgm.fname <- function() {
  sc <- sys.calls()
  i <- length(sc)
  listofnames <- NULL
  while (i>1) { 
    i <- i-1
    fname <- as.character(sc[[i]][1])
    listofnames <- c(listofnames, fname)
    if (substring(fname,1,8)=="InitErgm") {
      return(fname)
    }
  }
  # Didn't find InitErgm... in the list of functions
  return(NULL)
}

# zerowarnings is now deprecated; it's been replaced by extremewarnings
zerowarnings <- function(gs) {
  out <- gs==0
  if(any(out)) {
    cat(" Warning:  These coefs set to -Inf due to obs val=0:\n ")
    cat(paste(names(gs)[out], collapse=", "), "\n")
  }
  out
}

extremewarnings <- function(gs, minval=0, maxval=NULL) {
  out <- rep(FALSE, times=length(gs))
  if (!is.null(minval))  {
    tmp <- (gs <= minval)
    if (any(tmp)) {
      out <- out | tmp
      cat(" Warning:  These coefficients will be -Inf:\n")
      cat(paste(names(gs)[tmp], collapse=", "), "\n")
    }
  }
  if (!is.null(minval))  {
    tmp <- (gs >= maxval)
    if (any(tmp)) {
      out <- out | tmp
      cat(" Warning:  These coefficients will be +Inf:\n")
      cat(paste(names(gs)[tmp], collapse=", "), "\n")
    }
  }
  out
}


