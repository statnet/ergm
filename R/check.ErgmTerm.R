#====================================================================================
# This file contains the following 6 files that help check the validity of ergm terms
#       <check.ErgmTerm>               <get.InitErgm.fname>
#       <assignvariables>              <zerowarnings>
#       <check.ErgmTerm.summarystats>  <extremewarnings>
#====================================================================================




######################################################################################
# The <check.ErgmTerm> function ensures for the <InitErgmTerm.X> function that the
# term X:
#   1) is applicable given the 'directed' and 'bipartite' attributes of the given
#      network
#   2) has an appropiate number of arguments
#   3) has correct argument types if arguments where provided
#   4) has default values assigned if defaults are available
# by halting execution if any of the first 3 criteria are not met
#
# --PARAMETERS--
#  nw           : the network that term X is being checked against  
#  arglist      : the list of arguments for term X
#  directed     : whether term X requires a directed network (T or F); default=NULL
#  bipartite    : whether term X requires a bipartite network (T or F); default=NULL
#  varnames     : the vector of names of the possible arguments for term X;
#                 default=NULL 
#  vartypes     : the vector of types of the possible arguments for term X;
#                 default=NULL 
#  defaultvalues: the list of default values for the possible arguments of term X;
#                 default=list()
#  required     : the logical vector of whether each possible argument is required;
#                 default=NULL
#
# --RETURNED--
#   out: a list of the values for each possible argument of term X; user provided 
#        values are used when given, default values otherwise.
#
######################################################################################

check.ErgmTerm <- function(nw, arglist, directed=NULL, bipartite=NULL,
                           varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL) {
  fname <- get.InitErgm.fname() # From what InitErgm function was this called?
  fname <- sub('.*[.]', '', fname) # truncate up to last '.'
  message <- NULL
  if (!is.null(directed) && directed != (dnw<-is.directed(nw))) {
    #directed != (dnw<-eval(expression(nw$gal$dir),parent.frame()))) {
    message <- paste("networks with directed==", dnw, sep="")
  }
  if(is.null(bnw<- nw %n% "bipartite")) bnw <- 0
  if (!is.null(bipartite) && bipartite != (bnw > 0)) {
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
  #  c(.conflicts.OK=TRUE,out)
  out
}




############################################################################
## As of ergm version 2.2, the assignvariables function is deprecated.
## This is so that there are no "mysterious" variable assignments in the
## InitErgmTerm functions, which is prehaps better programming style
## and which also prevents a whole raft of warnings when using
## R CMD build.
## However, assignvariables will still function so as not to break 
## code.  It will simply produce a warning when called.
##############################################################################
# The <assignvariables> function following creates a variable out of each of
# the elements in the given list *in the calling environment*.  In this way, it
# sort of works like "attach" but without creating a new environment and
# without all of the headaches that "attach" can give because the variables
# it creates are not in the correct frame.
##############################################################################

assignvariables <- function(a) {
  cat("The assignvariables function has been deprecated.  Please modify\n",
       "the", get.InitErgm.fname(), "function so that it does not rely on\n",
       "this function.  For instance, instead of 'assignvariables(a)'\n",
       "followed by using 'attrname' throughout, use 'a$attrname' throughout.\n ")
  if(length(a)>0)
    for(i in 1:length(a)) 
      assign(names(a)[i], a[[i]], envir=parent.frame())
}


check.ErgmTerm.summarystats <- function(nw, arglist, ...) {
  fname <- get.InitErgm.fname() # From what InitErgm function was this called?
  Initfn <- get(fname,envir=.GlobalEnv)
  outlist <- Initfn(nw, arglist, drop=FALSE, ...)
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



####################################################################
# The <extremewarnings> function checks and returns whether the
# global statistics are extreme, in terms of being outside the
# range of 'minval' and 'maxval'; warning messages are printed if
# any statistic is extreme
#
# --PARAMETERS--
#   gs    : the vector of global statistics returned by
#           <check.ErgmTerm.summarystats> or <ergm.getglobalstats>
#   minval: the value at/below which statistics are deemed extreme;
#           default=0
#   maxval: the value at/above which statistics are deemed extreme;
#           default=NULL
#
# --RETURNED--
#   out: a logical vector of whether each statistic was extreme
#
####################################################################

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


