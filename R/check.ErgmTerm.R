#  File R/check.ErgmTerm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
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
#   1.5) is not applied to a directed bipartite network
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
#  nonnegative  : whether term X requires a network with only nonnegative weights; default=FALSE
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

check.ErgmTerm <- function(nw, arglist, directed=NULL, bipartite=NULL, nonnegative=FALSE,
                           varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL, response=NULL) {
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
  if (is.directed(nw) && bnw > 0) {
    message <- "directed bipartite networks"
  }
  if (is.null(message) && nonnegative && any(nw %e% response < 0)){
    message <- "networks with negative dyad weights"
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
        if (all(sapply(strsplit(vartypes[m],",",fixed=TRUE)[[1]], function(vartype) !eval(call(paste("is.",vartype,sep=""),arglist[[i]]))))) {
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
        if (all(sapply(strsplit(vartypes[i],",",fixed=TRUE)[[1]], function(vartype) !eval(call(paste("is.",vartype,sep=""),arglist[[i]]))))) {
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
    if (substring(fname,1,8)=="InitErgm" || substring(fname,1,10)=="InitWtErgm") {
      return(fname)
    }
  }
  # Didn't find Init[Wt]Ergm... in the list of functions
  return(NULL)
}
