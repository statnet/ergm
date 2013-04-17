#  File R/ergm.check.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
#============================================================================
# This file contains the following 4 functions for checking ergm arguments
#          <ergm.checkargs>           <ergm.checkdirected>
#          <ergm.checkbipartite>      <ergm.checkdegeneracy>
#============================================================================





######################################################################################
# The <ergm.checkargs> function ensures for the <InitErgm.X> function that the
# term X:
#   1) has an appropiate number of arguments
#   2) has correct argument types if arguments where provieded
#   3) has default values assigned for non-required arguments
# by halting execution if either of the first 2 criteria are not met
#
# --PARAMETERS--
#  fname        : the name of the model term as a character string
#  arglist      : the list of arguments for term X
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
#   a vector containing the following 2 variables:
#     .conflicts.OK: always TRUE
#     out          : a list of the values for each possible argument of term X;
#                    user provided values are used when given, default values otherwise.
#
######################################################################################

ergm.checkargs <- function(fname, arglist, varnames=NULL, vartypes=NULL,
                           defaultvalues=list(), required=NULL) {
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



#################################################################################
# The <ergm.checkdirected> function halts execution for the <InitErgm> functions
# with an error message if the given model term cannot be used with the network
# because of its bipartite state
#
# --PARAMETERS--
#   fname            : the name of the model term as a character string
#   nw.bipartiteflag : whether the network is bipartite (T or F)
#   requirement      : whether the term requires a bipartite network (T or F)
#   extramessage     : additional messages to attach to the warning;
#                      default value = ""
#
# --RETURNED--
#   NULL
#
#################################################################################

ergm.checkbipartite <- function(fname, nw.bipartiteflag, requirement,
                               extramessage="") {
  if (!nw.bipartiteflag && requirement)
    stop(paste(fname, "model term may not be used with an non-bipartite network.",
               extramessage), call.=FALSE)
  if (nw.bipartiteflag && !requirement)
    stop(paste(fname, "model term may not be used with a bipartite network.",
               extramessage), call.=FALSE)
}




#################################################################################
# The <ergm.checkdirected> function halts execution for the <InitErgm> functions
# with an error message if the given model term cannot be used with the network
# because of its state as (un)directed
#
# --PARAMETERS--
#   fname           : the name of the model term as a character string
#   nw.directedflag : whether the network is directed (T or F)
#   requirement     : whether the term requires a directed network (T or F)
#   extramessage    : additional messages to attach to the warning;
#                     default value = ""
#
# --RETURNED--
#   NULL
#
#################################################################################

ergm.checkdirected <- function(fname, nw.directedflag, requirement,
                               extramessage="") {
  if (!nw.directedflag && requirement)
    stop(paste(fname, "model term may not be used with an undirected network.",
               extramessage), call.=FALSE)
  if (nw.directedflag && !requirement)
    stop(paste(fname, "model term may not be used with a directed network.",
               extramessage), call.=FALSE)
}





#################################################################################
# The <ergm.checkdegeneracy> function checks for degeneracy by looking for
# variability in the stats matrix 
#
# --PARAMETERS--
#   statsmatrix :  the matrix of summary sample statistics
#   verbose     :  whether the degeneracy warning should be printed (T or F)
#
#
# --IGNORED PARAMETERS--
#   statsmatrix.obs : the matrix of summary sample statistics from the
#                      observation process network; default=NULL
#
#
# --RETURNED--
#   degen: whether the ergm model is degenerate, in the sense that
#          there was no varibility in the statsmatrix (T or F)
#
#################################################################################

ergm.checkdegeneracy <- function(statsmatrix, statsmatrix.obs=NULL, verbose=FALSE) {
 degen <- FALSE
 novar <- apply(statsmatrix,2,var)<1e-6
 if(all(novar)){
  if(verbose){
    warning("All the MCMC sample statistics are the same.\n", call.=FALSE)
    print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
  }
  degen <- TRUE
 }
#  ########### from ergm.estimate:
#   # First, perform simple check to see if observed statistics
#   # (all zeros) lie outside the range of the MCMC sample statistics.
#  outofbounds <- (apply(statsmatrix,2,max)*
#                  apply(statsmatrix,2,min)) > epsilon
#  if (any(outofbounds)){
#    ergm.marquardt()
#  }
#
#  ############# from ergm.statseval:  CHECK FOR NO VARIANCE
#  # below is what was done; not necessarily something to keep
#  novar=rep(FALSE,length(init)) # This is a hack
#  init[!novar] <- l$coef
#  l$coef <- init
#  init[!novar] <- l$MCMCtheta
#  covar <- 0*diag(length(l$coef)) # initialize to zero matrix
#  covar[!novar,!novar] <- l$covar
#
#
#  ############## from ergm.mainfitloop:
#  if(sum(z$newg) > 12000){
#    cat("Returned network is too full. Retaining the original.\n")
#    newnetwork <- g
#  }
degen
}
