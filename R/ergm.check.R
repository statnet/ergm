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

ergm.checkbipartite <- function(fname, nw.bipartiteflag, requirement,
                               extramessage="") {
  if (!nw.bipartiteflag && requirement)
    stop(paste(fname, "model term may not be used with an non-bipartite network.",
               extramessage), call.=FALSE)
  if (nw.bipartiteflag && !requirement)
    stop(paste(fname, "model term may not be used with a bipartite network.",
               extramessage), call.=FALSE)
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


ergm.checkdegeneracy <- function(statsmatrix, statsmatrix.miss=NULL, verbose=FALSE) {
#
 degen <- FALSE
 novar <- apply(statsmatrix,2,var)<1e-6
 if(all(novar)){
  if(verbose){
    warning("All the MCMC sample statistics are the same.\n", call.=FALSE)
    print(apply(statsmatrix,2,summary.statsmatrix.ergm),scipen=6)
  }
  degen <- TRUE
 }
# 
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
#  novar=rep(FALSE,length(theta0)) # This is a hack
#  theta0[!novar] <- l$coef
#  l$coef <- theta0
#  theta0[!novar] <- l$MCMCtheta
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
