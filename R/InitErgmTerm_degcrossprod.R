#===========================================================================
# This file contains the following 4 new, easier-to-write degree-based
# ergm-term  initialization functions (each prepended with "InitErgmTerm"):
#          <degcrossprod>      <adegcor>
#          <degcor>            <rdegcor>
#====================================================================




##############################################################################
# The <InitErgmTerm.X> functions initialize each ergm term, X, by
#   1) checking the validity of X and its arguments via <check.ErgmTerm> and
#   2) setting appropiate values for each of the components in the returned list
# X is initialized for inclusion into a model that is specified by formula F and
# built via <ergm.getmodel>
# 
# --PARAMETERS--
#   nw        : the network given in formula F
#   arglist   : the arguments given with term X in formula F
#   drop      : whether coefficients with 0 values should be
#               dropped; default=TRUE
#
# --IGNORED PARAMETERS--
#   ... : ignored, but necessary to accomodate other arguments
#         passed by <ergm.getmodel>
#
# --RETURNED--
#   a list of term-specific elements required by the C changestats
#   functions and other R rountines; the first two components of this
#   list are required*, the remaining components are optional:
#     *name      : the name of term X; this is used to locate the C function
#                  calculating the change statistics for X, which will be
#                  'name' prepended with "d_"; for example if X=absdiff,
#                  'name'="absdiff", and the C function is "d_absdiff"
#     *coef.names: the vector of names for the coefficients (parameters)
#                  as they will be reported in the output
#     inputs     : the vector of (double-precision numeric) inputs that the 
#                  changestat function called d_<name> will require
#                  (see WHAT THE C CHANGESTAT FUNCTION RECEIVES below);
#                  this MUST be a vector!; thus, if the inputs are  matrices,
#                  they must be "flattened" to vectors; if they are categorical
#                  character-valued variables, they must be converted to numbers;
#                  optionally, 'inputs' may have an attribute named
#                  "ParamsBeforeCov",which is the number that used to be the
#                  old Element 1 and is needed for backwards compatability
#                  (see the old <InitErgm> for details); default=NULL
#     dependence : whether the addition of term X to the model makes the model
#                  into a dyadic dependence model (T or F); if all terms have
#                  'dependence' set FALSE, the model is assumed to be a
#                  dyadic independence model; default=TRUE
#
##############################################################################################


InitErgmTerm.degcrossprod<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 
  ### Construct the list to return
  list(name="degcrossprod",                            #name: required
       coef.names = "degcrossprod",                    #coef.names: required
       inputs=2*summary(nw ~ edges),
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}


InitErgmTerm.degcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="degcor",                            #name: required
       coef.names = "degcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}


InitErgmTerm.adegcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="adegcor",                            #name: required
       coef.names = "adegcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}


InitErgmTerm.rdegcor<-function (nw, arglist, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=FALSE) 

  deg=summary(nw ~ sociality(base=0))
  el=as.matrix.network.edgelist(nw)
  deg1<-deg[el[,1]]
  deg2<-deg[el[,2]]
  alldeg<-c(deg1,deg2)
  sigma2<-(sum(alldeg*alldeg)-length(alldeg)*(mean(alldeg)^2))
  ### Construct the list to return
  list(name="rdegcor",                            #name: required
       coef.names = "rdegcor",                    #coef.names: required
       inputs=sigma2,
       dependence = TRUE # So we don't use MCMC if not necessary
       )
}
