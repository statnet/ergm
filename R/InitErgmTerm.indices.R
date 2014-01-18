# This term is not meaningful for modelling, but it has the property
# that its change score for dyad (i,j) is always (i,j). It is used by
# ergmMPLE() to get a covariate matrix with each dyad identified.
InitErgmTerm.indices<-function(nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())
  ### Process the arguments
  list(name="indices",                                        #required
       coef.names = c("tail","head"), #required
       dependence = FALSE # So we don't use MCMC if not necessary
       )
}
