#  File R/InitErgmTerm.indices.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
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
