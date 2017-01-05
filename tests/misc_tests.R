#  File tests/misc_tests.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################

# tests for miscelaneous things like specific bug fixes


#I just found a bug where the the deviance was incorrectly
#computed/reported for all of dyad-independence models. 
#  bug  https://statnet.csde.washington.edu/trac/ticket/1182
library(ergm)
data(florentine)
# If FALSE then the bug exists.
if(ergm(flomarriage ~ edges)$mle.lik < -60){
  stop("deviance was incorrectly computed/reported for all of dyad-independence models")
}