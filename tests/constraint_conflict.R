#  File ergm/tests/constraint_conflict.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
library(ergm)
data(florentine)

#if(!inherits(try(
efit <- ergm(flomarriage~edges, constraints=~edges)#),
#   "try-error"))
#  stop("Should have had an error here.")
