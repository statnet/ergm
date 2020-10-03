#  File tests/constraint_conflict.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
library(ergm)
data(florentine)

efit <- ergm(flomarriage~edges, constraints=~edges)

stopifnot("The specified model's sample space constraint holds statistic(s) edges  constant. They will be ignored." %in% names(warnings()))
