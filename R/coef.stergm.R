#  File ergm/R/coef.stergm.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
coef.stergm <- function(object, ...){list(formation=object$formation.fit$coef,
                                          dissolution=object$dissolution.fit$coef)}
coefficients.stergm <- coef.stergm
