#  File R/InitErgmTerm.blockop.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
## Creates a submodel that ignores any edges not within the
## blocks.

InitErgmTerm.NodematchFilter <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrname"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  attrname <- a$attrname
  c(list(name="on_blockdiag_net", submodel=m, auxiliaries=trim_env(~.blockdiag.net(attrname), "attrname")),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('NodematchFilter',a$attrname)))
}

