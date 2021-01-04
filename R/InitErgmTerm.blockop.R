## Creates a submodel that ignores any edges not within the
## blocks.

InitErgmTerm.NodematchFilter <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "attrname"),
                      vartypes = c("formula", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  m <- ergm_model(a$formula, nw,...)

  attrname <- a$attrname
  c(list(name="on_blockdiag_net", submodel=m, auxiliaries=trim_env(~.blockdiag.net(attrname), "attrname")),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap('NodematchFilter',a$attrname)))
}

