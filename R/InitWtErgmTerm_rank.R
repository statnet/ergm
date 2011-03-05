InitWtErgmTerm.deference<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)

  list(name="d_deference",
       coef.names="deference",
       inputs=NULL,
       dependence=TRUE)
}

# nodeicov for ranks is implemented in InitWtErgmTerm.R

InitWtErgmTerm.nonconformity<-function(nw, arglist, response, drop=TRUE, ...) {
  a <- check.ErgmTerm(nw, arglist, directed=TRUE,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)

  list(name="d_nonconformity",
       coef.names="nonconformity",
       inputs=NULL,
       dependence=TRUE)
}
