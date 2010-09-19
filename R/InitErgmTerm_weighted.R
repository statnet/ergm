InitErgmTerm.sum<-function(nw, arglist, drop=TRUE, response=NULL, ...) {
  if(is.null(response)) stop('Statistic "sum" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="sum",
       coef.names="sum",
       inputs=NULL,
       dependence=FALSE)
}

InitErgmTerm.nonzero<-function(nw, arglist, drop=TRUE, response=NULL, ...) {
  if(is.null(response)) stop('Statistic "nonzero" requires a weighted network with a response= argument.')
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  list(name="nonzero",
       coef.names="nonzero",
       inputs=NULL,
       dependence=FALSE)
}
