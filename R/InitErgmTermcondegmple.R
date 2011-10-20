InitErgmTerm.conddegmple<-function (nw, arglist, ...) {
  ### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("type"),
                      vartypes = c("character"),
                      defaultvalues = list("tetrad"),
                      required = c(FALSE))
  ### Process the arguments
  type<-a$type
  ### Construct the list to return
  out=list(name="conddegmple",
       coef.names="conddegmple",
       inputs=NULL,
       emptynwstats=NULL,
       dependence=TRUE)
  out
}
