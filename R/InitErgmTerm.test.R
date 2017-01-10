InitErgmTerm.test.abs.edges.minus.5<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="test_abs_edges_minus_5", coef.names="test_abs_edges_minus_5", dependence=TRUE,
       minval = 0, maxval = max(5,network.dyadcount(nw,FALSE)-5), conflicts.constraints="edges")
}
