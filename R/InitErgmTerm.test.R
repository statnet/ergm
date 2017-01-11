InitErgmTerm.test.abs.edges.minus.5<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = "summary",
                      vartypes = "logical",
                      defaultvalues = list(TRUE),
                      required = FALSE)
  
  list(name=if(a$summary) "test_abs_edges_minus_5" else "test_abs_edges_minus_5_no_s",
       coef.names=if(a$summary) "test_abs_edges_minus_5" else "test_abs_edges_minus_5_no_summary", dependence=TRUE, emptynwstats = 5,
       minval = 0, maxval = max(5,network.dyadcount(nw,FALSE)-5), conflicts.constraints="edges")
}

