InitWtErgmTerm.test.abs.sum.minus.5<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = "summary",
                      vartypes = "logical",
                      defaultvalues = list(TRUE),
                      required = FALSE)
  
  list(name=if(a$summary) "test_abs_sum_minus_5" else "test_abs_sum_minus_5_no_s",
       coef.names=if(a$summary) "test_abs_sum_minus_5" else "test_abs_sum_minus_5_no_summary", dependence=TRUE, emptynwstats = 5,
       minval = 0, maxval = +Inf, conflicts.constraints="sum")
}

