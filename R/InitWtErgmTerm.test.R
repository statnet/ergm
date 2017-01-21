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

InitWtErgmTerm..sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("mode"),
                      vartypes = c("character"),
                      defaultvalues = list("numeric"),
                      required = c(FALSE))

  mode <- match.arg(a$mode, c("numeric"))
  name <- switch(mode,
                 numeric = "_dsociomatrix")
  
  list(name=name,
       coef.names=c(), dependence=FALSE)
}

InitWtErgmTerm.sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("mode"),
                      vartypes = c("character"),
                      defaultvalues = list("numeric"),
                      required = c(FALSE))

  mode <- match.arg(a$mode, c("numeric"))
  name <- switch(mode,
                 numeric = "dsociomatrix")
  n <- network.size(nw)
  tails <- rep(1:n,n)
  heads <- rep(1:n,each=n)
  list(name=name,
       coef.names=paste(tails,heads,sep="."), dependence=FALSE,
       auxiliaries = ~.sociomatrix(mode))
}

