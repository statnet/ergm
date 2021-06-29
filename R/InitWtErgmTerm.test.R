#  File R/InitWtErgmTerm.test.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2021 Statnet Commons
################################################################################
InitWtErgmTerm.test.abs.sum.minus.5<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("summary","aux"),
                      vartypes = c("logical","logical"),
                      defaultvalues = list(TRUE, FALSE),
                      required = c(FALSE,FALSE))
  
  list(name=if(a$aux) "test_abs_sum_minus_5_aux"
            else if(a$summary) "test_abs_sum_minus_5"
            else "test_abs_sum_minus_5_no_s",
       coef.names=if(a$aux) "test_abs_sum_minus_5_aux"
                  else if(a$summary) "test_abs_sum_minus_5"
                  else "test_abs_sum_minus_5_no_summary", dependence=TRUE, emptynwstats = 5,
       minval = 0, maxval = +Inf, conflicts.constraints="sum",
       auxiliaries=if(a$aux) trim_env(~.sum))
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
       auxiliaries = trim_env(~.sociomatrix(mode), "mode"))
}

InitWtErgmTerm..sum<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_sum",
       coef.names=c(), dependence=FALSE)
}
