#  File R/InitErgmTerm.test.R in package ergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2017 Statnet Commons
#######################################################################
InitErgmTerm..sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("mode"),
                      vartypes = c("character"),
                      defaultvalues = list("integer"),
                      required = c(FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "_isociomatrix")
  
  list(name=name,
       coef.names=c(), dependence=FALSE)
}

InitErgmTerm..discord.sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "mode"),
                      vartypes = c("network", "character"),
                      defaultvalues = list(nw, "integer"),
                      required = c(FALSE, FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "_discord_isociomatrix")
  
  list(name=name,
       coef.names=c(), dependence=FALSE,
       auxiliaries = ~ .discord.net(a$x, implementation="Network"))
}

InitErgmTerm..discord.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "implementation"),
                      vartypes = c("network,matrix", "character"),
                      defaultvalues = list(NULL, "DyadSet"),
                      required = c(TRUE, FALSE))

  impl <- match.arg(a$implementation, c("Network","DyadSet"))
  
  list(name=paste0("_discord_net_",impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..intersect.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "assume_all_toggles_in_list", "implementation"),
                      vartypes = c("network,matrix", "logical", "character"),
                      defaultvalues = list(NULL, FALSE, "DyadSet"),
                      required = c(TRUE, FALSE, FALSE))
  
  impl <- match.arg(a$implementation, c("Network","DyadSet"))

  list(name=paste0(if(a$assume_all_toggles_in_list) "_intersect_net_toggles_in_list_" else "_intersect_net_", impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..union.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "implementation"),
                      vartypes = c("network,matrix", "character"),
                      defaultvalues = list(NULL, "DyadSet"),
                      required = c(TRUE, FALSE))

  impl <- match.arg(a$implementation, c("Network","DyadSet"))
  
  list(name=paste0("_union_net_",impl),
       coef.names=c(),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       dependence=FALSE)
}
