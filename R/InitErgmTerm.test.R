#  File R/InitErgmTerm.test.R in package ergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
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

InitErgmTerm.sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("mode"),
                      vartypes = c("character"),
                      defaultvalues = list("integer"),
                      required = c(FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "isociomatrix")
  n <- network.size(nw)
  tails <- rep(1:n,n)
  heads <- rep(1:n,each=n)
  list(name=name,
       coef.names=paste(tails,heads,sep="."), dependence=FALSE,
       auxiliaries = trim_env(~.sociomatrix(mode),"mode"))
}

InitErgmTerm.discord.sociomatrix<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x","mode"),
                      vartypes = c("network","character"),
                      defaultvalues = list(nw,"integer"),
                      required = c(FALSE,FALSE))

  mode <- match.arg(a$mode, c("integer"))
  name <- switch(mode,
                 integer = "discord_isociomatrix")
  n <- network.size(nw)
  tails <- rep(1:n,n)
  heads <- rep(1:n,each=n)
  x <- a$x
  list(name=name,
       coef.names=paste(tails,heads,sep="."), dependence=FALSE,
       auxiliaries = trim_env(~.discord.sociomatrix(x,mode), c("x","mode")),
       emptynwstats=as.matrix(a$x,matrix.type="adjacency"))
}

InitErgmTerm.discord.inter.union.net <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x", "implementation"),
                      vartypes = c("network,matrix", "character"),
                      defaultvalues = list(NULL, "DyadSet"),
                      required = c(TRUE, FALSE))

  impl <- match.arg(a$implementation, c("Network","DyadSet"))
  x <- a$x
  nedges <- network.edgecount(x)

  list(name=paste0("disc_inter_union_net_", impl),
       coef.names=c("Diun","dIun","diUn","Diun2","dIun2","diUn2"),
       auxiliaries = trim_env(~ .discord.net(x, implementation=impl) + .intersect.net(x, implementation=impl) + .union.net(x, implementation=impl), c("x","impl")),
       inputs=to_ergm_Cdouble(a$x, prototype=nw),
       emptynwstats=c(nedges, 0, nedges, nedges^2, 0, nedges^2),
       dependence=TRUE)
}

InitErgmTerm..edges_times<-function(nw, arglist, ..., times) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="_edges_times", coef.names=paste("edges_times",times), dependence=FALSE,
       minval = min(times*network.dyadcount(nw,FALSE),0), maxval = max(times*network.dyadcount(nw,FALSE),0), conflicts.constraints="edges", inputs=c(times))
}
