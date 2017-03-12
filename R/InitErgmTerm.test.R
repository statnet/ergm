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
       auxiliaries = ~.sociomatrix(mode))
}

InitErgmTerm..discord.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x"),
                      vartypes = c("network,matrix"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  list(name="_discord_net",
       coef.names=c(),
       inputs=ergm.Cprepare.el(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..intersect.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x"),
                      vartypes = c("network,matrix"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  list(name="_intersect_net",
       coef.names=c(),
       inputs=ergm.Cprepare.el(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm..union.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x"),
                      vartypes = c("network,matrix"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  list(name="_union_net",
       coef.names=c(),
       inputs=ergm.Cprepare.el(a$x, prototype=nw),
       dependence=FALSE)
}

InitErgmTerm.discord.inter.union.net <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("x"),
                      vartypes = c("network"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  nedges <- network.edgecount(a$x)
  
  
  list(name="disc_inter_union_net",
       coef.names=c("Diun","dIun","diUn","Diun2","dIun2","diUn2"),
       auxiliaries = ~ .discord.net(a$x) + .intersect.net(a$x) + .union.net(a$x),
       inputs=ergm.Cprepare.el(a$x, prototype=nw),
       emptynwstats=c(nedges, 0, nedges, nedges^2, 0, nedges^2),
       dependence=FALSE)
}
