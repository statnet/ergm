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
                      varnames = c("x", "assume_all_toggles_in_list"),
                      vartypes = c("network,matrix", "logical"),
                      defaultvalues = list(NULL, FALSE),
                      required = c(TRUE, FALSE))
  
  list(name=if(a$assume_all_toggles_in_list) "_intersect_net_toggles_in_list" else "_intersect_net",
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

## Exports a network that tracks the current network but excludes ties outside of the blocks specified by attrname.

InitErgmTerm..blockdiag.net <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attrname"),
                      vartypes = c("character"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  ### Process the arguments
    ### Process the arguments
  nodecov <-
    if(length(a$attrname)==1)
      get.node.attr(nw, a$attrname)
    else{
      do.call(paste,c(sapply(a$attrname,function(oneattr) get.node.attr(nw,oneattr),simplify=FALSE),sep="."))
    }
  u <- sort(unique(nodecov))
  nodecov <- match(nodecov,u)

  list(name="_blockdiag_net", coef.names = c(), inputs=nodecov, dependence=FALSE)
}
