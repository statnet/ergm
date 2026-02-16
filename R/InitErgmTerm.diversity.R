#  File R/InitErgmTerm.diversity.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' @templateVar name nodecovrange
#' @title Range of covariate values for neighbors of a node
#' @description This term adds a single network statistic equalling
#'   the sum over the nodes of the range over of its neighbors'
#'   values.
#'
#' @details This is a network analogue of the statistic introduced by
#'   \insertCite{HoBl23m;textual}{ergm}.
#' 
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept undirected
#' @concept quantitative nodal attribute
InitErgmTerm.nodecovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodecovrange")
  list(name="nodecovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name nodeocovrange
#' @title Range of covariate values for out-neighbors of a node
#'
#' @usage
#' # binary: nodeocovrange(attr)
#'
#' @inherit nodecovrange-ergmTerm
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept quantitative nodal attribute
InitErgmTerm.nodeocovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodeocovrange")
  list(name="nodeocovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name nodeicovrange
#' @title Range of covariate values for in-neighbors of a node
#'
#' @usage
#' # binary: nodeicovrange(attr)
#'
#' @inherit nodecovrange-ergmTerm
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept quantitative nodal attribute
InitErgmTerm.nodeicovrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=TRUE, bipartite=NULL,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric")
  coef.names <- nodecov_names(nodecov, "nodeicovrange")
  list(name="nodeicovrange", coef.names=coef.names, inputs=c(nodecov))
}


#' @templateVar name b1covrange
#' @title Range of covariate values for neighbors of a mode-1 node
#'
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @inherit nodecovrange-ergmTerm
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept bipartite
#' @concept quantitative nodal attribute
InitErgmTerm.b1covrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip="b2")
  coef.names <- nodecov_names(nodecov, "b1covrange")
  list(name="b1covrange", coef.names=coef.names, inputs=c(nodecov))
}



#' @templateVar name b2covrange
#' @title Range of covariate values for neighbors of a mode-2 node
#'
#' @usage
#' # binary: nodecovrange(attr)
#'
#' @inherit nodecovrange-ergmTerm
#' @template ergmTerm-attr
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept bipartite
#' @concept quantitative nodal attribute
InitErgmTerm.b2covrange<-function (nw, arglist, ...) {
### Check the network and arguments to make sure they are appropriate.
  a <- check.ErgmTerm(nw, arglist, directed=NULL, bipartite=TRUE,
                      varnames = c("attr"),
                      vartypes = c(ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
### Process the arguments
  nodecov <- ergm_get_vattr(a$attr, nw, accept="numeric", bip="b1")
  coef.names <- nodecov_names(nodecov, "b2covrange")
  list(name="nodeicovrange", coef.names=coef.names, inputs=c(nodecov))
}



.nodefactordistinct_impl <- function(deg, dir, bip, nw, arglist, ..., degname=deg){
  a <- check.ErgmTerm(nw, arglist, directed=dir, bipartite=bip,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))

  attr <- a$attr
  levels <- a$levels

  nodecov <- if(NVL(bip, FALSE)) ergm_get_vattr(attr, nw, bip = c(b1="b2",b2="b1")[deg]) else ergm_get_vattr(attr, nw)
  attrname <- attr(nodecov, "name")
  u <- ergm_attr_levels(levels, nodecov, nw, levels = sort(unique(nodecov)))

  if (length(u)==0) { # Get outta here!  (can happen if user passes attribute with one value)
    return()
  }
  #   Recode to numeric
  nodepos <- match(nodecov,u,nomatch=0)
  ### Construct the list to return
  inputs <- c(max(nodepos), nodepos)
  list(name=paste0(degname, "factordistinct"),                                        #required
       coef.names = paste(paste0(deg, "factordistinct"), paste(attrname,collapse="."), sep="."), #required
       iinputs = inputs,
       minval = 0
       )
}

#' @templateVar name nodefactordistinct
#' @title Number of distinct neighbor types
#' @description This term adds a single network statistic to the
#'   model, counting, for each node, the number of distinct values of
#'   the attribute found among its neighbors.
#'
#' @details This is a network analogue of the statistic introduced by
#'   \insertCite{HoBl23m;textual}{ergm}.
#'
#' @usage
#' # binary: nodefactordistinct(attr, levels=TRUE)
#'
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept undirected
#' @concept categorical nodal attribute
InitErgmTerm.nodefactordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                      dep.inform = list(FALSE, FALSE))
  .nodefactordistinct_impl("node", NULL, NULL, nw, arglist)
}


#' @templateVar name nodeofactordistinct
#' @title Number of distinct out-neighbor types
#'
#' @usage
#' # binary: nodeofactordistinct(attr, levels=TRUE)
#'
#' @inherit nodefactordistinct-ergmTerm
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.nodeofactordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                        dep.inform = list(FALSE, FALSE))
  .nodefactordistinct_impl("nodeo", TRUE, NULL, nw, arglist)
}


#' @templateVar name nodeifactordistinct
#' @title Number of distinct in-neighbor types
#'
#' @usage
#' # binary: nodeifactordistinct(attr, levels=TRUE)
#'
#' @inherit nodefactordistinct-ergmTerm
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept directed
#' @concept categorical nodal attribute
InitErgmTerm.nodeifactordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                        dep.inform = list(FALSE, FALSE))
  .nodefactordistinct_impl("nodei", TRUE, NULL, nw, arglist)
}


#' @templateVar name b1factordistinct
#' @title Number of distinct neighbor types for the first node
#'
#' @usage
#' # binary: b1factordistinct(attr, levels=TRUE)
#'
#' @inherit nodefactordistinct-ergmTerm
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept bipartite
#' @concept categorical nodal attribute
InitErgmTerm.b1factordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                        dep.inform = list(FALSE, FALSE))
  .nodefactordistinct_impl("b1", NULL, TRUE, nw, arglist)
}


#' @templateVar name b2factordistinct
#' @title Number of distinct neighbor types for the second mode
#'
#' @usage
#' # binary: b2factordistinct(attr, levels=TRUE)
#'
#' @inherit nodefactordistinct-ergmTerm
#' @template ergmTerm-attr
#' @templateVar explain this optional argument controls which levels of the attribute
#'   should be included and which should be excluded.
#' @template ergmTerm-levels-doco
#'
#' @template ergmTerm-general
#'
#' @references \insertAllCited{}
#' @concept bipartite
#' @concept categorical nodal attribute
InitErgmTerm.b2factordistinct<-function (nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "levels"),
                        vartypes = c(ERGM_VATTR_SPEC, ERGM_LEVELS_SPEC),
                        defaultvalues = list(NULL, NULL),
                        required = c(TRUE, FALSE),
                        dep.inform = list(FALSE, FALSE))
  .nodefactordistinct_impl("b2", NULL, TRUE, nw, arglist)
}
