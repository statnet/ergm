#  File ergm/R/InitErgmTerm.duration.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
InitErgmTerm.edges.ageinterval<-function(nw, arglist, role, ...) {
  if(!any(role %in% c("dissolution","target"))) stop("Term edges.ageinterval can only be used in a dissolution model or as a target statistic.")
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("from","to"),
                      vartypes = c("numeric","numeric"),
                      defaultvalues = list(NULL, Inf),
                      required = c(TRUE, FALSE))
     
  from<-a$from
  to<-a$to

  if(from<1) stop("An extant edge cannot have an \"age\" of less than 1.")
  list(name=if(role=="target") "edges_ageinterval_mon" else "edges_ageinterval",
       coef.names = paste("edges","age",from,"to",to,sep="."),
       inputs=c(from, if(to==Inf) 0 else to),
       dependence=FALSE)
}

InitErgmTerm.edge.ages<-function(nw, arglist, role, ...) {
  if(role!="target") stop("Term edge.ages can only be used as a target statistic.")
  a <- check.ErgmTerm(nw, arglist)
  
  list(name="edge_ages_mon",
       coef.names = "edge.ages",
       dependence=FALSE)
}

InitErgmTerm.edgecov.ages<-function(nw, arglist, role, ...) {
  if(role!="target") stop("Term edgecov.ages can only be used as a target statistic.")

  a <- check.ErgmTerm(nw, arglist, 
                      varnames = c("x", "attrname"),
                      vartypes = c("matrix,network", "character"),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, FALSE))
  
  ### Check the network and arguments to make sure they are appropriate.
  ### Process the arguments
  if(is.network(a$x))
    xm<-as.matrix(a$x,matrix.type="adjacency",a$attrname)
  else if(is.character(a$x))
    xm<-get.network.attribute(nw,a$x)
  else
    xm<-as.matrix(a$x)
  
  ### Construct the list to return
  if(!is.null(a$attrname)) {
    # Note: the sys.call business grabs the name of the x object from the 
    # user's call.  Not elegant, but it works as long as the user doesn't
    # pass anything complicated.
    cn<-paste("edgecov.ages", as.character(a$attrname), sep = ".")
  } else {
    cn<-paste("edgecov.ages", as.character(sys.call(0)[[3]][2]), sep = ".")
  }
  inputs <- c(NCOL(xm), as.double(xm))
  attr(inputs, "ParamsBeforeCov") <- 1
  list(name="edgecov_ages_mon", coef.names = cn, inputs = inputs, dependence=FALSE
       )
}
