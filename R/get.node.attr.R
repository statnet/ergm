#  File ergm/R/get.node.attr.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
#  get.node.attr takes the network.and a vertex attribute name, as
#  well as the name of the function calling get.node.attr, as arguments.
#  If there is a vertex attribute by that name, it returns the vector
#  of covariates; otherwise, it aborts and prints an error message.

# This is a kludge, which has been patched to bring it in line with the
# corrected class definitions.  -CTB

get.node.attr <- function(nw, attrname, functionname=NULL, numeric=FALSE) {  
  if (is.null(functionname)) {
    # Assume it's being called from InitErgm.* or InitErgmTerm.*
    # Otherwise,  get.InitErgm.fname() will return NULL
    functionname <- get.InitErgm.fname()
    functionname <- ifelse (is.null(functionname), "unknown function",
                        sub('.*[.]', '', functionname)) # truncate up to last '.'
  }
  if (!is.character(attrname) || length(attrname)>1)
    stop(paste("The argument", attrname, "passed to", functionname,
               "must be a single character string naming a nodal attribute."),
         call.=FALSE)
  #We'll assume that every vertex must have a value, so checking the 
  #first is reasonable.
  if (!any(attrname==names(nw$val[[1]])))
    stop(paste("Attribute", attrname, "named in", functionname,
               "model term is not contained in vertex attribute list."),
         call.=FALSE)
  #"[["(nw$val,attrname)
  out <- unlist(get.vertex.attribute(nw,attrname))
  if(numeric && !is.numeric(out)) {
    stop("The ", attrname, " attribute for the ", functionname, 
         " term is not numeric as required.", call.=FALSE)
  }
  out
}

