#  File ergm/R/get.node.attr.R
#  Part of the statnet package, http://statnet.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnet.org/attribution
#
#  Copyright 2012 the statnet development team
######################################################################
###############################################################################
# The <get.node.attr> function returns the vector of covariates for the given
# network and specified attribute if the attribute exists - execution will
# halt if the attribute is not correctly given as a single string or is not 
# found in the vertex attribute list; optionally <get.node.attr> will also 
# check that return vector is numeric, halting execution if not
###############################################################################

get.node.attr <- function(nw, attrname, functionname=NULL, numeric=FALSE) {  

# This is a kludge, which has been patched to bring it in line with the
# corrected class definitions.  -CTB

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

