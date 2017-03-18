%  File man/get.node.attr.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{get.node.attr}
\alias{get.node.attr}
\title{
Retrieve and check assumptions about vertex attributes (nodal covariates) in a network
}
\description{
The \code{get.node.attr} function returns the vector of nodal covariates for the given
network and specified attribute if the attribute exists - execution will
halt if the attribute is not correctly given as a single string or is not 
found in the vertex attribute list; optionally \code{get.node.attr} will also 
check that return vector is numeric, halting execution if not. The purpose is to validate assumptions before passing attribute data into an ergm term. 
}
\usage{
get.node.attr(nw, attrname, functionname = NULL, numeric = FALSE)
}
\arguments{
  \item{nw}{
a \code{\link{network}} object
}
  \item{attrname}{
the name of a nodal attribute, as a character string
}
  \item{functionname}{
the name of the calling function a character string; this is only used for
                the warning messages that accompany a halt
}
  \item{numeric}{
logical, whether to halt execution if the return vector is not
           numeric; default=FALSE
}
}

\value{
returns the vector of 'attrname' covariates for the vertices in the network
}

\seealso{
\code{\link[network]{get.vertex.attribute}} for a version without the checking functionality
}
\examples{
data(faux.mesa.high)
get.node.attr(faux.mesa.high,'Grade')
}
