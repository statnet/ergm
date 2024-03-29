% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R, R/InitWtErgmTerm.R
\name{absdiff-ergmTerm}
\alias{absdiff-ergmTerm}
\alias{InitErgmTerm.absdiff}
\alias{InitWtErgmTerm.absdiff}
\title{Absolute difference in nodal attribute}
\usage{
# binary: absdiff(attr,
#                 pow=1)

# valued: absdiff(attr,
#                 pow=1,
#                 form="sum")
}
\arguments{
\item{attr}{a vertex attribute specification (see Specifying Vertex attributes and Levels (\code{?nodal_attributes}) for details.)}

\item{pow}{power to which to take the absolute difference}

\item{form}{character how to aggregate tie values in a valued ERGM}
}
\description{
This term adds one network statistic to the model equaling
the sum of \code{abs(attr[i]-attr[j])^pow} for all edges \verb{(i,j)} in
the network.
}
\note{
\pkg{ergm} versions 3.9.4 and earlier used different arguments for this
term. See \code{\link{ergm-options}} for how to invoke the old behaviour.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "absdiff", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
\concept{quantitative nodal attribute}
\concept{undirected}
