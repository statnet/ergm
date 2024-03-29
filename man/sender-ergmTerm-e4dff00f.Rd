% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.R, R/InitWtErgmTerm.R
\name{sender-ergmTerm}
\alias{sender-ergmTerm}
\alias{InitErgmTerm.sender}
\alias{InitWtErgmTerm.sender}
\title{Sender effect}
\usage{
# binary: sender(base=1, nodes=-1)

# valued: sender(base=1, nodes=-1, form="sum")
}
\arguments{
\item{base}{deprecated}

\item{nodes}{specify which nodes' statistics should be included or excluded (see Specifying Vertex attributes and Levels (\code{?nodal_attributes}) for details)}

\item{form}{character how to aggregate tie values in a valued ERGM}
}
\description{
This term adds one network statistic for each node equal to the number of
out-ties for that node. This measures the activity of the node. The term for
the first node is omitted by default because of linear dependence that
arises if this term is used together with \code{edges} , but its coefficient
can be computed as the negative of the sum of the coefficients of all the
other actors. That is, the average coefficient is zero, following the
Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
1981).

For undirected networks, see \code{sociality} .
}
\note{
The argument \code{base} is retained for backwards compatibility and may be
removed in a future version. When both \code{base} and \code{nodes} are passed,
\code{nodes} overrides \code{base}.

This term can only be used with directed networks.
}
\seealso{
\code{\link{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "sender", "subsection")}
}
\concept{directed}
\concept{dyad-independent}
