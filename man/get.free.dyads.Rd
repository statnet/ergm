\name{get.free.dyads}
\alias{get.free.dyads}
\alias{get.miss.dyads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a network containing only edges meeting a specific criteria
}
\description{
\code{get.free.dyads} will create a network object with only edges that are not targeted by \link{ergm-constraints}. \code{get.miss.dyads} creates a network object with only edges that are missing/not present (not missing in the NA sense) in a network subject to constraints
}
\usage{
get.free.dyads(constraints)

get.miss.dyads(constraints, constraints.obs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{constraints}{
%%     ~~Describe \code{constraints} here~~
}
  \item{constraints.obs}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A network object containing the specified set of edges
}

\seealso{
\link{ergm-constraints}
}

\examples{

}
