\name{vector.namesmatch}
\alias{vector.namesmatch}
\title{
reorder vector v into order determined by matching the names of its elements to a vector of names
}
\description{
A helper function to reorder vector \code{v} (if named) into order specified by matching its names to the argument \code{names}
}
\usage{
vector.namesmatch(v, names, errname = NULL)
}
\arguments{
  \item{v}{
a vector (or list) with named elements, to be reorderd
} 
  \item{names}{
a character vector of element names, corresponding to names of \code{v}, specificying desired orering of \code{v}
}
  \item{errname}{
   optional, name to be reported in any error messages. default to \code{deparse(substitute(v))}
}
}
\details{
does some checking of appropriateness of arguments, and reorders v by matching its names to character vector \code{names}
}
\value{
returns \code{v}, with elements reordered
}

\note{
earlier versions of this function did not order as advertiased
}


\examples{
test<-list(c=1,b=2,a=3)
vector.namesmatch(test,names=c('a','c','b'))
}
