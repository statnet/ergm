\name{search.ergmTerms}
\alias{search.ergmTerms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Search the ergm-terms documentation for appropriate terms
}
\description{
Searches through the \code{\link{ergm.terms}} help page and prints out a list of terms appropriate for the specified network's structural constraints, optionally restricting by additional categories and keyword matches. 
}
\usage{
search.ergmTerms(keyword, net,  categories, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{keyword}{
optional character keyword to search for in the text of the term descriptions. Only matching terms will be returned. Matching is case insensitive.
}
  \item{net}{
a network object that the term would be applied to, used as template to determine directedness, bipartite, etc
}
  \item{categories}{
optional character vector of category tags to use to restrict the results (i.e. 'curved', 'triad-related')
}
\item{name}{
optional character name of a specific term to return
}
}
\details{
Uses \code{\link{grep}} internally to match keywords against the term description, so \code{keywords} is currently matched as a single phrase.  Category tags will only return a match if all of the specified tags are included in the term. 
}
\value{
prints out the name and short description of matching terms, and invisibly returns them as a list.  If \code{name} is specified, prints out the full definition for the named term. 
}

\author{
skyebend@uw.edu
}

\seealso{
See also \code{\link{ergm.terms}} for the complete documentation
}
\examples{
# find all of the terms that mention triangles
search.ergmTerms('triangle')

# two ways to search for bipartite terms:

# search using a bipartite net as a template
myNet<-network.initialize(5,bipartite=3)
search.ergmTerms(net=myNet)

# or request the bipartite category
search.ergmTerms(categories='bipartite')

# search on multiple categories
search.ergmTerms(categories=c('bipartite','dyad-independent'))

# print out the content for a specific term
search.ergmTerms(name='b2factor')

}



