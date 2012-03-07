%  File ergm/man/edgelist.ergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{edgelist.ergm}
\alias{edgelist.ergm}
\alias{edgelist.ergm.default}
\alias{edgelist.ergm.network}
\alias{edgelist.ergm.matrix}
\title{Return edgelist in standard ergm format}
\description{
  \code{\link{edgelist.ergm}} returns an edgelist for a network in a format
  that is expected by many of the routines of the ergm package.
}
\usage{
\method{edgelist.ergm}{default}(x,\dots)
\method{edgelist.ergm}{network}(x, \dots)
\method{edgelist.ergm}{matrix}(x, directed=TRUE, check.uniqueness=TRUE,
                               check.sorted=TRUE, \dots)
}
\arguments{
 \item{x}{an \R object. Either a network or a matrix of some sort.  }
 \item{directed}{logical:  Are the edges directed?  If x is a square 
 matrix with NROWS(x)!=2, equal to !all(x==t(x))}
 \item{check.uniqueness}{An edgelist should have unique rows.  By default,
 therefore, uniqueness is checked.  However, if the input is known to 
 have unique rows already, then setting this option to FALSE will save
 computing time}
 \item{check.sorted}{An edgelist should have its rows sorted in dictionary
 order.  By default, therefore, the matrix is checked to see whether it is
 sorted.  However, if the input is known to be sorted already, then
 setting this option to FALSE will save computing time.}
 \item{\dots}{Additional
 arguments, to be passed to lower-level functions in the future.}
}

\details{
  The function takes a network or an (adjacency or edgelist) matrix, then
  returns an edgelist (of class "matrix") 
  in standard format.  If x is a network, the value returned
  equals \code{edgelist.ergm(as.matrix(x, "edgelist"), directed=is.directed(x))}.
  
  The standard format is as follows:  

1. The matrix has two columns

2. No row has two identical entries

3. Each row is unique

4. The rows are in dictionary order:  They are sorted by the first column,
  then by the second in case of ties
  
5. If \code{directed=TRUE}, the element in the first column is always
  smaller than the element in the second column (otherwise, the entries in that
  row are switched before sorting).
}
\value{
A matrix with two rows, in the format described under "Details"
      }

\seealso{\code{\link[network]{as.matrix.network}}}
\keyword{models}
