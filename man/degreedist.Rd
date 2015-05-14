\name{degreedist}
\alias{degreedist}


\title{
Computes and Returns the Degree Distribution Information for a Given Network
}
\description{
 The \code{degreedist} function computes and returns the degree distribution (number of vertices in the network with each degree value) for a given network.
 
 The \code{degreedistfactor} function returns the cross table of the degree
 distribution for a network and a given factor (vertex attribute name)
}
\usage{
degreedist(g, print = TRUE)

}
\arguments{
  \item{g}{
a \code{network} object
}
  \item{print}{
logical, whether to print the degree distribution; default=TRUE
}

}
\details{
Calculates the degrees using the appropriate ergm terms for for network type. 
}
\value{


  \itemize{
    \item{if directed}{  -- a matrix of the distributions of in and out degrees;
                         this is row bound and only contains degrees for which
                       one of the in or out distributions has a positive count}
    \item{if bipartite}{ -- a list containing the degree distributions of b1 and b2}
    \item{otherwise}    {-- a vector of the positive values in the degree
                       distribution}
   }


}



\examples{
data(faux.mesa.high)
degreedist(faux.mesa.high)
}
