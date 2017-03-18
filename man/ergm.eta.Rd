%  File man/ergm.eta.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.eta}
\alias{ergm.eta}
\alias{ergm.etagrad}
\alias{ergm.etagradmult}
\alias{ergm.etamap}

\title{
Operations with 'eta' vector of canonical parameter values from ergm model 
}
\description{
The \code{ergm.eta} function calculates and returns eta, mapped from
theta using the etamap object created by \code{ergm.etamap}.

The \code{ergm.etagrad} function caculates and returns the gradient of eta
 mapped from theta using the etamap object created by \code{ergm.etamap}. If the
 gradient is only intended to be a multiplier for some vector, the more
 efficient \code{ergm.etagradmult} is recommended.
 
The \code{ergm.etagradmult} function calculates and returns the product of the
 gradient of eta with a vector v
 
 The \code{ergm.etamap} function takes a model object and creates a mapping
 from the model parameters, theta, to the canonical (linear) eta parameters;
 the mapping is carried out by \code{ergm.eta}

}
\usage{
ergm.eta(theta, etamap)

ergm.etagrad(theta, etamap)

ergm.etagradmult(theta, v, etamap)

ergm.etamap(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{the curved model parameters }
  \item{etamap}{the list of values that constitutes the theta-> eta mapping and is returned by \code{ergm.etamap}}
  \item{v}{a vector of the same length as the vector of mapped eta parameters}
  \item{model}{ model object, as returned by \code{\link{ergm.getmodel}}}
}
\details{
  This function is only important in the case of curved exponential family models, i.e., those in which the parameter of interest (theta) is not a linear function of the sufficient statistics (eta) in the exponential-family model. In non-curved models, we may assume without loss of generality that eta(theta)=theta.  
  
  A succinct description of how eta(theta) is incorporated into an ERGM is given by equation (5) of Hunter (2007).  See Hunter and Handcock (2006) and Hunter (2007) for further details about how eta and its derivatives are used in the estimation process.
}
\value{
\itemize{
\item for \code{ergm.eta}: \code{eta}  the canonical eta parameters as mapped from theta
\item for \code{ergm.etagrad}: \code{etagrad} a matrix of the gradient of eta 
\item{ for \code{ergm.etagradmult}: \code{ans} the vector that is the product of the gradient of eta and \code{v}; infinite values are replaced by (+-)10000}
\item{ for \code{ergm.etamap} the theta -> eta mapping given by a list of the following:
\itemize{
     \item{canonical  : a numeric vector whose ith entry specifies whether
                 the ith component of theta is canonical (via non-
                  negative integers) or curved (via zeroes)}
     \item{offsetmap  : a logical vector whose ith entry tells whether the
                  ith coefficient of the canonical parameterization
                  was "offset", i.e fixed }
     \item{offset     : a logical vector whose ith entry tells whether the
                  ith model term was offset/fixed}
     \item{offsettheta: a logical vector whose ith entry tells whether the
                  ith curved theta coeffient was offset/fixed;}
     \item{curved     : a list with one component per curved EF term in the
                  model containing
        \itemize{        
         \item{from    : the indices of the curved theta parameter that are
                   to be mapped from}
         \item{to      : the indices of the canonical eta parameters to be
                   mapped to}
         \item{map     : the map provided by <InitErgmTerm>}
         \item{gradient: the gradient function provided by \link{InitErgmTerm} }
         \item{cov     : the eta covariance ??, possibly always NULL (no
                   <Init> function creates such an item)}
          }
                   
      }             
     \item{etalength  : the length of the eta vector}

}
}
}
}
\references{
  \itemize{
    \item Hunter, D. R. and M. S. Handcock (2006). Inference in curved
  exponential family models for networks. \emph{Journal of Computational
  and Graphical Statistics}, 15: 565--583.

    \item Hunter, D. R. (2007). Curved exponential family models for social
  networks. \emph{Social Networks}, 29: 216--230.
  }
}


\seealso{
\code{\link{ergm-terms}}
}
