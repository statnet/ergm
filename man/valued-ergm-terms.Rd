\name{valued-ergm-terms}
\alias{valued-ergm-terms}
\alias{valued.ergm.terms}
\alias{valued-terms-ergm}
\alias{valued.terms.ergm}
\alias{valued-terms}
\alias{valued.terms}
\alias{InitWtErgmTerm.absdiff}
\alias{InitWtErgmTerm.absdiffcat}
\alias{InitWtErgmTerm.atleast}
\alias{InitWtErgmTerm.cyclicalties}
\alias{InitWtErgmTerm.cyclicalweights}
\alias{InitWtErgmTerm.deference}
\alias{InitWtErgmTerm.edgecov}
\alias{InitWtErgmTerm.greaterthan}
\alias{InitWtErgmTerm.inconsistency}
\alias{InitWtErgmTerm.ininterval}
\alias{InitWtErgmTerm.mutual}
\alias{InitWtErgmTerm.nodecov}
\alias{InitWtErgmTerm.nodecovar}
\alias{InitWtErgmTerm.nodefactor}
\alias{InitWtErgmTerm.nodeicov}
\alias{InitWtErgmTerm.nodeicovar}
\alias{InitWtErgmTerm.nodeifactor}
\alias{InitWtErgmTerm.nodeisqrtcovar}
\alias{InitWtErgmTerm.nodemain}
\alias{InitWtErgmTerm.nodeocov}
\alias{InitWtErgmTerm.nodeocovar}
\alias{InitWtErgmTerm.nodeofactor}
\alias{InitWtErgmTerm.nodeosqrtcovar}
\alias{InitWtErgmTerm.nodesqrtcovar}
\alias{InitWtErgmTerm.nonconformity}
\alias{InitWtErgmTerm.nonzero}
\alias{InitWtErgmTerm.sum}
\alias{InitWtErgmTerm.transitiveties}
\alias{InitWtErgmTerm.transitiveweights}
\alias{absdiff}
\alias{absdiffcat}
\alias{atleast}
\alias{cyclicalties}
\alias{cyclicalweights}
\alias{deference}
\alias{edgecov}
\alias{greaterthan}
\alias{inconsistency}
\alias{ininterval}
\alias{mutual}
\alias{nodecov}
\alias{nodecovar}
\alias{nodefactor}
\alias{nodeicov}
\alias{nodeicovar}
\alias{nodeifactor}
\alias{nodeisqrtcovar}
\alias{nodemain}
\alias{nodeocov}
\alias{nodeocovar}
\alias{nodeofactor}
\alias{nodeosqrtcovar}
\alias{nodesqrtcovar}
\alias{nonconformity}
\alias{nonzero}
\alias{sum}
\alias{transitiveties}
\alias{transitiveweights}
\docType{methods}
\title{Terms used in Valued Exponential Family Random Graph Models}
\description{
  The function \code{\link{ergm}} is used to fit exponential 
  random graph models, in which
  the probability of a given valued network, \eqn{y}, on a set of nodes is 
  \eqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}, where
  \eqn{h(y)} is the reference measure,
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same 
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  This page describes the terms (elements of \eqn{g(y)}) for networks
  with weighted edges that are included with the
  \code{\link[=ergm-package]{ergm}} package. Other packages may add new
  terms. For binary ERGM terms, and for general information about
  specifying ERGMs, go \link[=ergm-terms]{here} and \code{\link{ergm}}.
  }
%\usage{
%ergm-terms()
%}
%\arguments{
%}
\section{Possible terms to represent network statistics}{

  A substantial number of binary ERGM statistics --- particularly the
  ones with dyadic indepenence --- have simple generalizations to valued
  ERGMs, and have been adapted in
  \code{\link[=ergm-package]{ergm}}. They have the same form as their
  binary ERGM counterparts listed in the \link[=ergm-terms]{binary ERGM
  terms document}, with an additional argument: \code{form}, which, at
  this time, has two possible values: \code{"sum"} (the default) and
  \code{"nonzero"}. The former creates a statistic of the form
  \eqn{\sum_{i,j} x_{i,j} y_{i,j}}, where \eqn{y_{i,j}} is the value of
  dyad \eqn{(i,j)} and \eqn{x_{i,j}} is the term's covariate associated
  with it. The latter computes the binary version, with the edge
  considered to be present if its value is not 0. The terms generalized
  in this way include the following: \code{\link{absdiff}},
  \code{\link{absdiffcat}}, \code{\link{edgecov}},
  \code{\link{nodefactor}}, \code{\link{nodeofactor}},
  \code{\link{nodeifactor}},
  \code{\link{nodecov}}/\code{\link{nodemain}}, \code{\link{nodeicov}},
  and \code{\link{nodeocov}}.

  In addition, the following terms act like their binary ERGM
  counterparts, considering an edge to be present if its value is not 0:
  \code{\link{edges}}, \code{\link{transitiveties}},
  \code{\link{cyclicalties}}.
  
  The following terms are specific to valued ERGMs. For further
  discussion, see Krivitsky (2012):
  \describe{
    
    \item{\code{atleast(threshold=0)}}{\emph{Number of ties with values
	greater than or equal to a threshold}
      Adds one statistic equaling to the number of ties whose values
      equal or exceed \code{threshold}.
    }
    
    \item{\code{ininterval(lower=-Inf, upper=+Inf, open=c(TRUE,TRUE))}}{\emph{Number of ties
	whose values are in an interval}
      Adds one statistic equaling to the number of ties whose values
      are between \code{lower} and
      \code{upper}. Argument \code{open}
      is a \code{logical} vector of length 2 that controls whether
      the interval is open (exclusive) on the lower and on the upper
      end, respectively.
    }
    
    \item{\code{greaterthan(threshold=0)}}{\emph{Number of dyads with
	values strictly greater than a threshold:}
      Adds one statistic equaling to the number of ties whose values
      exceed \code{threshold}.
    }
    
    \item{\code{sum(pow=1)}}{\emph{Sum of dyad values (optionally taken
	to a power):} This term adds one statistic equal to the sum of
      dyad values taken to the power \code{pow}, which defaults to 1.
    }
    
    \item{\code{nodecovar}}{\emph{Uncentered covariance of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} (y_{i,j}y_{i,k}+y_{k,j}y_{k,j})}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic.
    }

    \item{\code{nodeicovar}}{\emph{Uncentered covariance of in-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{k,j}y_{k,j}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.
    }

    \item{\code{nodeocovar}}{\emph{Uncentered covariance of out-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{i,j}y_{i,k}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.
    }
    
    \item{\code{nodesqrtcovar(center=TRUE)}}{\emph{Covariance of square roots of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k}
	(\sqrt{y_{i,j}}\sqrt{y_{i,k}}+\sqrt{y_{k,j}}\sqrt{y_{k,j}})} if
      \code{center=FALSE}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic. If \code{center=FALSE} (the default), the statistic is
	instead \eqn{\sum_{i,j,k}
	((\sqrt{y_{i,j}}-\bar{\sqrt{y}})(\sqrt{y_{i,k}}-\bar{\sqrt{y}})+(\sqrt{y_{k,j}}-\bar{\sqrt{y}})(\sqrt{y_{k,j}}-\bar{\sqrt{y}}))},
	where \eqn{\bar{\sqrt{y}}} is the mean of the square root of
	dyad values.
    }

    \item{\code{nodeisqrtcovar}}{\emph{Uncentered covariance of square roots of in-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{k,j}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.
    }

    \item{\code{nodeosqrtcovar}}{\emph{Uncentered covariance of square roots of out-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{i,k}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.
    }

    \item{\code{mutual(form="min",threshold=0)}}{\emph{Mutuality:} This term
      adds one statistic equal to \eqn{\sum_{i<j} m(y_{i,j},y_{j,i})},
      where \eqn{m} is determined by \code{form} argument: \code{"min"}
      for \eqn{\min(y_{i,j},y_{j,i})}, \code{"nabsdiff"} for
      \eqn{-|y_{i,j},y_{j,i}|}, \code{"product"} for
      \eqn{y_{i,j}y_{j,i}}, and \code{"geometric"} for
      \eqn{\sqrt{y_{i,j}}\sqrt{y_{j,i}}}. See Krivitsky (2012) for a
      discussion of these statistics. \code{form="threshold"} simply
      computes the binary \code{\link[=ergm-terms]{mutuality}} after
      thresholding at \code{threshold}.
    }

    \item{\code{transitiveweights(twopath="min",combine="max",affect="min")}}{\emph{Transitive
	weights:} This statistic implements the transitive weights
	statistic defined by Krivitsky (2012), Equation 13. The
	currently implemented options for \code{twopath} is the minimum
	of the constituent dyads (\code{"min"}) or their geometric mean
	(\code{"geomean"}); for \code{combine}, the maximum of the
	2-path strengths (\code{"max"}) or their sum (\code{"sum"}); and
	for \code{affect}, the minimum of the focus dyad and the
	combined strength of the two paths (\code{"min"}) or their
	geometric mean (\code{"geomean"}). For each of these options,
	the first (and the default) is more stable but also more
	conservative, while the second is more sensitive but more likely
	to induce a multimodal distribution of networks.
      }

      \item{\code{cyclicalweights(twopath="min",combine="max",affect="min")}}{\emph{Cyclical
	weights:} This statistic implements the cyclical weights
	statistic, like that defined by Krivitsky (2012), Equation 13,
	but with the focus dyad being \eqn{y_{j,i}} rather than
	\eqn{y_{i,j}}. See
	\code{\link[=valued-ergm-terms]{transitiveweights}} for the
	arguments.
      }
    }
  }

\references{
  \itemize{
    \item Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
    Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
    1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}
    
    \item Morris M, Handcock MS, and Hunter DR. Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects. \emph{Journal of Statistical Software}, 2008, 24(4), 1-24. \url{http://www.jstatsoft.org/v24/i04}

    }
    
}
\seealso{ergm, network, \%v\%, \%n\%, sna,
  summary.ergm, print.ergm}
% TODO: Write an example.
     
\keyword{models}
