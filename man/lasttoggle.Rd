%  File man/lasttoggle.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{lasttoggle}
\alias{lasttoggle}
\alias{last.toggle}
\alias{last-toggle}
\alias{to.matrix.lasttoggle}
\alias{to.lasttoggle.matrix}
\alias{ergm.el.lasttoggle}
\title{Storing last toggle information in a network}
\description{
  An informal extension to \code{\link{network}} objects allowing some
  limited temporal information to be stored.
}
\details{

  WARNING: THIS DOCUMENTATION IS PROVIDED AS A COURTESY, AND THE API
  DESCRIBED IS SUBJECT TO CHANGE WITHOUT NOTICE, DOWN TO COMPLETE
  REMOVAL. NOT ALL FUNCTIONS THAT COULD SUPPORT IT DO. USE AT YOUR OWN
  RISK.
  
  While \code{\link[networkDynamic]{networkDynamic}} provides a flexible, consistent
  method for storing dynamic networks, the \code{C} routines of
  \code{\link[=ergm-package]{ergm}} and
  \code{\link[tergm:tergm-package]{tergm}} required a simpler and more
  lightweight representation.

  This representation consisted of a single integer representing the
  time stamp and an integer vector of length to
  \code{\link{network.dyadcount}(nw)} --- the number of potential ties
  in the network, giving the last time point during which each of the
  dyads in the network had changed.

  Though this is an API intended for internal use, some functions,
  like \code{\link[tergm]{stergm}} (for EGMME),
  \code{\link[tergm:simulate.stergm]{simulate}}, and
  \code{\link[=summary.formula]{summary}} can be passed networks with
  this information using the following \code{\link{network}} (i.e.,
  \code{\link{\%n\%}}) attributes:
  \describe{
    \item{\code{time}}{the time stamp associated with the network}
    \item{\code{lasttoggle}}{a vector of length
      \code{\link{network.dyadcount}(nw)}, giving the last change time
      associated with each dyad. See the source code of
      \code{\link[=ergm-package]{ergm}} internal functions
      \code{to.matrix.lasttoggle}, \code{ergm.el.lasttoggle}, and
      \code{to.lasttoggle.matrix} for how they are serialized.}
  }

  For technical reasons, the \code{\link[tergm:tergm-package]{tergm}} routines
  treat the \code{lasttoggle} time points as shifted by \eqn{-1}.

  Again, this API is subject to change without notice.
}