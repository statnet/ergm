%  File man/mcmc.list_utils.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{colMeans.mcmc.list}
\alias{colMeans.mcmc.list}
\alias{sweep.mcmc.list}
\title{
utility operations for mcmc.list objects
}
\description{
Adaptations of \code{\link{colMeans}} and \code{\link{sweep}} for working with \code{\link[coda]{mcmc.list}} objects. The function mcmc.list is used to represent parallel runs of the same chain, with different starting values and random seeds. 

}
\usage{
colMeans.mcmc.list(x, ...)

sweep.mcmc.list(x, STATS, FUN = "-", check.margin = TRUE, ...)
}
\arguments{
  \item{x}{
a \code{\link{mcmc.list}} object
}
 \item{STATS}{
 the summary statistic which is to be swept out
 }
 \item{FUN}{
 the function to be used to carry out the sweep (default is '-', subtraction)
 }
 \item{check.margin}{
  logical. If TRUE (the default), warn if the length or dimensions of STATS do not match the specified dimensions of x. Set to FALSE for a small speed gain when you know that dimensions match.
 }
  \item{\dots}{
additional arguments to \code{\link{colMeans}} or \code{\link{sweep}}
}
}
\details{
\code{colMeans.mcmc.list} converts the mcmc list into a matrix and then runs \code{colMeans} on it

\code{sweep.mcmc.lists} modifies the values of the chains by computing some aggregate summary statistic function (\code{STATS}) on the entire chain and applying it to each value via \code{FUN} (default is subtraction).  See \code{\link{sweep}} for additional examples.
}
\value{
\code{colMeans.mcmc} returns a vector with length equal to the number of mcmc chains in \code{x} with the mean value for each chain.
\code{sweep.mcmc.lists} returns an appropriately modified version of \code{x}
}
