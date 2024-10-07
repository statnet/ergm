%  File man/macros/ergm-spec.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2024 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\ERGMspec}{In an exponential-family random graph model (ERGM), the probability or density of a given network, \eqn{y \in Y}{y in sample space Y}, on a set of nodes is \deqn{h(y) \exp[\eta(\theta) \cdot g(y)] / \kappa(\theta),}{h(y) exp{eta(theta).g(y)} / c(theta),} where \eqn{h(y)} is the reference distribution (particularly for valued network models), \eqn{g(y)} is a vector of network statistics for \eqn{y}, \eqn{\eta(\theta)} is a natural parameter vector of the same length (with \eqn{\eta(\theta)\equiv\theta}{eta(theta)=theta} for most terms), \eqn{\cdot}{.} is the dot product, and \eqn{\kappa(\theta)}{c(theta)} is the normalizing constant for the distribution. A complete ERGM specification requires a list of network statistics \eqn{g(y)} and (if applicable) their \eqn{\eta(\theta)} mappings provided by a formula of \code{\link{ergmTerm}}s; and, optionally, sample space \eqn{\mathcal{Y}}{Y} and reference distribution \eqn{h(y)} information provided by \code{\link{ergmConstraint}}s and, for valued ERGMs, by \code{\link{ergmReference}}s.}