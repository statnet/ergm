% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ergm-package.R
\docType{package}
\name{ergm-package}
\alias{ergm-package}
\title{ergm: Fit, Simulate and Diagnose Exponential-Family Models for Networks}
\description{
An integrated set of tools to analyze and simulate networks based on exponential-family random graph models (ERGMs). 'ergm' is a part of the Statnet suite of packages for network analysis. See Hunter, Handcock, Butts, Goodreau, and Morris (2008) \doi{10.18637/jss.v024.i03} and Krivitsky, Hunter, Morris, and Klumb (2023) \doi{10.18637/jss.v105.i06}.
}
\details{
\insertNoCite{HuHa08e,KrHu23e}{ergm}

For a complete list of the functions, use \code{library(help="ergm")} or
read the rest of the manual. For a simple demonstration, use
\code{demo(packages="ergm")}.

When publishing results obtained using this package, please cite the
original authors as described in \code{citation(package="ergm")}.

All programs derived from this package must cite it. Please see the
file \code{LICENSE} and \url{https://statnet.org/attribution}.

Recent advances in the statistical modeling of random networks have had an
impact on the empirical study of social networks. Statistical exponential
family models (Strauss and Ikeda 1990) are a generalization of the Markov
random network models introduced by \insertCite{FrSt86m;textual}{ergm}, which in turn
derived from developments in spatial statistics \insertCite{Be74s}{ergm}. These models
recognize the complex dependencies within relational data structures.  To
date, the use of stochastic network models for networks has been limited by
three interrelated factors: the complexity of realistic models, the lack of
simulation tools for inference and validation, and a poor understanding of
the inferential properties of nontrivial models.

This manual introduces software tools for the representation, visualization,
and analysis of network data that address each of these previous
shortcomings.  The package relies on the \code{\link[network:network]{network}}
package which allows networks to be represented in . The
\CRANpkg{ergm} package implements maximum likelihood
estimates of ERGMs to be calculated using Markov Chain Monte Carlo (via
\code{\link[=ergm]{ergm()}}). The package also provides tools for simulating networks
(via \code{\link[=simulate.ergm]{simulate.ergm()}}) and assessing model goodness-of-fit (see
\code{\link[=mcmc.diagnostics]{mcmc.diagnostics()}} and \code{\link[=gof.ergm]{gof.ergm()}}).

A number of Statnet Project packages extend and enhance
\CRANpkg{ergm}. These include
\CRANpkg{tergm} (Temporal ERGM), which provides
extensions for modeling evolution of networks over time;
\CRANpkg{ergm.count}, which facilitates
exponential family modeling for networks whose dyadic measurements are
counts; and
\pkg{ergm.userterms}, available on GitHub at \url{https://github.com/statnet/ergm.userterms}, which
allows users to implement their own ERGM terms.

For detailed information on how to download and install the software, go to
the \CRANpkg{ergm} website: \url{https://statnet.org}. A
tutorial, support newsgroup, references and links to further resources are
provided there.
}
\references{
\insertAllCited{}

Admiraal R, Handcock MS (2007).  \CRANpkg{networksis}: Simulate
bipartite graphs with fixed marginals through sequential importance
sampling.  Statnet Project, Seattle, WA.  Version 1,
\url{https://statnet.org}.

Bender-deMoll S, Morris M, Moody J (2008).  Prototype Packages for Managing
and Animating Longitudinal Network Data: \pkg{dynamicnetwork} and
\pkg{rSoNIA}.  \emph{Journal of Statistical Software}, 24(7).
\doi{10.18637/jss.v024.i07}

Boer P, Huisman M, Snijders T, Zeggelink E (2003).  StOCNET: an open
software system for the advanced statistical analysis of social networks.
Groningen: ProGAMMA / ICS, version 1.4 edition.

Butts CT (2007).  \CRANpkg{sna}: Tools for Social Network Analysis.  R package
version 2.3-2. \url{https://cran.r-project.org/package=sna}

Butts CT (2008).  \CRANpkg{network}: A Package for Managing Relational Data in .
\emph{Journal of Statistical Software}, 24(2).
\doi{10.18637/jss.v024.i02}

Butts C (2015). \CRANpkg{network}: Classes for Relational Data. The Statnet
Project (\url{https://statnet.org}). R package version 1.12.0,
\url{https://cran.r-project.org/package=network}.

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
\CRANpkg{statnet} Tutorial.  \emph{Journal of Statistical Software}, 24(8).
\doi{10.18637/jss.v024.i08}

Goodreau SM, Kitts J, Morris M (2008b).  Birds of a Feather, or Friend of a
Friend? Using Exponential Random Graph Models to Investigate Adolescent
Social Networks.  \emph{Demography}, 45, in press.

Handcock, M. S. (2003) Assessing Degeneracy in Statistical Models of Social
Networks, Working Paper #39, Center for Statistics and the Social Sciences,
University of Washington.
\url{https://csss.uw.edu/research/working-papers/assessing-degeneracy-statistical-models-social-networks}

Handcock MS (2003b).  \CRANpkg{degreenet}: Models for Skewed Count Distributions
Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.0,
\url{https://statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
\CRANpkg{statnet}: Software Tools for the Statistical Modeling of Network Data.
Statnet Project, Seattle, WA.  Version 3, \url{https://statnet.org}.

Hunter, D. R. and Handcock, M. S. (2006) Inference in curved exponential
family models for networks, \emph{Journal of Computational and Graphical
Statistics}, 15: 565-583

Krivitsky PN, Handcock MS (2007).  \CRANpkg{latentnet}: Latent position and
cluster models for statistical networks.  Seattle, WA.  Version 2,
\url{https://statnet.org}.

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128.
\doi{10.1214/12-EJS696}

Morris M, Handcock MS, Hunter DR (2008).  Specification of
Exponential-Family Random Graph Models: Terms and Computational Aspects.
\emph{Journal of Statistical Software}, 24(4).
\doi{10.18637/jss.v024.i04}

Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
networks. \emph{Journal of the American Statistical Association}, 85,
204-212.
}
\seealso{
\code{\link{ergmTerm}}, \code{\link{ergmConstraint}}, \code{\link{ergmReference}},
\code{\link{ergmHint}}, and \code{\link{ergmProposal}} for indices of model
specification and estimation components visible to the \CRANpkg{ergm}'s API at any given time.
}
\author{
\strong{Maintainer}: Pavel N. Krivitsky \email{pavel@statnet.org} (\href{https://orcid.org/0000-0002-9101-3362}{ORCID})

Authors:
\itemize{
  \item Mark S. Handcock \email{handcock@stat.ucla.edu}
  \item David R. Hunter \email{dhunter@stat.psu.edu}
  \item Carter T. Butts \email{buttsc@uci.edu}
  \item Steven M. Goodreau \email{goodreau@u.washington.edu}
  \item Martina Morris \email{morrism@u.washington.edu}
}

Other contributors:
\itemize{
  \item Li Wang \email{lxwang@gmail.com} [contributor]
  \item Kirk Li \email{kirkli@u.washington.edu} [contributor]
  \item Skye Bender-deMoll \email{skyebend@u.washington.edu} [contributor]
  \item Chad Klumb \email{cklumb@gmail.com} [contributor]
  \item Michał Bojanowski \email{michal2992@gmail.com} (\href{https://orcid.org/0000-0001-7503-852X}{ORCID}) [contributor]
  \item Ben Bolker \email{bbolker+lme4@gmail.com} [contributor]
  \item Christian Schmid \email{songhyo86@gmail.com} [contributor]
  \item Joyce Cheng \email{joyce.cheng@student.unsw.edu.au} [contributor]
  \item Arya Karami \email{a.karami@unsw.edu.au} [contributor]
  \item Adrien Le Guillou \email{git@aleguillou.org} (\href{https://orcid.org/0000-0002-4791-418X}{ORCID}) [contributor]
}

}
\keyword{models}
