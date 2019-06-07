%  File man/Getting.Started.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2003-2019 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{Getting.Started}
\alias{Getting.Started}
\docType{package}
\title{
Getting Started with "ergm": Fit, simulate and diagnose exponential-family models for networks
}
\description{
\code{\link{ergm}} is a collection of functions to plot, fit, diagnose, and simulate from random 
graph models. For a list of functions type: help(package='ergm')

For a complete list of the functions, use \code{library(help="ergm")} 
or read the rest of the manual. For a simple demonstration,
use \code{demo(packages="ergm")}.

When publishing results obtained using this package the original authors
are to be cited as given in \code{citation("ergm")}:

Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau,
and Martina Morris. 2003 
\emph{ergm: Fit, simulate and diagnose exponential-family models for networks}
\cr
\url{https://statnet.org}.

All published work derived from this package must cite it.
For complete citation information, use\cr
\code{citation(package="ergm")}.
}
\details{
Recent advances in the statistical modeling of random networks have had an
impact on the empirical study of social networks. Statistical exponential
family models (Strauss and Ikeda 1990) are a generalization of the Markov
random network models introduced by Frank and Strauss (1986), which in turn
derived from developments in spatial statistics (Besag, 1974). These models
recognize the complex dependencies within relational data structures.  To
date, the use of stochastic network models for networks has been limited by
three interrelated factors: the complexity of realistic models, the lack of
simulation tools for inference and validation, and a poor understanding of the
inferential properties of nontrivial models.

This manual introduces software tools for the representation, visualization,
and analysis of network data that address each of these previous shortcomings.
The package relies on the \code{\link[network]{network}} package which allows networks to be
represented in \R. The \code{\link{ergm}} package allows maximum likelihood estimates of
exponential random network models to be calculated using Markov Chain Monte
Carlo.  The package also provides tools for plotting networks, simulating
networks and assessing model goodness-of-fit.

For detailed information on how to download and install the software,
go to the \code{\link{ergm}} website:
\url{https://statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Mark S. Handcock \email{handcock@stat.ucla.edu},\cr
David R. Hunter \email{dhunter@stat.psu.edu},\cr
Carter T. Butts \email{buttsc@uci.edu},\cr
Steven M. Goodreau \email{goodreau@u.washington.edu},\cr
Pavel N. Krivitsky \email{krivitsky@stat.psu.edu}, and\cr
Martina Morris \email{morrism@u.washington.edu}

Maintainer: David R. Hunter \email{dhunter@stat.psu.edu}
}
\references{
Admiraal R, Handcock MS (2007).
 \pkg{networksis}: Simulate bipartite graphs with fixed
  marginals through sequential importance sampling.
 Statnet Project, Seattle, WA.
 Version 1, \url{https://statnet.org}.

Bender-deMoll S, Morris M, Moody J (2008).
 Prototype Packages for Managing and Animating Longitudinal
  Network Data: \pkg{dynamicnetwork} and \pkg{rSoNIA}.
 \emph{Journal of Statistical Software}, 24(7).
 \url{https://www.jstatsoft.org/v24/i07/}.

Besag, J., 1974, Spatial interaction and the statistical analysis
of lattice systems (with discussion), \emph{Journal of the Royal Statistical
Society, B}, 36, 192-236.

Boer P, Huisman M, Snijders T, Zeggelink E (2003).
 StOCNET: an open software system for the advanced statistical
  analysis of social networks.
 Groningen: ProGAMMA / ICS, version 1.4 edition.

%Butts CT (2006).
% \pkg{netperm}: Permutation Models for Relational Data.
% Version 0.2, \url{https://erzuli.ss.uci.edu/R.stuff}.

Butts CT (2007).
 \pkg{sna}: Tools for Social Network Analysis.
 R package version 2.3-2. \url{https://cran.r-project.org/package=sna}.

Butts CT (2008).
 \pkg{network}: A Package for Managing Relational Data in \R.
 \emph{Journal of Statistical Software}, 24(2).
 \url{https://www.jstatsoft.org/v24/i02/}.

Butts C (2015).
 \pkg{network}: The Statnet Project (https://statnet.org). R package version 1.12.0, \url{https://cran.r-project.org/package=network}.

Frank, O., and Strauss, D.(1986). Markov graphs. \emph{Journal of the American
Statistical Association}, 81, 832-842. 

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{https://www.jstatsoft.org/v24/i08/}.

Goodreau SM, Kitts J, Morris M (2008b).
 Birds of a Feather, or Friend of a Friend? Using Exponential
  Random Graph Models to Investigate Adolescent Social Networks.
 \emph{Demography}, 45, in press.

Handcock, M. S. (2003)
    Assessing Degeneracy in Statistical Models of Social Networks,
    Working Paper \#39, 
Center for Statistics and the Social Sciences,
University of Washington.
\url{https://www.csss.washington.edu/Papers/wp39.pdf}

Handcock MS (2003b).
 \pkg{degreenet}: Models for Skewed Count Distributions Relevant
  to Networks.
 Statnet Project, Seattle, WA.
 Version 1.0, \url{https://statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003a).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 Statnet Project, Seattle, WA.
 Version 2, \url{https://statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
 \pkg{statnet}: Software Tools for the Statistical Modeling of
  Network Data.
 Statnet Project, Seattle, WA.
 Version 2, \url{https://statnet.org}.

Hunter, D. R. and Handcock, M. S. (2006)
   Inference in curved exponential family models for networks,
   \emph{Journal of Computational and Graphical Statistics},
   15: 565-583.

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{https://www.jstatsoft.org/v24/i03/}.

Krivitsky PN, Handcock MS (2007).
 \pkg{latentnet}: Latent position and cluster models for
  statistical networks.
 Seattle, WA.
 Version 2, \url{https://statnet.org}.

Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{https://www.jstatsoft.org/v24/i04/}.

Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
networks. \emph{Journal of the American Statistical Association}, 85, 204-212. 
}
\keyword{ package }
\keyword{ models }
