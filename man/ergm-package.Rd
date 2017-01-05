%  File man/ergm-package.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-package}
\alias{ergm-package}
\docType{package}
\title{
Fit, Simulate and Diagnose Exponential-Family Models for Networks
}
\description{
\code{\link[=ergm-package]{ergm}} is a collection of functions to plot,
fit, diagnose, and simulate from exponential-family random 
graph models (ERGMs). For a list of functions type: \code{help(package='ergm')}

For a complete list of the functions, use \code{library(help="ergm")} 
or read the rest of the manual. For a simple demonstration,
use \code{demo(packages="ergm")}.

When publishing results obtained using this package, please cite the
original authors as described in \code{citation(package="ergm")}.

All programs derived from this package must cite it.
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

This manual introduces software tools for the representation,
visualization, and analysis of network data that address each of these
previous shortcomings.  The package relies on the
\code{\link[network]{network}} package which allows networks to be
represented in \R. The \code{\link[=ergm-package]{ergm}} package
implements maximum likelihood estimates of ERGMs to be calculated using
Markov Chain Monte Carlo (via \code{\link{ergm}}). The package also
provides tools for simulating networks (via
\code{\link{simulate.ergm}}) and assessing model
goodness-of-fit (see \code{\link{mcmc.diagnostics}} and
\code{\link{gof.ergm}}).

A number of Statnet Project packages extend and enhance
\code{\link[=ergm-package]{ergm}}. These include
\code{\link[tergm:tergm-package]{tergm}} (Temporal ERGM), which provides
extensions for modeling evolution of networks over time;
\code{\link[ergm.count:ergm.count-package]{ergm.count}}, which
facilitates exponential family modeling for networks whose dyadic
measurements are counts; and
\code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}},
which allows users to implement their own ERGM terms.

For detailed information on how to download and install the software,
go to the \code{\link[=ergm-package]{ergm}} website:
\url{statnet.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Mark S. Handcock \email{handcock@stat.ucla.edu},\cr
David R. Hunter \email{dhunter@stat.psu.edu},\cr
Carter T. Butts \email{buttsc@uci.edu},\cr
Steven M. Goodreau \email{goodreau@u.washington.edu},\cr
Pavel N. Krivitsky \email{krivitsky@stat.psu.edu}, and\cr
Martina Morris \email{morrism@u.washington.edu}

Maintainer: Pavel N. Krivitsky \email{krivitsky@stat.psu.edu}
}
\references{
Admiraal R, Handcock MS (2007).
 \pkg{networksis}: Simulate bipartite graphs with fixed
  marginals through sequential importance sampling.
 Statnet Project, Seattle, WA.
 Version 1, \url{statnet.org}.

Bender-deMoll S, Morris M, Moody J (2008).
 Prototype Packages for Managing and Animating Longitudinal
  Network Data: \pkg{dynamicnetwork} and \pkg{rSoNIA}.
 \emph{Journal of Statistical Software}, 24(7).
 \url{http://www.jstatsoft.org/v24/i07/}.

Besag, J., 1974, Spatial interaction and the statistical analysis
of lattice systems (with discussion), Journal of the Royal Statistical
Society, B, 36, 192-236.

Boer P, Huisman M, Snijders T, Zeggelink E (2003).
 StOCNET: an open software system for the advanced statistical
  analysis of social networks.
 Groningen: ProGAMMA / ICS, version 1.4 edition.

% commented because server is down so paper can't be found
%Butts CT (2006).
% \pkg{netperm}: Permutation Models for Relational Data.
% Version 0.2, \url{http://erzuli.ss.uci.edu/R.stuff}.

Butts CT (2007).
 \pkg{sna}: Tools for Social Network Analysis.
 R package version 2.3-2. \url{https://cran.r-project.org/package=sna}

Butts CT (2008).
 \pkg{network}: A Package for Managing Relational Data in \R.
 \emph{Journal of Statistical Software}, 24(2).
 \url{http://www.jstatsoft.org/v24/i02/}.

Butts C (2015). \pkg{network}: Classes for Relational Data. The Statnet Project (\url{http://www.statnet.org}). R package version 1.12.0, \url{https://cran.r-project.org/package=network}.

Frank, O., and Strauss, D.(1986). Markov graphs. Journal of the American
Statistical Association, 81, 832-842. 

Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
 A \pkg{statnet} Tutorial.
 \emph{Journal of Statistical Software}, 24(8).
 \url{http://www.jstatsoft.org/v24/i08/}.

Goodreau SM, Kitts J, Morris M (2008b).
 Birds of a Feather, or Friend of a Friend? Using Exponential
  Random Graph Models to Investigate Adolescent Social Networks.
 \emph{Demography}, 45, in press.

Handcock, M. S. (2003)
    \emph{Assessing Degeneracy in Statistical Models of Social Networks},
    Working Paper \#39, 
Center for Statistics and the Social Sciences,
University of Washington.
\url{www.csss.washington.edu/Papers/wp39.pdf}

Handcock MS (2003b).
 \pkg{degreenet}: Models for Skewed Count Distributions Relevant
  to Networks.
 Statnet Project, Seattle, WA.
 Version 1.0, \url{statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003a).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 Statnet Project, Seattle, WA.
 Version 3, \url{statnet.org}.

Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
 \pkg{statnet}: Software Tools for the Statistical Modeling of
  Network Data.
 Statnet Project, Seattle, WA.
 Version 3, \url{statnet.org}.

Hunter, D. R. and Handcock, M. S. (2006)
   Inference in curved exponential family models for networks,
   \emph{Journal of Computational and Graphical Statistics},
   15: 565-583
   

Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
 \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
 \emph{Journal of Statistical Software}, 24(3).
 \url{http://www.jstatsoft.org/v24/i03/}.

Krivitsky PN, Handcock MS (2007).
 \pkg{latentnet}: Latent position and cluster models for
  statistical networks.
 Seattle, WA.
 Version 2, \url{http://statnet.org}.

Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
1100-1128. \doi{10.1214/12-EJS696}

 
Morris M, Handcock MS, Hunter DR (2008).
 Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
 \emph{Journal of Statistical Software}, 24(4).
 \url{http://www.jstatsoft.org/v24/i04/}.

Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
networks \emph{Journal of the American Statistical Association}, 85, 204-212. 
}
\keyword{ package }
\keyword{ models }
