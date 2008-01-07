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
are to be cited as:

Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau,
and Martina Morris. 2003 
\emph{ergm: An R package for the Statistical Modeling of Social Networks}
\cr
\url{http://statnetproject.org}.

All programs derived from this package must cite it.
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
\url{http://statnetproject.org}.
A tutorial, support newsgroup, references and links to further resources are provided there.
}
\author{
Mark S. Handcock \email{handcock@stat.washington.edu},\cr
David R. Hunter \email{dhunter@stat.psu.edu},\cr
Carter T. Butts \email{buttsc@uci.edu},\cr
Steven M. Goodreau \email{goodreau@u.washington.edu}, and\cr
Martina Morris \email{morrism@u.washington.edu}

Maintainer: Mark S. Handcock \email{handcock@stat.washington.edu}
}
\references{
Besag, J., 1974, Spatial interaction and the statistical analysis
of lattice systems (with discussion), Journal of the Royal Statistical
Society, B, 36, 192-236.

Frank, O., and Strauss, D.(1986). Markov graphs. Journal of the American
Statistical Association, 81, 832-842. 

Handcock, M. S., Hunter, D. R., Butts, C. T.,  Goodreau,
S. M., and Morris, M. (2003), 
\emph{statnet: Software tools for the Statistical Modeling of Network Data.},\cr
URL \url{http://statnetproject.org}

Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
networks. Journal of the American Statistical Association, 85, 204-212. 
}
\keyword{ package }
\keyword{ models }
