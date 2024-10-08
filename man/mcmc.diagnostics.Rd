% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.diagnostics.ergm.R
\name{mcmc.diagnostics}
\alias{mcmc.diagnostics}
\alias{mcmc.diagnostics.default}
\alias{mcmc.diagnostics.ergm}
\title{Conduct MCMC diagnostics on a model fit}
\usage{
mcmc.diagnostics(object, ...)

\method{mcmc.diagnostics}{ergm}(
  object,
  center = TRUE,
  esteq = TRUE,
  vars.per.page = 3,
  which = c("plots", "texts", "summary", "autocorrelation", "crosscorrelation", "burnin"),
  compact = FALSE,
  ...
)
}
\arguments{
\item{object}{A model fit object to be diagnosed.}

\item{\dots}{Additional arguments, to be passed to plotting functions.}

\item{center}{Logical: If \code{TRUE}, center the samples on the observed
statistics.}

\item{esteq}{Logical: If \code{TRUE}, for statistics corresponding to
curved ERGM terms, summarize the curved statistics by their
negated estimating function values (evaluated at the MLE of any curved
parameters) (i.e., \eqn{\eta'_{I}(\hat{\theta})\cdot (g_{I}(Y)-g_{I}(y))}
for \eqn{I} being indices of the canonical parameters in
question), rather than the canonical (sufficient) vectors of the
curved statistics relative to the observed (\eqn{g_{I}(Y)-g_{I}(y)}).}

\item{vars.per.page}{Number of rows (one variable per row) per
plotting page.  Ignored if \CRANpkg{latticeExtra} package is not
installed.}

\item{which}{A character vector specifying which diagnostics to
plot and/or print. Defaults to all of the below if meaningful:
\describe{

\item{\code{"plots"}}{Traceplots and density plots of sample values for all statistic or estimating function elements.}

\item{\code{"texts"}}{Shorthand for the following text diagnostics.}

\item{\code{"summary"}}{Summary of network statistic or estimating function elements as produced by \code{\link[coda:summary.mcmc]{coda::summary.mcmc.list()}}.}

\item{\code{"autocorrelation"}}{Autocorrelation of each of the network statistic or estimating function elements.}

\item{\code{"crosscorrelation"}}{Cross-correlations between each pair of the network statistic or estimating function elements.}

\item{\code{"burnin"}}{Burn-in diagnostics, in particular, the Geweke test.}

} Partial matching is supported. (E.g., \code{which=c("auto","cross")}
will print autocorrelation and cross-correlations.)}

\item{compact}{Numeric: For diagnostics that print variables in
columns (e.g. correlations, hypothesis test p-values), try to
abbreviate variable names to this many characters and round the
numbers to \code{compact - 2} digits after the decimal point; 0 or
\code{FALSE} for no abbreviation.}
}
\description{
This function prints diagnistic information and creates simple diagnostic
plots for MCMC sampled statistics produced from a fit.
}
\details{
A pair of plots are produced for each statistic:a trace of the sampled
output statistic values on the left and density estimate for each variable
in the MCMC chain on the right.  Diagnostics printed to the console include
correlations and convergence diagnostics.

For \code{\link[=ergm]{ergm()}} specifically, recent changes in the
estimation algorithm mean that these plots can no longer be used
to ensure that the mean statistics from the model match the
observed network statistics. For that functionality, please use
the GOF command: \code{gof(object, GOF=~model)}.

In fact, an \code{\link[=ergm]{ergm()}} output object contains the sample of
statistics from the last MCMC run as element \verb{$sample}. If
missing data MLE is fit, the corresponding element is named
\verb{$sample.obs}. These are objects of \code{\link[coda:mcmc]{mcmc}} and can be used
directly in the \CRANpkg{coda} package to assess MCMC
convergence.

More information can be found by looking at the documentation of
\code{\link[=ergm]{ergm()}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{mcmc.diagnostics(ergm)}: 

}}
\examples{

\dontrun{
#
data(florentine)
#
# test the mcmc.diagnostics function
#
gest <- ergm(flomarriage ~ edges + kstar(2))
summary(gest)

#
# Plot the probabilities first
#
mcmc.diagnostics(gest)
#
# Use coda directly
#
library(coda)
#
plot(gest$sample, ask=FALSE)
#
# A full range of diagnostics is available
# using codamenu()
#
}

}
\references{
Raftery, A.E. and Lewis, S.M. (1995).  The number of iterations, convergence
diagnostics and generic Metropolis algorithms.  In Practical Markov Chain
Monte Carlo (W.R. Gilks, D.J. Spiegelhalter and S. Richardson, eds.).
London, U.K.: Chapman and Hall.
}
\seealso{
\code{\link[=ergm]{ergm()}}, \CRANpkg{network} package, \CRANpkg{coda} package,
\code{\link[=summary.ergm]{summary.ergm()}}
}
\keyword{models}
