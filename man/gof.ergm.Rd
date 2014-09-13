%  File man/gof.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{gof}
\alias{gof.ergm}
\alias{gof.default}
\alias{gof.formula}
\alias{gof.ergm}
\alias{gof}
\title{Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph Model}
\description{
  \code{\link{gof}} calculates \eqn{p}-values for geodesic
  distance, degree, and reachability summaries to
  diagnose the goodness-of-fit of exponential family random graph
  models.  See \code{\link{ergm}} for more information on these models.
}
\usage{
\method{gof}{default}(object,\dots)
\method{gof}{formula}(object, 
                      \dots, 
                      coef=NULL,
                      GOF=NULL,
                      constraints=~.,
                      control=control.gof.formula(),
                      verbose=FALSE)
\method{gof}{ergm}(object, 
                   \dots,
                   coef=NULL, 
                   GOF=NULL,
                   constraints=NULL,
                   control=control.gof.ergm(),
                   verbose=FALSE)
}
\arguments{
  \item{object}{an \R object. Either a formula or an \code{\link{ergm}} object. 
    See documentation for \code{\link{ergm}}.
  }
  \item{\dots}{Additional arguments, to be passed to lower-level functions
    in the future.
  }
  \item{coef}{
    When given either
    a formula or an object of class ergm, \code{coef} are the parameters
    from which the sample is drawn. By default set to a vector of 0.
  }
  \item{GOF}{formula; an \R formula object, of the form
    \code{~ <model terms>} specifying the
    statistics to use to diagnosis the goodness-of-fit of the model. 
    They do not need to be in the model formula specified in
    \code{formula}, and typically are not.
    Currently supported terms are the degree distribution (\dQuote{degree}
    for undirected graphs, or \dQuote{idegree} and/or \dQuote{odegree} for 
    directed graphs), geodesic distances (\dQuote{distance}), shared partner 
    distributions (\dQuote{espartners} and \dQuote{dspartners}), the triad
    census (\dQuote{triadcensus}), and the terms of the original model 
    (\dQuote{model}). The default formula for undirected networks is
    \code{~ degree + espartners + distance}, and the default formula for
    directed networks is \code{~ idegree + odegree + espartners + distance}.
  }
  \item{constraints}{A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    modeled. See the help for similarly-named argument in
    \code{\link{ergm}} for more information. For
    \code{gof.formula}, defaults to unconstrained. For \code{gof.ergm},
    defaults to the constraints with which \code{object} was fitted.
  }
  \item{control}{A list to control parameters, constructed using
    \code{\link{control.gof.formula}} or \code{\link{control.gof.ergm}}
    (which have different defaults).
  }
  \item{verbose}{Provide verbose information on the progress of the
    simulation.
  }
}

\details{
  A sample of graphs is randomly drawn from the specified model.
  The first argument is typically 
  the output of a call to \code{\link{ergm}} and the model
  used for that call is the one fit.

  A plot of the summary measures is plotted.
  More information can be found by looking at the documentation of
  \code{\link{ergm}}.
  
  For \code{gof.ergm} and \code{gof.formula}, default behavior depends on
  the directedness of the network involved; if undirected then degree,
  espartners, and distance are used as default properties to examine.  If
  the network in question is directed, \dQuote{degree} in the above is replaced
  by idegree and odegree.
}
\value{\code{\link{gof}}, \code{\link{gof.ergm}}, and \code{\link{gof.formula}}
         return an object of class \code{gofobject}.
       This is a list of the tables of statistics and \eqn{p}-values.
       This is typically plotted using \code{\link{plot.gofobject}}.
      }

\seealso{ergm, network, simulate.ergm, summary.ergm, plot.gofobject}
\examples{
\donttest{
data(florentine)
gest <- ergm(flomarriage ~ edges + kstar(2))
gest
summary(gest)

# test the gof.ergm function
gofflo <- gof(gest)
gofflo
summary(gofflo)

# Plot all three on the same page
# with nice margins
par(mfrow=c(1,3))
par(oma=c(0.5,2,1,0.5))
plot(gofflo)

# And now the log-odds
plot(gofflo, plotlogodds=TRUE)

# Use the formula version of gof
gofflo2 <-gof(flomarriage ~ edges + kstar(2), coef=c(-1.6339, 0.0049))
plot(gofflo2)
}
}
\keyword{models}
