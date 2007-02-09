\name{fauxhigh}
\title{Goodreau's ''Faux High School" as a ``network" object}
\alias{fauxhigh}
\description{
This data set represents a simulation of an in-school friendship 
network, based upon a model fit to data from one school community from the 
AddHealth Study, Wave I (Resnick et al., 1997).

The data set was constructed as follows:

A vector representing the sex of each student in the school was 
randomly re-ordered.  The same was done with the students' response to 
questions on race and grade.  These three attribute vectors were 
permuted independently. Missing values for each were randomly assigned.

The following \code{\link{ergm}} model was fit to the original data:

\code{
~ edges\
 + nodefactor(''Grade") + nodefactor(''Race") + nodefactor(''Sex")\
 + nodematch(''Grade",diff=T) + nodematch(''Race",diff=T)\
 + nodematch(''Sex",diff=F) + gwdegree(1.0,fixed=T)\
 + gwesp(1.0,fixed=T) + gwdsp(1.0,fixed=T)
}

The resulting model fit was then applied to a network with actors 
possessing the permuted attributes and with the same number of edges as in the 
original data.

The processes for handling missing data and defining the race attribute 
are described in Hunter, Goodreau \& Handcock (2005).
}

\seealso{network, plot.network, ergm}
\references{

Hunter D.R., Goodreau S.M. and Handcock M.S. (2005). \emph{Goodness of Fit 
of Social Network Models}, Working Paper \#46, 
Center for Statistics and the Social Sciences,
University of Washington.
\url{www.csss.washington.edu/Papers/wp47.pdf}

Resnick M.D., Bearman, P.S., Blum R.W. et al. (1997). \emph{Protecting adolescents
from harm. Findings from the National Longitudinal Study on Adolescent Health},
\emph{Journal of the American Medical Association}, 278: 823-32.

}

\keyword{datasets}
\section{Licenses and Citation}{
If the source of the data set does not specified otherwise,
this data set is protected by the Creative Commons License
\url{http://creativecommons.org/licenses/by-nc-nd/2.5/}.

When publishing results obtained using this data set the original authors should be cited.  In addition this package should be cited as:

Mark S. Handcock, David R. Hunter, Carter T. Butts, Steven M. Goodreau,
and Martina Morris. 2003 
\emph{statnet: An R package for the Statistical Modeling of Social Networks}
\cr
\url{http://www.csde.washington.edu/statnet}.
}
