%  File man/florentine.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{florentine}
\title{Florentine Family Marriage and Business Ties Data as a ``network" object}
\alias{florentine}
\description{
This is a data set of marriage and business ties among
Renaissance Florentine families. The data is originally from
Padgett (1994) via \code{UCINET} and stored as a \code{\link[network]{network}} object.

Breiger \& Pattison (1986), in their discussion of local
role analysis, use a subset of data on the social relations
among Renaissance Florentine families (person aggregates)
collected by John Padgett from historical documents. The two
relations are business ties (\code{\link{flobusiness}} - specifically, recorded
financial ties such as loans, credits and joint partnerships)
and marriage alliances (\code{\link{flomarriage}}).

As Breiger \& Pattison point out, the original data are
symmetrically coded. This is acceptable perhaps for marital
ties, but is unfortunate for the financial ties (which are
almost certainly directed). To remedy this, the financial
ties can be recoded as directed relations using some external
measure of power - for instance, a measure of wealth. Both graphs
provide vertex information on (1) \code{wealth} each family's net wealth in 1427
(in thousands of lira); (2) \code{priorates} the number of priorates (seats on
the civic council) held between 1282- 1344; and (3) \code{totalties} the total
number of business or marriage ties in the total dataset of
116 families (see Breiger \& Pattison (1986), p 239).

Substantively, the data include families who were locked
in a struggle for political control of the city of Florence
around 1430. Two factions were dominant in this struggle:
one revolved around the infamous Medicis (9), the other around
the powerful Strozzis (15).
}
\usage{
 data(florentine)
}
\source{Padgett, John F. 1994. Marriage and Elite Structure
in Renaissance Florence, 1282-1500. Paper delivered to the
Social Science History Association.}
\seealso{flo, network, plot.network, ergm}
\references{
Wasserman, S. and Faust, K. (1994) \emph{Social Network Analysis: Methods and
Applications}, Cambridge University Press, Cambridge, England.

Breiger R. and Pattison P. (1986). \emph{Cumulated social roles:
The duality of persons and their algebras}, Social Networks,
8, 215-256.
}
\keyword{datasets}
