%  File ergm/man/sampson.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{sampson}
\alias{sampson}
\alias{samplike}
\title{Cumulative network of positive affection within a monastery as a ``network" object}
\description{
Sampson (1969) recorded
the social interactions among a group of monks
while resident as an experimenter on vision, and collected
numerous sociometric rankings.  
During his stay, a political
``crisis in the cloister" resulted in the expulsion of four monks
(Nos. 2, 3, 17, and 18) and the voluntary departure of several
others - most immediately, Nos. 1, 7, 14, 15, and 16. (In the
end, only 5, 6, 9, and 11 remained).
Of particular interest
is the data on positive affect relations (``liking"), in which each
monk was asked if they had positive relations to each of the
other monks.

The data were gathered at three times to capture changes in group sentiment
over time.
They represent three time points in the period during which a
new cohort entered the monastery near the end of the study but
before the major conflict began. 
Each member ranked only his top three choices on ``liking."
(Some subjects offered tied ranks for their top four choices).
A tie from monk A to monk B exists if A
nominated B as one of his three best friends at that that time point.

\code{samplike} is the time-aggregated graph. 
It is the cumulative tie for ``liking" over the three
periods.  For this, a tie from monk A to monk B exists if A
nominated B as one of his three best friends at any of the
three time points.

The graph is stored as an \code{\link[network]{network}} objects.
It has three vertex attributes:

\describe{
\item{group}{Groups of novices as classified by Sampson: 
"Loyal", "Outcasts", and "Turks". There is also an 
interstitial group not represented here.}
\item{cloisterville}{An indicator if attendance the minor
seminary of ``Cloisterville" before coming to the monastery.}
\item{vertex.names}{The given names of the novices.}
}

This data set is standard in the social network analysis
literature, having been modeled by Holland and Leinhardt
(1981), Reitz (1982), Holland, Laskey and Leinhardt (1983),
and Fienberg, Meyer, and Wasserman (1981), Hoff, Raftery, and Handcock
(2002), etc.
This is only a small piece of the data collected by Sampson.

This dataset was updated for version 2.5 (March 2012) to add the \code{cloisterville} variable and refine the names. This information is from
de Nooy, Mrvar, and Batagelj (2005).
The original vertex names were: 
Romul_10, Bonaven_5, Ambrose_9, Berth_6, Peter_4, Louis_11, Victor_8, Winf_12, John_1, Greg_2, Hugh_14, Boni_15, Mark_7, Albert_16, Amand_13, Basil_3, Elias_17, Simp_18.
} \usage{
 data(sampson)
}
\source{
 Sampson, S.~F. (1968), \emph{A novitiate in a period of change: An
 experimental and case study of relationships,} 
 Unpublished Ph.D. dissertation, Department of Sociology, Cornell University.

 \url{http://vlado.fmf.uni-lj.si/pub/networks/data/esna/sampson.htm}
}
\seealso{florentine, network, plot.network, ergm}
\references{
 White, H.C., Boorman, S.A. and Breiger, R.L. (1976).
 \emph{Social structure from multiple networks. I. Blockmodels of roles and
 positions.}
 American Journal of Sociology, 81(4), 730-780.

 Wouter de Nooy, Andrej Mrvar, Vladimir Batagelj (2005) \emph{Exploratory Social Network Analysis with Pajek},
 Cambridge: Cambridge University Press
}
\keyword{datasets}
