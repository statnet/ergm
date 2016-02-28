%  File man/anova.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{approx.hotelling.diff.test}
\title{Approximate Hotelling T^2-Test for One Sample Means}
\alias{approx.hotelling.diff.test}
%\alias{anova.ergm.null}

\docType{methods}

\description{
A multivariate hypothesis test of difference between a vector of sample means and a vector of true means (\code{mu}). The null hypothesis assumes no differences for each pair of sample mean and true mean, the alternative hypothesis assumes not all pairs are equal. Smaller p-value would reject the null hypothesis.  The distribution under the null hypothesis is an approximation to Hotelling's T^2 distribution, hence the test is an approximated test. 
}

\section{Hotelling's T^2 distribution}{
A multivariate method that is the multivariate counterpart of Student's t and which also forms the basis for certain multivariate control charts is based on Hotelling's T2 distribution, which was introduced by Hotelling (1947).}
 
\keyword{models}


\references{

Hotelling, H. (1947). Multivariate Quality Control. In C. Eisenhart, M. W. Hastay, and W. A. Wallis, eds. Techniques of Statistical Analysis. New York: McGraw-Hill.
}
