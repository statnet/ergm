%  File man/ergm.mple.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.mple}
\alias{ergm.mple}
\alias{ergm.pl}

\title{
Find a maximizer to the psuedolikelihood function
}
\description{
The \code{ergm.mple} function finds a maximizer to the psuedolikelihood function (MPLE). It is the default method for finding the ERGM starting coefficient values. It is normally called internally the ergm process and not directly by the user. Generally \code{\link{ergmMPLE}} would be called by users instead.

\code{ergm.pl} is an even more internal workhorse function that prepares many of the components needed by \code{ergm.mple} for the regression rountines that are used to find the MPLE estimated ergm. It should not be called directly by the user. 
}
\usage{
ergm.mple(Clist, Clist.miss, m, init = NULL, MPLEtype = "glm", family = "binomial",
         maxMPLEsamplesize = 1e+06, save.glm = TRUE, theta1 = NULL, conddeg = NULL, 
         control = NULL, MHproposal = NULL, verbose = FALSE, ...)

ergm.pl(Clist, Clist.miss, m, theta.offset=NULL,
                    maxMPLEsamplesize=1e+6,
                    conddeg=NULL, control, MHproposal,
                    verbose=FALSE)
}                    

\arguments{
  \item{Clist}{
a list of parameters used for fitting and returned by \code{\link{ergm.Cprepare}}
}
  \item{Clist.miss}{
the corresponding 'Clist' for the network of missing
                  edges returned by \code{\link{ergm.design}}
}
  \item{m}{
the model, as returned by \code{\link{ergm.getmodel}}
}
  \item{init}{
a vector a vector of initial theta coefficients
}
  \item{MPLEtype}{
the method for MPL estimation as "penalized", "glm" or  "logitreg"; default="glm"    
}
  \item{family}{
the family to use in the R native routine \code{\link{glm}}; only
                    applicable if "glm" is the 'MPLEtype'; default="binomial"
}
  \item{maxMPLEsamplesize}{
the sample size to use for endogenous sampling in the psuedo-likelihood computation; default=1e6
}
  \item{save.glm}{
whether the mple fit and the null mple fit should be
                      returned (T or F); if false, NULL is returned for both;
                    default==TRUE
}
  \item{theta1}{
the independence theta; if specified and non-NULL, this is
                     ignored except to return its value in the returned ergm;
                     default=NULL, in which case 'theta1' is computed
}
  \item{conddeg}{
an indicator of whether the MPLE should be conditional
                      on degree; non-NULL values indicate yes, NULL no;
                     default=NULL.
}
  \item{control}{
a list of MCMC related parameters; recognized components
                      include:
         samplesize : the number of networks to sample
         Clist.miss : see 'Clist.miss' above; some of the code uses this Clist.miss,
}
  \item{MHproposal}{
an MHproposal object, as returned by \code{\link{MHproposal}}
}
  \item{verbose}{
whether this and the C routines should be verbose (T or F);
                     default=FALSE
}

 \item{theta.offset}{a logical vector specifying which of the model
                    coefficients are offset, i.e. fixed}
  \item{\dots}{
additional parameters passed from within; all will be
                     ignored
}
}
\details{
According to Hunter et al. (2008):
"The maximizer of the pseudolikelihood may thus easily be found (at least in principle) by using logistic regression as a computational device."
In order for this to work, the predictors of the logistic regression model
must be calculated.  These are the change statistics as described
in Section 3.2 of Hunter et al. (2008), put into matrix form
so that each pair of nodes is one row whose values are the vector of change
statistics for that node pair.  The ergm.pl function computes these change
statistics and the ergm.mple function implements the logistic regression using
R's glm function.  Generally, neither ergm.mple nor ergm.pl should be called
by users if the logistic regression output is desired; instead, use the \code{\link{ergmMPLE}} function.

In the case where the ERGM is a dyadic independence model, the MPLE is the same as the MLE.  However, in general this is not the case and, as van Duijn et al. (2009) warn, the statistical properties of MPLEs in general are somewhat mysterious.

MPLE values are used even in the case of dyadic dependence models as starting points for the MCMC algorithm.
}
\value{
\code{ergm.mple} returns an ergm object as a list containing several items; for details see the return list in the \code{\link{ergm}}

\code{ergm.pl} returns a list containing:
\itemize{
\item{  xmat     : the compressed and possibly sampled matrix of change
             statistics}
\item{  zy       : the corresponding vector of responses, i.e. tie values}
\item{  foffset  : ??}
\item{  wend     : the vector of weights for 'xmat' and 'zy'}
\item{  numobs   : the number of dyads }
\item{  xmat.full: the 'xmat' before sampling; if no sampling is needed, this
               is NULL}
\item{  zy.full  : the 'zy' before  sampling; if no sampling is needed, this
                is NULL}
\item{  foffset.full     : ??}
\item{  theta.offset     : a numeric vector whose ith entry tells whether the
                        the ith curved coefficient?? was offset/fixed; -Inf
                       implies the coefficient was fixed, 0 otherwise; if
                        the model hasn't any curved terms, the first entry
                       of this vector is one of
                           log(Clist$nedges/(Clist$ndyads-Clist$nedges))
                           log(1/(Clist$ndyads-1))
                        depending on 'Clist$nedges'}
\item{   maxMPLEsamplesize: the 'maxMPLEsamplesize' inputted to \code{ergm.pl}}
}
}
\references{
Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris and Martina (2008). 
"ergm: A Package to Fit, Simulate and Diagnose Exponential-Family Models for Networks." _Journal of Statistical Software_, *24*(3), pp. 1-29. \url{http://www.jstatsoft.org/article/view/v024i03}

van  Duijn  MAJ,  Gile  K,  Handcock  MS  (2009).  "Comparison  of  Maximum  Pseudo  Likelihood and Maximum Likelihood Estimation of Exponential Family Random Graph Models."
_Social Networks_, *31*, pp. 52-62.
}

\seealso{
\code{\link{ergmMPLE}}, \code{\link{ergm}},\code{\link{control.ergm}}
}


