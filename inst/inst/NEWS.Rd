\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

\section{Changes in version 2.4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The "news" feature has been implemented. Type
  \code{news(package="ergm")} to see the version's updates.
      
      \item A function \code{\link{ergm.allstats}} was added to
      calculate all possible network statistics for a given ERGM and
      network. The number of possible networks grows on the order of \eqn{2^{n^2}}, so this should only be used for small networks.  Related function \code{\link{ergm.exact}} calculates exact log-likelihood for a particular canonical parameter and a particular ergm.

      \item A \code{cutoff} argument has been added to specify maximum shared partner
      count or degree to be evaluated to the following model terms:
      \code{\link{gwb1degree}}, \code{\link{gwb2degree}}, \code{\link{gwdegree}}, \code{\link{gwdsp}}, \code{\link{gwesp}}, \code{\link{gwidegree}}, \code{\link{gwnsp}},
      \code{\link{gwodegree}}. This significantly speeds up fitting
      curved ERGMs on large, sparse networks.

      \item The \code{\link{logLik}} method for
      \code{\link{ergm}} objects, with the attendant \code{\link{AIC}}
      and \code{\link{BIC}} methods, has been implemented.

      \item Calculation of loglikelihoods for MCMC-based ergm fits has
      been reengineered to use bridge sampling for dyad-dependent
      models. (See \code{\link{logLik.ergm}},
      \code{\link{ergm.bridge.llr}}, and
      \code{\link{ergm.bridge.dindstart.llk}}.)

      \item A generic function \code{\link{is.dyad.independent}} to test
      whether an ERGM fit or formula is dyad-independent has been added.
      
      \item A generic function \code{\link{fix.curved}} to convert
      curved exponential family models to variants in which the curved
      parameters are fixed has been added. See the function
      documentation for usage scenarios.

      \item A generic function \code{\link{enformulate}} function to
      work around a current limitation in the \code{\link{ergm}}
      function to substitute curved parameters from a model fit into the
      model formula, allowing them to be used as initial values (the
      \code{theta0} argument). See the function documentation for usage
      scenarios.

      \item The following new models terms have been added:  \code{\link{degcrossprod}} and \code{\link{degcor}}.
  
      \item For version 2.4-3, a new dataset has been added:  \code{\link{kapferer}}, the well-known
      tailor shop network dataset of Kapferer (1972).
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Log-likelihood calculation for dyad-dependent models has
    been fixed.

      \item Some errors in estimating parameters in curved exponential
    family models have been fixed.
      
      \item Bug in is.inCH in 2.4-2 fixed in 2.4-3.
    }
  }
}


