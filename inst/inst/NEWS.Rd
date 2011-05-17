\name{NEWS}
\title{NEWS file for the ergm package}
\encoding{UTF-8}

\section{Changes in version 2.4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Implemented "news" feature so users can type 
      \code{news(package="ergm")}.
      
      \item Added "ergm.allstats" function that calculates all possible network
      statistics for a given ergm and network.  The number of possible networks
      grows extremely fast, so this should only be used for small networks.  
      Related function \code{ergm.exact} calculates exact log-likelihood for 
      a particular canonical parameter and a particular ergm.

      \item Added "cutoff" argument to the following model terms:
      gwb1degree, gwb2degree, gwdegree, gwdsp, gwesp, gwidegree, gwnsp, gwodegree

      \item Implemented new "logLik.ergm" method to use existing logLik generic in R.

      \item Re-engineered calculation of loglikelihoods for MCMC-based ergm fits;
       added bridge sampler that can be used for these calculations in the case of dyad-dependent models
       (see ergm.bridge.llr and ergm.bridge.dindstart.llk functions)

      \item Added is.dyad.independent function to test whether an ERGM fit or formula is dyad-independent.
      
      \item Added fix.curved function to convert curved exponential family models to variants in which
            the curved parameters are fixed.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Identified and fixed errors in calculating loglikelihood values.

      \item Identified and fixed errors in estimating parameters in curved exponential family models.
    }
  }
}


