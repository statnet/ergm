%  File inst/NEWS.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.
\section{Changes in version 3.1.2}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{\link{ergmMPLE}} can now output an array of change
      scores associated with each dyad. See the help page for more
      details.
    } 
  }
  
  \subsection{BUG FIXES}{
    \itemize{
      \item Some Metropolis-Hastings proposal functions would sometimes
      return incorrect acceptance probabilities when combined with the
      \code{\link[=ergm-constraints]{bd}} constraint. This has been
      fixed.

      \item In \code{\link{simulate.ergm}} and
      \code{\link{simulate.formula}}, some combinations of
      \code{control$parallel}, \code{statsonly}, and \code{sequential}
      arguments would cause an error. This has been fixed. Note that the
      current implementation of parallel simulation when
      \code{statsonly=FALSE} or \code{sequential=FALSE} is very
      inefficient.

      \item A namespace bug in the likelihod calculation code, reported
      by Matthew Hoover, has been fixed.

      \item MPLE for some degree-based constrained has been reenabled.
      
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item Argument \code{fitmodel} to \code{\link{ergmMPLE}} has been
      deprecated in favor of \code{output}.
    }
  }

  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      
      \item Default \code{SAN.burnin} control parameter has been
      increated to 100,000.
      
    }
  }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }

}

\section{Changes in version 3.1.1}{
 
  \subsection{BUG FIXES}{
    \itemize{
      \item Long-standing bugs in MCMC standard error calculation have
      been fixed. In particular, MCMC standard errors are now calculated
      for offset models, and they are no longer overestimated (by too
      much).
      
      \item \code{\link{san}} now behaves sensibly when the network exceeds the
      number of edges allocated.

      \item \code{\link{degrange}} terms now initialize correctly for their
      respective network types.
      
      \item Some spurious warnings are no longer printed.
      
      \item A variety of minor errors has been fixed.
    
    }
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link[=ergm-package]{ergm}}'s implementation and tests now work with
      \code{\link[=network-package]{network}} version 1.8. (In
      particular, \code{bipartite=0} is now treated as a network with no
      actors.)
    }
  }
}

\section{Changes in version 3.1}{
 
  \subsection{NEW FEATURES}{
    \itemize{ \item Infrastructure
      has been added for fitting valued ERGMs of Krivitsky (2012). The
      \dQuote{user interface} is subject to change, however. (The
      implementation of reference measures for counts can be found in \R package
      \code{\link{ergm.count}})

      \item \code{\link[=ergm-package]{ergm}} itself now includes
      rudimentary implementations for discrtete uniform and continuous
      uniform reference measures. See \code{help('ergm-references')}.

      \item A large number of ERGM terms for counts have been
      implemented and documented. See \code{help('ergm-terms')}.

      \item A number of new binary ERGM terms have been
      implemented as well. Particularly, \code{\link{degrange}},
      \code{\link{odegrange}}, \code{\link{idegrange}},
      \code{\link{b1degrange}}, and \code{\link{b2degrange}}, which
      generalize \code{\link{degree}} and \code{\link{concurrent}} by
      allowing the user to specify contiguous ranges of degree values.

      \item A rudimentary \dQuote{degeneracy} heuristic has been
      implemented and enabled by default, immediately stopping MC MLE
      optimization if simulating from a parameter configuration of
      interest encounters produces networks with too many edges. See
      \code{MCMLE.density.guard} argument in \code{\link{control.ergm}}
      documentation for details. Note that this is intended mainly to
      keep the fitting from bogging down on bad models, and it should
      not be used as a formal degeneracy diagnostic.

      \item \code{\link[=ergm-package]{ergm}} now keeps track of of when
      it is and when it's not possible to evaluate the likelihood of the
      null model (a model with no terms), and takes into account whether
      a model is constrained, has missing dyads, etc.. This can also
      improve the MPLE.
      
      A new generic, \code{\link{logLikNull}} has been added that
      returns the null model likelihood.

      For models where this is not possible, the
      null likelihood (and deviance, etc.) is defined to be 0.
      
      \item A new constraint, \code{blockdiag} has been added, for
      modeling networks with only within-block ties allowed. This is
      useful for, e.g., fitting an ERGM to several networks jointly. See
      \code{? 'ergm-constraints'} for available constraints.
      
    } 
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item The implementation of dynamic network models has been moved to a
      new package, \code{\link[=tergm-package]{tergm}}.

      \item Sample space constraints \code{\link{degrees}},
      \code{\link{idegrees}}, \code{\link{odegrees}},
      \code{\link{b1degrees}}, \code{\link{b2degrees}} have
      been reimplemented or fixed.

      \item Certain functions are now only exported if they haven't
      already be defined (e.g., by
      \code{\link[=latentnet-package]{latentnet}}).

      \item Documentation and citation information have been  cleaned up
      and reorganized. In particular, see \code{? 'ergm-constraints'},
      \code{? 'ergm-terms'} for an article on terms, and \code{?
      'ergm-references'} for reference measures implemented.
    }
  }

  \subsection{C-LEVEL FACILITIES}{
    \itemize{
      
      \item \code{\link[=ergm-package]{ergm}} now exports header and
      stub files for the \code{LinkingTo:} API, and package
      \code{ergm.userterms} now uses it. Please see
      \code{inst/include/README} before using it, however.
      
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item Dynamic network modeling code has been moved to
      \code{\link[=tergm-package]{tergm}} package.

      \item Function \code{statnet.edit} has been removed.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Some bugs introduced by the developers relying on a certain
      behavior of the \code{\link{.C}} API which had changed have been
      fixed. 

      \item Some undocumented terms have been fixed and documented,
      including \code{\link{odegreepopularity}},
      \code{\link{idegreepopularity}}, \code{\link{degreepopularity}},
      and \code{\link{smalldiff}}.
    
    }
  }
}

\section{Changes in version 3.0}{
 
  \subsection{NEW FEATURES}{
    \itemize{
      \item Facilities for fitting and simulating Separable Temporal
      ERGMs (STERGMs) of Krivitsky (2009) and, Krivitsky and Handcock
      (2010) have been added. See \code{\link{stergm}} and
      \code{\link[=simulate.network]{simulate}}.
      
      \item \code{\link{ergm}} MCMLE will now test if the optimization
      has converged to the MLE and stop if it detects convergence.

      \item \code{\link{ergm}} Now checks if a constraint specified in
      the \code{constraints} argument fixes any statistics in the
      model (e.g., \code{ergm(y~edges, constraints=~edges)}). It issues a
      warning, and fixes their parameters at 0.

      \item \code{\link{simulate}} functions now take an argument
      \code{monitor}, a one-sided formula for additional statistics to
      be returned. This allows one to easily simulate from an
      \code{\link{ergm}} or \code{\link{stergm}} fit statistics that
      were not in the fitted model.

      \item \code{\link{ergm}} can now optionally detect when a burn-in
      is insufficient and retry it until it is sufficient. See
      \code{\link{control.ergm}} argument \code{MCMC.burnin.retries} for
      details.

      \item \code{\link{ergm}} MCMLE now stores \dQuote{history} of
      \eqn{\theta} and simulated statistics.

      \item \code{\link{ergm}} extreme statistics detection (if
      \code{control.ergm(drop=TRUE)}) can now detect both extremely high
      and extremely low values, setting the coefficients to +Inf and
      -Inf respectively. Direction of \dQuote{dropping} is returned.

      \item An initial method \code{"zeros"} has been, to use a vector
      of zeros as starting parameters.

      \item A number of new ERGM terms have been added.
      
    }
  }
  
  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{ergm}}, \code{\link[=simulate.ergm]{simulate}},
      and other functions now enforce a conceptual separation between
      \emph{what} to fit, which is controlled by arguments to the
      function and \emph{how} to fit it, which is controlled by the
      \code{control} lists, generated by \code{\link{control.ergm}} and
      \code{\link{control.simulate}}, respectively. (Old arguments for
      \code{ergm} and \code{simulate} still work with a warning, and
      they may be removed later.)
      
      \item The data structure returned by
      \code{\link[=simulate.ergm]{simulate}} is now a list of networks
      with \code{\link{attr}}-style attributes describing how it was
      generated. The list has print and summary methods, and it can be
      used in the left-hand side of a
      \code{\link[=summary.formula]{summary}} formula.

      \item Matrices of sample statistics returned by \code{\link{ergm}}
      are \code{\link[coda]{mcmc}} objects (or
      \code{\link[coda]{mcmc.list}} objects, for parallel computing).
      
      \item Arguments \code{theta} have been renamed to
      \code{coef}. Arguments \code{meanstats} have been renamed to
      \code{target.stats}. Many control parameters have also been
      renamed and reorganized. (Old parameters still work with a
      warning, and they may be removed later.)
      
      \item A number of arguments to \code{\link{control.ergm}},
      \code{\link{control.simulate}}, and others, have been renamed.

      \item \code{\link{ergm}} now uses bridge sampling evaluate the
      log-likelihood by default. (Before, the default was to not
      evaluate it.)

      \item \code{\link[=mcmc.diagnostics]{mcmc.diagnostics}} has been
      rewritten:
      \itemize{
	\item Raftery-Lewis diagnostic has been  replaced with several
	others, including tests of (non)convergence.
	\item For fits with missing data, diagnostics for both the
	constrained and the unconstrained sample are printed and plotted.
	\item Plots now use \code{\link{lattice}} graphics, if available.
      }

      \item Handling of networks with many edges has been improved
      somewhat.

      \item \code{\link{gof.ergm}} and \code{\link{gof.formula}} now set
      the statistics used to diagnose the fit adaptively, depending on
      the directedness of the network being modeled.
    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item The following control parameters have been removed or
      renamed:
      \code{maxedges}, \code{prop.weights}, \code{prop.args},
      \code{packagename}, \code{summarizestats},    \code{nr.maxit},
      \code{nr.reltol}, \code{maxNumDyadTypes},
      \code{steplength}, \code{initialfit}, \code{style},
      \code{obs.MCMCsamplesize}, \code{obs.interval}, \code{obs.burnin},
      \code{compress}, \code{metric}, \code{force.mcmc},
      \code{adaptive.trustregion}, \code{adaptive.epsilon},
      \code{mcmc.precision}, \code{method}, \code{MPLEtype},
      \code{check.degeneracy}, \code{maxMPLEsamplesize},
      \code{MPLEsamplesize}, \code{phase1_n}, \code{nsubphases},
      \code{RobMon.phase1n_base}, \code{RobMon.phase2n_base},
      \code{RobMon.phase2sub}, \code{RobMon.init_gain},
      \code{RobMon.phase3n}, \code{trustregion}, \code{stepMCMCsize},
      \code{steppingmaxit}, \code{gridsize}, \code{sequential},
      \code{returnMCMCstats}, \code{calc.mcmc.se}, \code{hessian},
      \code{prop.weights}, \code{prop.args}, and \code{packagenames}.
    
      \item The following arguments to functions have been moved to the
      corresponding control functions, removed or renamed:
      \code{burnin}, \code{MCMCsamplesize}, \code{interval},
      \code{maxit}, \code{seed}, and \code{theta0}.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed some errors in log-likelihood calculation, especially
      in models with missing data and/or offset terms.

      \item Parallel processing has been made more robust. In
      particular, \code{\link[=ergm-package]{ergm}} now checks to make
      sure that the version it finds on slave nodes is the same as the
      version on the master. Also, although \code{\link{ergm}} will
      start its own cluster as needed, it will also detect if an MPI
      cluster has already been started and use that.

      \item Several bugs relating to handling models with with offset
      terms have been fixed.

      \item The bounded degree (\code{bd}) ERGM would sometimes tabulate
      degrees incorrectly. This has been fixed.

      \item Miscellaneous robustification and bug fixes.
    }
  }
}

\section{Changes in version 2.4}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The \dQuote{news} feature has been implemented. Type
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

      \item The following new models terms have been added:
      \code{\link{degcrossprod}} and \code{\link{degcor}}.

      \item Tests for change statistics have been added. 
  
      \item A new dataset has been added:  \code{\link{kapferer}}, the well-known
      tailor shop network dataset of Kapferer (1972).
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item The package has been \dQuote{standardized} to represent
      directed ties as \code{(tail,head)}.   
    }
  }
  
  \subsection{BUG FIXES}{
    \itemize{
      \item Log-likelihood calculation for dyad-dependent models has
    been fixed.

      \item Some errors in estimating parameters in curved exponential
      family models have been fixed.

      \item \code{\link[=gof.ergm]{gof}} now automatically selects the
      right goodness-of-fit statistics for the network type (directed
      vs. undirected) given.
      
      \item Bug in is.inCH in 2.4-2 fixed in 2.4-3.
    }
  }

    \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item \code{as.directed} function has been removed.

      \item All model terms included in the package have been
      reimplemented using the new \code{InitErgmTerm} framework.
    }
  } 
}
