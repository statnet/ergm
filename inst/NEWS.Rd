%  File inst/NEWS.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.


\section{Changes in version 3.6.1}{

  %% \subsection{NEW FEATURES}{
  %%   \itemize{
  %%     \item{
  %%     }
  %%   }
  %% }
   
  
   \subsection{BUG FIXES}{
     \itemize{
       
       \item{The implementation of the bounded degree
	 (\code{\link[ergm-constraints]{bd}}) constraint had a bug that
	 distorted somewhat the distribution of networks simulated (and
	 therefore estimated) away from the boundary of the constraint.}
       
       \item{Vertex labels in the \code{\link{sampson}} and
	 \code{\link{samplike}} datasets have been corrected. (Thanks to
	 Dr. Linda Tan Siew Li for pointing this error out.)}

       \item{Use of the deprecated \code{\link{unix.time}} has been
       changed to \code{\link{system.time}}.}
       
       \item{Miscellaneous fixes to the documentation and robustifications.}

     }
   }

  
   %% \subsection{OTHER USER-VISIBLE CHANGES}{
   %%   \itemize{
   %%   }
   %% }

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


\section{Changes in version 3.6.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
	A new suite of geometrically-weighted family of terms has been
	added that allow the user to specify what directions of shared
	partnerships are of interest. The new terms are,
	\itemize{
	  \item{(Directed) geometrically weighted dyadwise shared partners (\code{\link{dgwdsp}})}
	  \item{(Directed) geometrically weighted edgewise shared partners (\code{\link{dgwesp}})}
	  \item{(Directed) geometrically weighted non-edgewise shared
	    partners (\code{\link{dgwnsp}})}
	}
	In addition, terms for specific shared partner counts
	(\code{\link{ddsp}}, \code{\link{desp}}, and \code{\link{dnsp}}) are also
	implemented. See \code{\link{ergm-terms}} for help.
      }
      
      \item{New terms (binary and valued) \code{\link{b1cov}} and \code{\link{b2cov}} have
	been implemented, analogous to \code{\link{nodecov}},
	\code{\link{nodeicov}}, and \code{\link{nodeocov}} but for
	bipatite networks.}
      
      \item{New valued terms \code{\link{atleast}},
	\code{\link{equalto}}, and \code{\link{smallerthan}} have
	been implemented. A number of terms that were binary-only before
	have been implemented for valued ERGMs.}

      \item{
	\code{ergm} package now works when loaded (e.g., using
	\code{\link{requireNamespace}}) but not attached (e.g., using
	\code{\link{library}}).
      }
    }
  }
   
  
   \subsection{BUG FIXES}{
     \itemize{
       \item{A number of implementation bugs in the valued ERGM terms
	 have been fixed.}
       \item{Parallel processing when using \code{ergm.userterms} and similar was
	 accidentaly broken in the previous release.}
       \item{Calculation of Hummel et. al. (2010) steplength for missing data MLE was
	 incorrect; a slower but more accurate algorithm is now used.}
       \item{Miscellaneous bugfixes and robustifications.}
       \item{\code{\link{gof.formula}} and others now work with missing
	 data.}
       \item{Miscelaneous cleanups and bug fixes.}       
     }
   }

  
   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
       \item{\code{\link{robust.inverse}} has been deprecated in favor
	 of \code{\link[MASS]{ginv}}.}
       \item{\code{\link{simulate.formula}} now does parallel processing
	 more efficiently.}
     }
   }
}


\section{Changes in version 3.5.1}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Added a new generic, \code{\link{is.curved}} that tests if a
      given ERGM or model is curved.

      \item Block-diagonal sampler now handles bipartite networks.

      \item Code for parallel calculation has been optimized to provide
      a bigger gain from multiple cores/processors/nodes.

      \item Bridge sampling to evaluate the log-likelihood at the MLE
      now prints out a brief progress report by default.
   }
  }
   
  
   \subsection{BUG FIXES}{
     \itemize{
       \item MPLE is no longer used if the model is curved and target
       statistics are passed, as that caused problems.

       \item MCMC standard error calculation, "dropping" of extreme
       terms, main MC MLE loop, and other places have been made more
       robust to weird combinations of estimating problems (e.g., curved
       terms + target statistics + offset terms + poor starting value).

       \item \code{\link[ergm-terms]{threepath}} has been deprecated in
       favor of \code{threetrail} to more accurately reflect what is
       being computed.

       \item The "naive" log-likelihood ratio metric now handles offsets
       more robustly.

       \item The null log-likelihood calculation prints messages instead
       of throwing warnings.
       
       \item The closed-form solution (see Hummel et al., 2010) for the "lognormal" metric had
       been inadvertantly disabled. It has been reenabled.

       \item Many documentation fixes and clarifications.
    }
   }

  
   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
       \item Skye Bender-deMoll has been added as a contributor to the
       package.
       
       \item Hummel et al.'s convex hull calculation and MCMC SE
       calculation have been robustified against highly correlated
       sample statistics by using PCA to rotate them into an orthogonal
       configuration.
       
       \item Some functions, particularly \code{ERRVL},
       \code{append.rhs.formula}, \code{vectors.namesmatch},
       \code{term.list.formula}, and \code{ergm.update.formula} have been
       moved to the \code{\link{statnet.common}} package.
       
       \item The \code{as.edgelist} family of functions have been moved
       to the \code{\link[network-package]{network}} package.
       
       \item \code{\link{mcmc.diagnostics.ergm}} now takes the
       \code{esteq=} argument. If \code{TRUE}, the diagnostics are
       performed on the simulated values of the estimating functions. If
       \code{FALSE}, they are performed on simulated values of canonical
       statistics. This obviates and replaces the the \code{curved=}
       argument.
       
       \item \code{\link{san}} when used with curved models now requires
       canonical statistics as targets.

     }
   }
}

\section{Changes in version 3.4}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item Added message to notify when ergm MCMLE does not converge after max iterations.
      \item Copied the combined Sampson monks dataset to ergm. See \code{?sampson}.
      \item ergm simulation, bridge, and logLik functions now "inherit"
    the parallel setting from the ergm object and generally pass it
    along when possible.
    \item A new estimate type and source for initial values, contrastive
    divergence, has been added. It is intended primarily for obtaining
    initial values for MCMLE for valued ERGMs and ERGMs with complex
    (e.g., degree and edge-count) constraints. For binary networks with
    dyad-independent constraints, MPLE still generally works better and
    remains the default.

    To obtain the contrastive divergence estimate, call
   \code{\link{ergm}} with \code{estimate="CD"}.
   }
  }
   
  
   \subsection{BUG FIXES}{
     \itemize{
      \item{Fixed issue where constraints was not appropriately compared to older version}
      \item{Fixed bug where the deviance of dyad-independence models was incorrectly computed/reported for all models}
      \item{Removed \code{Rglpk} package dependency; replaced it with \code{lpSolve} package dependency. See man page for \code{is.inCH} for details}
      \item{Fixed integer overflow bug in sample size reporting}
  
      \item{Fixed bug in saving target.stats vector. Offset terms in target.stats output are now NA instead of the SAN network stats. (issue #1049)}
      \item{Fixed ergm segfault when unknown node attributes are used in
      ergm formulas. Modified get.node.attr function so that it will
      check if specified attribute exists. (previous code only checked
      if network was bipartite). (issue #1140, reported by Philip
      Leifeld) }
           \item{Implemented solution to MCMLE steplength getting stuck at 0 for some models where initial networks are very far from target stats. If the Hummel step length is < control$MCMLE.steplength.min, stop with an informative error message, reporting that the the starting value is poor. If using target stats, try increasing SAN parameters, or your configuration may be impossible or unlikely. (This suggests a poor model for the observed data.) MCMLE.steplength.min=0.0001 by default. Estimation will stop with an error message when steplen is smaller than the min for 2 iterations. (issue #1089)}
       \item Replaced the integer function abs() by floating point function fabs() in MHproposals_degree.c.
       \item Various documentation corrections and fixes

    }
   }

  
   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
      \item{\code{\link{ergm}} with \code{estimate="MPLE"} will now
	refuse to estimate the MPLE for valued ERGMs, curved ERGMs, and ERGMs with
	certain complex constraints, as true MPLE
	for these is not yet implemented.}
       \item{ CRAN-requested code structure changes (generally not visable to user):
        \itemize{
        \item \code{require} and \code{library} calls changed to \code{requireNamespace} and namespace prefixes added to calls from other packages, 
        \item S3 methods explicitly registered in \code{NAMESPACE}.
        \item \code{packageDescription} is no longer used in package citation methods.
	\item Many formerly 'internal' functions have been documented in a
	semi-stable API for use by developers of other ergm-related packages.
  \item{ Many formerly undocumented 'internal' functions have been made internal by the \code{NAMESPACE}. As far as we know, this avoids functions in use by reverse-Depending packages.  These currently appear as commented items in the \code{NAMESPACE} file.
  }
  \item{ Many undocumented 'internal' functions that do not appear to be used have been deprecated and will likely be removed in the next release.  See \code{?ergm_deprecated} for a list.  If called, these functions will give a Deprecation warning and suggest an alternate function where appropriate. 
  
  }
        }
       }
     }
   }
}
\section{Changes in version 3.2.4}{

   \subsection{NEW FEATURES}{
     \itemize{
     \item{ergm MCMLE estimation, by default, uses the stepping algorithm by Hummel et al (2012) to determine step length, and terminates estimation when the step length = 1 for two consecutive iterations. The old method (from version 3.1) is still available by setting control parameter \code{MCMLE.termination="Hotelling"}. }
     \item{As a result of using Hummel stepping algorithm, the final MCMC sample average may be different from the observed statistics. This is reflected in \code{mcmc.diagnostics()}: these plots can no longer be used to ensure that the mean statistics from the model match the observed network statistics. For that functionality, please use the GOF command: \code{gof(object, GOF=~model)}.}
     \item{Default ergm MCMC parameters changed to \code{MCMC.interval=1024, MCMC.samplesize=1024, MCMC.burnin=MCMC.interval*16}. The final sample size in the Hummel stepping algorithm is boosted by a factor of \code{MCMLE.last.boost}}
     \item{Default values for a number of other control parameters have changed, see CHANGES below.}
     \item{Parallel functionality: ergm can take user-created clusters as the parallel control parameter. This is the recommended method for using ergm on a high-performance computing cluster. See \link{ergm-parallel}. Functionality now implemented via the \code{parallel} package, direct dependence on \code{snow} package is removed}
     \item{New helper function \code{search.ergmTerms()} function which through the \code{\link{ergm.terms}} help page and prints out a list of terms appropriate for the specified network's structural constraints, optionally restricting by additional categories and keyword matches.}
     \item{The estimated covariance matrix of the ergm estimate can be extracted using \code{vcov} function.}
     \item{New EXPERIMENTAL Effective Sample Size (ESS)  feature to adaptively determine MCMC length needed to get the specified effective size. Set \code{MCMC.effectiveSize} to non-NULL to enable this functionality; 50 is a reasonable value. Includes new \code{ergm.control} parameters
        \itemize{
          \item    \code{MCMC.effectiveSize=NULL}
          \item    \code{MCMC.effectiveSize.damp=10}
          \item    \code{MCMC.effectiveSize.maxruns=1000}
          \item    \code{MCMC.effectiveSize.base=1/2}
          \item    \code{MCMC.effectiveSize.points=5}
          \item    \code{MCMC.effectiveSize.order=0.1}
          \item    \code{MCMLE.MCMC.max.ESS.frac}
        }
      }
    \item SAN.control now includes an \code{MPLE.max.dyad.types} parameter which defaults to the value of the ergm.control parameter of the same name.
    \item { Addition of the \code{StdNormal} reference distribution. See \code{?ergm.references}}
    \item{ New ergm terms:
      \itemize{
       \item \code{b1mindegree, b2mindegree} Minimum degree for the designated mode in a bipartite network
       \item \code{b1nodematch, b2nodematch} Nodal attribute-based homophily effect for the designated mode in a bipartite network.
      }
    }
    \item { Additional arguments added to the \code{network.update} function to control which attributes and properties are updated.}
    \item {  new \code{esteq} argument for \code{simulate.formula} and \code{simulate.ergm}:  Default is FALSE. If TRUE, compute the sample estimating equations of an ERGM. If the model is linear, all non-offset statistics are passed. If the model is curved, the score estimating equations (3.1) by Hunter and Handcock (2006) are given instead.}
    \item{ New and expanded tests of the ergm functions for use with R CMD CHECK tools. }
    \item{ Addition of \code{is.durational} methods for models and formulas. These functions test whether the terms in an ERGM model or formula is duration dependent or not. If the formula or model does not include any terms that need information about the duration of existing ties, the ergm process can use more efficient internal data structures. This is mostly in support of \code{tergm} functionality.} 
    \item{New ergm constraint terms: \code{fixallbut(free.dyads)}, \code{fixedas(present,absent)}. The former one can specify the set of dyads (as an edgelist or a network object) to be fixed, the latter one can specify two sets of dyads (as edgelists or networks objects), to be fixed present and fixed absent respectively. See \code{?ergm.constraints}}
    \item{Potential speed/memory usage improvement: sampling algorithms in both R can C levels are optimized.}
    \item{ergm now imports the \code{network} package via NAMESPACE to avoid internal namespace collisions when used with \code{iGraph}}
  }
  }
     
   
   
   \subsection{CHANGES}{
    \itemize{
      \item {MCMC control parameters changes in control.ergm
        \itemize{ 
          \item default for \code{MCMC.interval} increased from 100 to 1000
          \item default for \code{MCMC.burnin} changed from 10000 to mcmc.interval*16
          \item default for \code{MCMC.samplesize} decreased from 10000 to 1024
        }
      }
      \item{ MCMLE control parameter changes in control.ergm
        \itemize{
          \item a new control argument \code{MCMLE.termination} has been added to specify the criterion used for terminating MCMLE estimation.
          \item a new \code{MCMLE.metric} argument \code{'logtaylor'} has been added
          \item control parameter \code{MCMCLE.steplength.margin} has been added
          \item default value for \code{MCMLE.MCMC.precision} has been increased from 0.05 to 0.005
          \item default value for \code{MCMLE.steplength} has been changed from 0.5 to 1
          
        }
      }
      \item { ergm.control parameter \code{SA.nsubphases} default changed from MCMLE.maxit to 4}
      \item { control.gof parameter \code{MCMC.burnin} default increased from 1000 to 10000}
      \item { similarly, the default for the \code{MCMC.burnin} control parameter for \code{control.simulate}, \code{control.simulate.formula} and \code{control.simulate.formula.ergm} increased from 1000 to 10000}
      }
  }
  
   \subsection{BUG FIXES}{
    \itemize{
    \item{Networks where some statistics are at the lowest value, leading to some terms being -Inf, was not giving the correct estimated coefficients. Now fixed.}
    \item{Bug fix in Penalized-likelihood MPLE estimation}
    \item{Bug fix in MPLE estimation when ergm edgecov term has an attribute as a matrix.}
    \item{Bug fix in ergm nodecov term when attributes are specified only for the second mode in a bipartite network}
    \item{Various compatibility improvements to avoid warnings from R 3.1.2 and 3.1.3}
    \item{control parameter \code{MPLE.max.dyad.types} is now passed to \code{logLik.ergm()} and \code{ergm.bridge.dindstart.llk()}.
fixes issue reported by Sebastian Daza}
    \item{Fixed a fencepost bug in \code{ergm.bridge.dindstart.llk()}'s construction of the offset coefficient vector in the dyad-independent submodel reported by Philip Leifeld}
    \item{Adjusted package initialization code to avoid throwing error reported by Philip Leifeld when ergm imported by another package instead of attached to workspace.}
    \item {Fixed a bug with \code{ergm.pl} not dropping the appropriate rows when there are Infs. When network stats are 0, the glm model should now be specified correctly. Reported by Kevin Lewis}
    \item{Many other fixes...}
    }
   }

   \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item{ The ergm control parameter \code{MCMLE.method} has been deprecated.  By default, ergm uses \code{trust}, and falls back to \code{optim} with Nelder-Mead method when \code{trust} fails.}
      \item The ergm control parameter \code{MCMLE.adaptive.epsilon} has been removed.
      \item The ergm control parameters \code{MCMC.burnin.retries}, \code{MCMC.burnin.check.last} and \code{MCMC.burnin.check.alpha} have been removed. 
      \item \code{MCMLE.adaptive.epsilon} parameter has been removed
      \item The \code{faux.high} dataset was removed as it was identical to \code{faux.mesa.high}
    }
   }

  
   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
      \item{The "iterations" output in an ergm object now returns the number of MCMLE iterations taken by the estimation.}
      \item{New help page for ergm MetropolisHastings Proposal methods.}
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

\section{Changes in version 3.1.3}{

  % \subsection{NEW FEATURES}{
  %   \itemize{
  %   } 
  % }
  
  \subsection{BUG FIXES}{
    \itemize{
      \item An array bounds issue was fixed in the \code{C}
      code. It has not affected any previous output.

      \item Some spurious convergence warnings are no longer printed by
      \code{\link{summary.ergm}}.

      \item Bugs in the valued Uniform and Discrete Uniform
      Metropolis-Hastings proposals for missing dyads have been fixed.

      \item Term \code{\link{edgecov}} invoked with the network
      attribute syntax now checks for the presence of the network
      attribute in question.

    }
  }

  \subsection{DEPRECATED & DEFUNCT}{
    \itemize{
      \item Calling \code{\link{ergm}}, with \code{estimate="MPLE"} for
      a curved exponential family model (e.g., \code{gwesp} with
      \code{fixed=FALSE}) is now an error. (Previously, it treated the
      curved term as fixed, which was misleading.)
    }
  }

  
  % \subsection{OTHER USER-VISIBLE CHANGES}{
  %   \itemize{
      
  %   }
  % }

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

      \item MPLE for some degree-based constraints has been reenabled.
      
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
      increased to 100,000.
      
    }
  }
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
      rudimentary implementations for discrete uniform and continuous
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
      been re-implemented or fixed.

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
