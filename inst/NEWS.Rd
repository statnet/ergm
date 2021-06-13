%  File inst/NEWS.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm-package]{ergm}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.

%% \section{Changes in version X.Y.Z}{

%%   \subsection{NEW FEATURES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{BUG FIXES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{OTHER USER-VISIBLE CHANGES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }

  %% \subsection{C-LEVEL FACILITIES}{
  %%   \itemize{
  %%   }
  %% }

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

%% }

\section{Changes in version 4.0.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        Terms, such as \code{\link{nodecov}}, that take matrix-valued functions of nodal attributes now generate coefficient names from the matrix column names more sensibly.
      }
      \item{
        The \code{levels2} argument in terms \code{\link{mm}} and \code{\link{nodemix}} now accepts an argument in the form of a logical mixing matrix, selecting which mixing matrix cells should be included.
      }
      \item{
        \code{\link{ergmMPLE}()} now has an \code{expand.bipartite=} argument, which, for \code{output='array'}, controls whether bipartite networks' adjacency matrices are represented as rectangular or as square. Other behaviors under \code{output='array'} have also been made more consistent.
      }
      \item{
        In \code{\link{gof}} functions, bipartite networks now have their own goodness-of-fit statistics, \code{\link{b1degree}} and \code{\link{b2degree}}.
      }
      \item{
        Formula or function-type nodal attribute specifications can now be more complex. In particular, \code{\link{magrittr-package}}'s \code{\link{\%>\%}} operator can now be used on one of the sides of an \code{\link{mm}} formula.
      }
      \item{
        \code{\link{nodal_attributes}} documentation and vignette now have more examples.
      }
      \item{
        Nonidentifiability detection code now checks for nonvarying statistics first and handles them appropriately. See \code{\link{control.ergm}} for the parameters that control this handling.
      }
      \item{
        Control lists can now be specified as plain lists with sublists allocated by argument name. Function \code{\link{snctrl}} can provide completion in base \R, RStudio, and elsewhere.
      }
      \item{
        ERGM now checks for whether an MPLE exists and warns the user if not. (Thanks to Christian Schmidt (PSU).)
      }
      \item{
        Arguments to \code{\link{ergm}} terms can now be specified a little bit more flexibly; e.g., one can specify an unnamed argument after a named argument, as one can in an R function call.
      }

      \item{
	Curved ERGMs are no longer \dQuote{second class citizens}. Curved parameters can now be passed via
	\code{coef=} and \code{init=} arguments rather than embedded in the formula.
      }

      \item{
        A new type of term, an \emph{operator} is introduced, taking one or more \code{\link{ergm}} formulas as arguments and transforming their inputs and outputs allowing a wide variety of models to be specified. These include interactions for dyad-independent effects (\code{TERM1:TERM2} and \code{TERM1*TERM2}), \code{B()}inary terms in valued ERGMs, \code{Label()}, \code{F()}ilter, \code{Offset()}, \code{Symmetrize()}, \code{Sum()}, \code{S()}ubgraph, \code{Curve()}, \code{NodematchFilter()}, \code{Log()}, \code{Exp()}, \code{Prod()}, and others.
      }

      \item{
        A new type of MPLE, \code{"logitreg"}, invoked when the model is curved.
      }

      \item{
        A new term, \code{\link{attrcov}}, can cumulate statistics from multiple cells in the mixing matrix.
      }

      \item{
        Termes \code{\link{transitiveties}} and \code{\link{cyclicalties}} now work for undirected networks (where they count undirected transitive ties).
      }

      \item{
        \code{\link{is.valued}} is a family of methods for testing if a particular aspect of an ERGM is valued.
      }

      \item{
        Long MCMC, SAN, etc. runs can now be interrupted (i.e., with CTRL-C) even from inside the C code.
      }

      \item{
        \code{\link{update.network}} can now update the network using a variety of inputs, including edgelists and \code{\link{ergm_state}}s.
      }

      \item{
        The MPLE initial method can now be used to estimate curved ERGMs without specifying the initial curved parameter. The \code{initialfit=} \code{InitErgmTerm} argument and \code{as.initialfit=} argument of \code{\link{ergmMPLE}} have been removed. Similarly, \code{\link{fix.curved}} and \code{\link{enformulate.curved}} have been deprecated and may be removed in a future release.
      }

      \item{
	A mechanism for specifying \dQuote{observational constraints} (\code{obs.constraints})
	has been added to \code{\link{ergm}}, allowing more complex observation processes than
	just missing data to be estimated.

	Currently, this can be used for estimating arbitrary ERGMs on networks whose edges have been randomly flipped with known probabilities (say, to preserve privacy) (Karwa, Krivitsky, and Slavcovic, 2017) using the \code{\link{dyadnoise}} or observed egocentrically with alters identified (Handcock and Gile, 2010) using the \code{\link{egocentric}} constraint.

        An argument \code{observational=} to \code{\link{simulate.ergm}} can specify whether the model or the observational constraint from the \code{\link{ergm}} fit should be used.
      }

      \item{
        The user can now specify \dQuote{hints} on a formula to \code{MCMC.prop=} and \code{obs.MCMC.prop=} arguments passed to \code{\link{control.ergm}} that can allow the MCMC sampler to select an optimal proposal for the process.

        Currently, the available hints are \code{\link{sparse}} structure, implying that the network is sparse (the default), and \code{\link{strat}}, implying that the dyads are grouped in blocks which should be subject to stratified sampling.

        This is documented in the \emph{Proposal Lookup API} vignette.
      }

      \item{
        \code{\link{ergmMPLE}} output types now include \code{"dyadlist"}, which returns the same output as \code{"matrix"}, but without compressing the duplicate rows.
      }

      \item{
        New constraint, \code{\link{blocks}} has been implemented, allowing dyads for some combinations of a factor to be prohibited from toggling (i.e., fixed, or set as structural 0s).
      }

      \item{
        An MCMC proposal \code{\link{BDStratTNT}} has been implemented, allowing for highly efficient MCMC sampling for some special cases of bounded degree and block-level stratification among the dyads.
      }

      \item{
        A number of model specification settings, including \code{constraints}, \code{obs.constraints}, and \code{response} can be attached to the network itself using the \code{\link{\%ergmlhs\%}} operator and are inherited by the respective arguments unless overridden.
      }

      \item{
        Terms \code{\link{nodecovar}}, \code{\link{nodeocovar}}, and similar have been improved. (Terms \code{\link{nodesqrtcovar}} terms had been folded into the former.)
      }

      \item{
        \code{\link{san}} can now take into account offsets (finite and infinite).
      }

      \item{
        A new default MCMLE stopping criterion, \code{"confidence"} has been implemented. It performs a multivarate equivalence test, stopping if it can reject the null hypothesis that the the parameter guess is outside of a certain tolerance region around the true parameter value. If the test fails to reject, it increases the target effective size and resumes the optimization.
      }

      \item{
        Data structures for edgelists returned from MCMC, SAN, and Godfather are now allocated dynamically. This means that control parameters \code{*.init.maxedges} and \code{*.max.maxeges} are no longer useful and has been removed. \code{MCMC.maxedges} (for \code{\link{control.ergm}} and \code{\link{control.simulate.ergm}}) is retained (with an infinite default) to be used to short-circuit a sampler that has gone degenerate.
      }

      \item{
        The \code{"Stochastic-Approxmation"} algorithm (a C-based implementation of Robbins-Monro) now fully supports curved ERGMs.
      }

      \item{
        The \code{\link{control.ergm}} parameters for \code{"Stochastic-Approximation"} and \code{"Robbins-Monro"} have been reorganised.
      }

      \item{
        Limited support for multithreaded (i.e., multicore) evaluation of change statistics has been implemented. It has to be enabled at compile time and is controlled by \code{\link{set.MT_terms}} and related functions, and \code{parallel.inherit.MT=} argument of \code{\link{control.ergm}} and others.
      }

      \item{
        Controls for bridge sampling (i.e., \code{\link{control.ergm}}.bridge) can now be set separately for the constrained (missing data) sample.
      }


      \item{
        The standard argument \code{response=} can now be a \code{rlang}-style formula containing an expression in terms of edge attributes for the modelled edge value. If the expression (or the attribute) evalutes to a logical (TRUE/FALSE) variable, a binary ERGM is estimated.
      }

      \item{
        Adaptive (effective size) sampling is now possible for the \code{\link{simulate.ergm}} family of functions as well.
      }

      \item{
        \code{\link{control.simulate.ergm}} control parameter \code{MCMC.scale} can now be used to set the MCMC burn-in and interval relative to those of the \code{\link{ergm}} fit.
      }

      \item{
        In \code{\link{ergm}} formulas, the \code{offset()} decorator now takes a second optional argument allowing only some of the parameters in the term to be fixed.
      }

      \item{
        \code{\link{ergm}} is now a bit more careful about saving formulas and environments, saving memory and disk space.
      }

      \item{
        An experimental lack-of-fit diagnostic \code{\link{rank_test.ergm}} has been added.
      }

      \item{
        The \pkg{\link[sna:sna-package]{sna}} function \code{\link{symmetrize}} is made into a generic and a method for \code{\link{network}} objects is implemented.
      }

      \item{
	The user term and proposal APIs has been extended significantly. A new API enables storage of information between calls that is updated with the network state, \emph{auxiliaries}: terms that make available information to other terms, and many other capabilities. See the \emph{Terms API} vignette for the current specification. The API remains source-compatible with the previous versions.

        The \R API for \code{\link{ergm_model}} takes an additional argument \code{extra.aux=}, which can be used to request auxiliaries by other functions, such as proposals.
      }

      \item{
        A number of auxiliaries has been implemented (see \verb{R/InitErgmTerm.auxnet.R}) and a convenient API for evaluating network statistics on auxiliary networks is provided.
      }

      \item{
        Log-likelihood estimates and derived quantities (deviance, AIC, BIC) now also report their standard errors due to bridge sampling. Setting the \code{bridge.target.se} control option to \code{\link{control.logLik.ergm}} and/or \code{\link{control.ergm.bridge}} will continue sampling and accumulating samples until the estimated standard error due to bridge sampling falls below it.
      }

      \item{
        \code{\link{san}} output now includes the sampled statistics of the network.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{

      \item{
        The package functions now work more reliably when loaded but not attached but rather invoked via \code{::}.
      }
      \item{
        The \code{\link{Dyads}} operator constraint now works correctly on undirected networks.
      }
      \item{
        Term initialization traceback printing has been robustified.
      }
      \item{
        Term \code{\link{nodeisqrtcovar}} now initializes without error. Thanks to Yu-Lun Chen (NYU) for reporting the issue.
      }
      \item{
        The MPLE code has been robustified for some corner cases (e.g., a model with 0 free parameters).
      }
      \item{
        Various fixes to \code{\link{predict.ergm}}. (Thanks to Michal Bojanowski.)
      }
      \item{
        In \code{\link{ddsp}}, \code{\link{dnsp}}, \code{\link{dgwdsp}}, and \code{\link{dgwnsp}}, the \code{"RTP"} (reciprocated transitive partners) implementation has been fixed.
      }
      \item{
        Term initialization error messages are now somewhat cleaner and more informative.
      }
      \item{
        Miscellaneous documentation fixes and CRAN warning fixes.
      }

    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{

      \item{
        \strong{IMPORTANT:} The \code{levels2=} argument of \code{\link{nodemix}} and \code{\link{mm}} terms now defaults to \code{-1} (i.e., exclude the first level combination). Old behaviour is available through the \code{version} \link{term option}.
      }

      \item{
        \strong{IMPORTANT:} In curved terms such as \code{\link{gwesp}} (with \code{fixed=FALSE}), the decay parameter in the term is now ignored; the initial decay parameter must be specified as a part of \code{\link{control.ergm}(init=...)}. Note that because MPLE for curved ERGMs has been implemented, it is optional in the first place.
      }

      \item{
        \code{\link{summary_formula}} family of functions now dispatch based on the \code{basis=} argument if present, with the formula being allowed to be one-sided.
      }

      \item{
        In \code{\link{control.ergm}}, arguments \code{SAN.control=} and \code{loglik.control=} have been shortened to \code{SAN=} and \code{loglik=}, respectively.
      }

      \item{
        Global ERGM options and ERGM term options are now properly documented in \code{\link[=ergm-options]{options?ergm}}.
      }

      \item{
        The \code{hammingmix} term has been deprecaed.
      }

      \item{
        The convex hull testing now uses the \code{\link{lpSolveAPI-package}} and is faster and more robust as a result. (Thanks to Christian Schmid (PSU).)
      }
      
      \item{
	\code{\link{enformulate.curved}} family of functions has been
	deprecated.
      }

      \item{
        \code{ergmfit$coef} element has been superseded by \code{ergmfit$coefficients} for consistency with \code{\link{glm}} and others. The former still works, with a deprecation warning. Using \code{coef(ergmfit)} is recommended.
      }

      \item{
        Arguments to \code{\link{ergm_MCMC_sample}} and \code{\link{ergm_MCMC_slave}} have changed, and similarly for SAN and CD.
      }

      \item{
        Shared partner counts in \code{gw*sp} and \code{*sp} terms are now cached, resulting in a significant speedup. This can be turned off with the \link[=ergm-options]{term option} \code{cache.sp}.
      }

      \item{
	The term API now supports box constraints on parameters or
	coefficients, which may be set via \code{minpar} and
	\code{maxpar} elements in the \code{InitErgmTerm}'s return list.
      }

      \item{
        In the valued ERGM API, the \code{response} setting is now carried as an \code{\%ergmlhs\%} setting \code{"response"}.
      }

      \item{
        Some likelihood metrics, particularly the Median one, should work better now.
      }

      \item{
        In the \code{InitErgmConstraint} API, \code{free_dyads} can now be a function with 0 arguments, that is only evaluated if needed and is memoized (using the \pkg{\link[=memoise:memoise-package]{memoise}} package).
      }

      \item{
        \code{\link{control.ergm}} parameter \code{MPLE.constraints.ignore} can be used to suppress calculations of complex RLEBDMs for large networks if the constraint is imposed in some other fashion.
      }

      \item{
        \code{\link{control.ergm}} \code{obs.MCMC.*} parameters can now be set relative to the unconstrained MCMC parameters via \code{obs.MCMC.mul}.
      }

      \item{
        A variety of improvements and tuning parameters for the Hummel step length algorithm (\code{\link{control.ergm}} parameters \code{MCMLE.steplength.*}).
      }

      \item{
        \code{\link{control.ergm}}'s \code{*.samplesize} can now be specified as a multiple of the number of free parameters in the model, via the \code{*.samplesize.per_theta} parameter.
      }

      \item{
        Adaptive MCMC has been significantly improved. Methods of Vats, Flegal, and Jones (2019) are used to determine the multivariate effective sample size, and automatic burn-in detection now works by fitting a broken stick linear model to see where the sample statistics level out, then uses the Geweke convergence diagnostic to verify. See help for \code{\link{control.ergm}} and other functons for the available tuning parameters.
      }

      \item{
        A number of defaults for control parameters have changed.
      }

      \item{
        The \code{pending_update_network}, \code{ergm.Cprepare}, and other APIs for moving between R and C have been replaced by the \code{ergm_state} API that uses the \code{\link{.Call}} interface.
      }

      \item{
        \code{\link{simulate.ergm}} now has not just the \code{\link{ergm_model}} methods but also \code{\link{ergm_state}} methods.
      }

      \item{
        The \code{lasttoggle} API created for the \pkg{\link[tergm:tergm-packag]{tergm}} package has been removed in favour of the \dQuote{extended state} API that allows client packages to implement their own states. Similarly, \code{is.durational} family of methods has been removed.
      }

      \item{
        Information used by \code{\link{print.summary.ergm}} to explain the type of fit is now set in \code{\link{ergm}}.
      }

      \item{
        The \code{\link{ergm_proposal_table}} API now allows a more flexible specification of which constraints a proposal \emph{can} enforce as opposed to which ones it \emph{does} (i.e., can't not) enforce.
      }

      \item{
        The \code{InitErgmConstraint} and \code{InitErgmReference} API now has the same calling convention as \code{Init*ErgmTerm}, which enables \code{\link{check.ErgmTerm}} to be used to preprocess and check their arguments.
      }

      \item{
        The \code{Init*Ergm*} functions now receive an optional env= argument containing the environment of the formula from which they were called.
      }

      \item{
        The \code{InitErgmConstraint} API now allows multiple nested constraints to be returned, to allow a proposal to address either a more general or a more specific one.
      }

      \item{
        The \code{do.sim=} aregument to the \code{\link{simulate.ergm}} family of functions has been deprecated in favor of a more consistent \code{return.args=}.
      }

      \item{
        Likelihood-change trust region (\code{trustregion}) code has been removed.
      }
    }
  }


  \subsection{C-LEVEL FACILITIES}{
    \itemize{
      \item{
        Almost all of the \R-C interface has been transition to use the \code{\link{.Call}} rather than the \code{\link{.C}} API. Some helper functions are provided in \code{\link{inst/include/ergm_Rutils.h}}.
      }

      \item{
        Functions \code{\link{ergm.eta}}, \code{\link{ergm.etagrad}}, and \code{\link{ergm.etagradmult}} have been reimplemented in C using the \code{\link{.Call}} interface. They are available in the C API in \verb{inst/include/ergm_etamap.h}.
      }

      \item{
        \verb{inst/include/ergm_changestat.h} and \verb{inst/include/ergm_wtchangestat.h} now provide more helpful helper macros.
      }

      \item{
        \code{InitErgmTerm} and \code{InitErgmProposal} functions can now pass a \code{\link{integer}} (as opposed to a \code{\link{numeric}}) vector to the C functions with an \code{outlist$iinputs} element.
      }

      \item{
	A general API for setting callback functions that are called whenever a network changes is provided. See the \emph{Network Callback API} vignette.
      }


      \item{
        The discordance network API has been removed in favour of the storage API.
      }

      \item{
        The \verb{khash} macro library (from \url{https://github.com/attractivechaos/klib/}) has been incorporated into the package with some modifications and is available to use, in particular, for \eqn{O(1)} storage of dyad information, particularly in applications where efficient interation over edges incident on a node is not needed. It can be found in \verb{/inst/include/ergm_khash.h}.

        A variety of helper data structures and functions make use of it and are also available, particularly \verb{/inst/include/ergm_dyad_hashmap.h}.
      }

      \item{
        Additional helper APIs and data structures have been provided in \verb{/inst/include/}.
      }

      \item{
        A \verb{DyadGen} API implemented on the C side by \verb{/inst/include/ergm_dyadgen.h} and on the R side by \code{\link{ergm_dyadgen_select}} automatically selects sampling arbitrary dyads, sampling from an edgelist, or sampling from an \link[=rlebdm]{RLEBDM}. It obviates a number of other proposals, such as \code{RLETNT} and \code{listTNT}.
      }

      \item{
        All header files exported by \pkg{ergm} have been prefixed with \code{ergm_}. (Old header file names now issue deprecation warnings.)
      }
    }
  }
}

\section{Changes in version 3.11.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        \code{\link{ergm}} now checks models for evidence of nonidentifiability (linear dependence between statistics) after MPLE and during the MCMLE and CD estimation. QR decomposition is used to attempt to identify which statistics are redundant.
      }
      \item{
        It is now possible to limit the systematic sample size taken by \code{\link[=ergm.mple]{MPLE}} before resorting to the case-control likelihood using the \code{MPLE.samplesize} control parameter (if initial fit is final) and by the \code{init.MPLE.samplesize} control parameter (if it is used only to initialize the MCMC-based estimation).
      }
      \item{
        A new \dQuote{operator} constraint \code{\link[=ergm-constraints]{Dyads}} takes \code{ergm} formulas and constrains dyads based on which dyads affect their change statistics.
      }
      \item{
        It is now possible to calculate model-predicted \strong{conditional} and \strong{unconditional} tie probabilities using  \code{\link{predict.ergm}} and \code{\link{predict.formula}}. (Contributed by Michal Bojanowski.)
      }
      \item{
        The nodal attribute specification mechamism (i.e., the \code{attr=} and \code{levels=} arguments to terms) now includes a number of convenience functions and predicates, such as \code{\link[=nodal_attributes]{LARGEST}} and \code{\link[=nodal_attributes]{SMALLEST}} to select or exclude the most/least frequent levels, and \code{\link[=nodal_attributes]{COLLAPSE_SMALLEST}} to collapse the smallest categories into one.
      }
      \item{
        \code{\link{ergmMPLE}} now takes a \code{constraints} argument and handles it and missing dyads correctly.
      }
      \item{
        \code{\link{san}} now has the ability to take offsets into account.
      }
      \item{
        \code{\link{mcmc.diagnostics.ergm}} now takes an argument \code{which=} to  to specify which output is wanted.
      }
      \item{
        The following new model terms have been implemented:
        \itemize{
          \item{
            \code{\link{b1dsp}}, \code{\link{b2dsp}}, \code{\link{gwb1dsp}}, and \code{\link{gwb1dsp}}, the bipartite versions of \code{\link{dsp}} and \code{\link{gwdsp}}.
          }
          \item{
            \code{\link{isolatededges}}, to count the number of edges that are in a single component.
          }
        }
      }
      \item{
        A vignette for the \link[=nodal_attributes]{nodal attributes UI} has been added.
      }
      \item{The \code{\link{simulate.ergm}} family of functions' \code{output=} argument can no only takes a function through which the simulated networks are mapped before returning.
      }
      \item{
        The \code{\link[ergm]{cycle}} term implementation has been greatly improved.
      }
      \item{
        \code{\link{ergm}} now takes a \code{basis=} argument as others.
      }
      \item{
        \code{\link{plot.gof}} now plots a blue rhombus to indicate the simulated mean.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
        A bug involving assumptions based on the first entry in the \code{lasttoggle} vector has been corrected.
      }
      \item{
        In \code{\link{logLik.ergm}}, degrees of freedom no longer count offsets.
      }
      \item{
        The \code{Dyad} type used in C code is now at least 64 bits on all systems, fixing some issues around large networks.
      }
      \item{
        A bug that made \code{\link{san}} fail to terminate when it had hit its target stats exactly has been corrected.
      }
      \item{
        The \code{\link{as.rlebdm.edgelist}} function now correctly handles the network size attribute attached to the edgelist.
      }
      \item{
        Fixed incorrect use of \code{MCMC.*} control parameters in \code{\link{san}}.
      }
      \item{
        Fixed the documentation for the \code{object} argument to \code{\link{ergm_proposal}}, and fixed a recursive default argument to \code{\link{ergm_proposal.character}}.
      }
      \item{
        Added a workaround for a limitation of \code{coda::densityplot.mcmc.list} when variable names are not unique, which caused a problem in \code{\link{mcmc.diagnostics}} pointed out by James Wells.
      }
      \item{
        The \code{levels} argument was removed from \code{starmix} terms, as it was not functioning correctly.
      }
      \item{
        A bug in \code{\link{ergmMPLE}} for bipartite networks has been fixed.
      }
      \item{
        Some fixes to \code{\link{approx.hotelling.diff.test}}.
      }
      \item{
        Contrastive divergence runtime diagnostics now work correctly.
      }
      \item{
        \code{\link{ergm.godfather}} now handles undirected valued networks correctly.
      }
      \item{
        Per CRAN warning, none of the internal \code{\link{glm}}-alike functions default to \code{data=sys.parent()} anymore.
      }
      \item{
        \code{\link{as.network.numeric}} now handles the \code{bip} argument consistently with \code{\link{network.initialize}}.
      }
      \item{
        For bipartite terms, \code{NA}s in nodal attributes are now tolerated if they are in the unused bipartition.
      }
      \item{
        \code{\link{ergm.bridge}} now correctly updates the state of the network after every bridge's burn-in. (This bug did not invalidate any prior results, but it made the log-likelihoods and derived quantities reported by \code{\link{summary.ergm}} less accurate.)
      }
      \item{
        \code{\link{ergmMPLE}} now handles offsets corectly when returning matrix and array output.
      }
      \item{
        Miscellaneous fixes to \code{\link{predict.ergm}} and \code{\link{ergmMPLE}}, particularly in handling offsets.
      }
      \item{
        Miscellaneous vignette and documentation fixes.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        Datasets \code{fauxhigh} and \code{fauxhigh.withnumericattributes} have been removed.
      }
      \item{
        Control parameters \code{*.Hummel.*} have been renamed to \code{*.steplength.*} for consistency.
      }
      \item{
        Parallel multisection search for the Hummel step length is now only enabled for missing data MLE by default, and can be disabled altogether with appropriate control parameters.
      }
      \item{
        Formulas are no longer attached to \code{\link{ergm_model}}, \code{\link{ergm_proposal}}, \code{\link{ergm_conlist}}, and others.
      }
      \item{
        The \code{\link[=simulate.ergm]{simulate}} family of methods has been restructured. In particular, it is now possible to call \code{simulate()} on an \code{\link{ergm_model}} object without reinitializing it.
      }
      \item{
        A number of improvements have been made to the nodal attribute API. In particular, \code{\link{ergm_get_vattr}} now has an \code{AsIs} method that makes minimal changes to the input before returning, \code{\link{ergm_get_vattr.function}} now propagates its arguments to the function it calls, and \code{\link{ergm_get_vattr}} now accepts a \code{"index"} type argument representing vertex indices (a logical of appropriate length, or a nonnegative or nonpositive integer vector).
      }
      \item{
        Missing nodal attributes now generate an error in model initialization.
      }
      \item{
        \code{\link{summary.ergm}} and its printing method have been cleaned up.
      }
      \item{
        A variety of miscellaneous optimizations, coleanups, and documentation improvements.
      }
      \item{
        The \code{InitErgmTerm} API's \code{\link{check.ErgmTerm}} now accepts \code{NA} or \code{""} for argument type, indicating that any type is accepted.
      }
      \item{
        In parallel computing, worker nodes now load packges \code{\link{.libPaths}} of the manager node.
      }
      \item{
        The following deprecated functions have been made defunct: \code{ergm.update.formula}, \code{remove.offset.formula}, \code{network.update}, \code{ergm.getmodel}, \code{ergm.getglobalstats}, \code{as.edgelist.compressed}, \code{as.network.uncompressed}, \code{standardize.network}, and \code{newnw.extract}.
      }
      \item{
        \code{san.ergm} has been deprecated, since it lacks a known use case.
      }
      \item{
        The \code{\link{simulate.formula}} dispatching mechanism has been rewritten to coexist with other pacakges that want to dispatch based on the LHS of the formula, particularly \code{lme4}. The \code{\link{simulate_formula}} API still works, but it may be phased out in the future.
      }
      \item{
        \code{\link{ergm.bridge}} has has been optimised to avoid repeated calls to model and proposal initialisation.
      }
    }
  }

  \subsection{C-LEVEL FACILITIES}{
    \itemize{
      \item{C export \code{ergm_constants.h} now exports the major and the minor version of \code{\link[ergm-package]{ergm}}.
      }
    }
  }
}

\section{Changes in version 3.10.4}{

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
	Deprecated term argument warnings are now only printed once per session.
      }
      \item{
	\code{mcparallel} antiwarnings are now printed less frequently.
      }
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item{
	\code{ergm} version 3.10.1 inadvertently depended on \R 3.6. The dependence has been reverted to 3.5.
      }
    }
  }
}

\section{Changes in version 3.10.1}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
	In all \code{gw*} terms, the \code{cutoff=} argument can now be controlled by a \code{gw.cutoff} term option.
      }

      \item{
	\code{\link[=nobs.ergm]{nobs}} method has now been implemented for \code{\link{ergm}} objects.
      }

      \item{
	\code{\link{simulate.ergm}} now takes a \code{basis=} argument (defaulting to NULL) that it passes to simulate.formula().
      }

      \item{
	\code{\link{ergm.godfather}} now takes an additional \code{changes.only=} argument, to return statistics relative to initial network.
      }

      \item{
	\code{\link{san}}'s internals have been singificantly updated. In particular,
	\itemize{

	  \item{SAN now performs true simulated annealing, and arguments are interpreted somewhat differently.}

	  \item{SAN's weighting of the statistics is now determined adaptively, in approximate proportion to the standard deviation of proposal distribution of these statistics.}

	}
      }

      \item{
	\code{\link{fix.curved}} and \code{\link{enformulate.curved}} now handle the \code{\link[=dgwesp]{dgw*sp}} terms.
      }

      \item{
	\code{\link{simulate.ergm}} and others now have an \code{output=} argument to provide finer control over the return format. It is to be preferred over the \code{statsonly=} argument.
      }

      \item{
	A number of terms (see terms help) now use the new \link[nodal_attributes]{nodal attributes UI}.
      }

      \item{
	Hummel steplength search can now be run in parallel.
      }

      \item{
	Terms \code{node*cov} and \code{b*cov} can now take multiple quantitative attributes or a function that returns a quantitative attribute matrix and produce a statistic for each.
      }

      \item{
	Terms \code{\link{b1sociality}} and \code{\link{b2sociality}} have been added.
      }

      \item{
	MCMLE estimation method can now checkpoint and resume optimization. See \code{\link{control.ergm}}'s \code{MCMLE.checkpoint=} and \code{MCMLE.resume=} arguments.
      }

      \item{
	A variety of documentation expansions, cleanup, and other improvements.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
	The default method of \code{\link{ergm.stopCluster}} had not
	been exported.
      }
      \item{
	\code{\link{summary.ergm}} fixed for CD estimation.
      }
      \item{
	\code{\link{diff}} term with \code{sign.action="abs"} no longer produces a warning on bipartite undirected networks. Thanks to Damien Georges for pointing out that this warning is unnecessary.
      }
      \item{
	More robust searching for initialization functions.
      }
      \item{
	Some spurious warnings fixed.
      }

      \item{
	\code{ergm} version is now stored in the fitted object; this should produce sensible warnings when, say, trying to compute a \code{\link{summary.ergm}} on a object produced by an earlier version.
      }

      \item{
	Disambiguation help pages have been added for concepts (such as \code{edges} and \code{hamming}) that have multiple roles (such as as an ERGM term and as an ERGM sample space constraint), removing duplication of aliases.
      }

      \item{
	Miscellaneous term and example documentation fixes. (Thanks to Christian Loehle and others.)
      }

    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
	Low-level function \code{\link{ergm_SAN_slave}} has now been split out of \code{\link{san}} like \code{\link{ergm_MCMC_slave}}.
      }

      \item{
	Lower-level \code{\link{simulate.ergm_model}} and \code{\link{san.ergm_model}} have been split out of \code{\link{simulate.formula}} and \code{\link{san.formula}}.
      }

      \item{
	\code{\link{gof.ergm}} now uses a longer \code{MCMC.interval} by default.
      }

      \item{
	Parallel processing \dQuote{API} now handles cluster starting and stopping more intelligently. In particular, \code{ergm.getCluster} now automatically sets up the calling function to stop the cluster when the calling function terminantes, as well as keeping track of which packages have already been loaded on the slave nodes.
      }

      \item{
	The vignette no longer tests egocentric code.
      }

      \item{
	Failure to allocate memory in the C routines now produces sensible errors.
      }

      \item{
	A semi-internal \code{\link{pending_udpate_network}} API can be used to carry around network information between \code{.C} calls without having to reconstruct the network.
      }

      \item{
	\code{\link{control.ergm}}'s \code{init.method=} argument now overrides the candidate list when selecting initialization methods.
      }

      \item{
	\code{\link{network.update}} is deprecated in favor of an \code{\link{update.network}} method.
      }

      \item{
	\code{\link{standardize.network}} has been deprecated, rendered obsolite by improvements to the \code{network} package method.
      }

      \item{
	\code{\link{simulate.ergm}} now defaults to starting the simulation from the \code{$newnetwork} element of the fit object. Thanks to Saeid Ghafouri for raising the issue and Sam Jenness, Steven Goodreau, and Martina Morris for suggesting a solution.
      }

      \item{
	\code{simulate} methods now attach an \code{attr(,"monitored")} to the returned statistics, indicating which statistics were generated from the \code{monitor=} argument.
      }

      \item{
	\code{\link{pending_update_network}} \dQuote{API} is more formalized now, including methods overriding those of \code{\link{network}}.
      }

      \item{
	\code{control$init.method} now overrides the candidate list when selecting initialization methods.
      }

      \item{
	The main vignette no longer contains egocentric inference code. Refer to the \code{ergm.ego} package.
      }
      \item{
	Improvements to \code{\link{plot.gof}} for large networks.
      }

      \item{
	\code{logLik.ergm(warn.dyads=)} control parameter has been replaced by \code{options(ergm.loglik.warn_dyads=)} option.
      }

      \item{
	A number of previously deprecated functions are now defunct; a
	number of previously defunct functions have been removed.
      }

      \item{SAN's documentation has been augmented with examples.}

      \item{
	\code{\link{sociality}} term's arguments \code{attrname} and \code{levels} have been deprecated.
      }

      \item{
	Most parallel calls are now a little bit more fault-tolerant, retrying on error.
      }

      \item{
	Chad Klumb has been added as a contributor to the package.
      }

    }
  }

  \subsection{C-LEVEL FACILITIES}{
    \itemize{
      \item{Some changes to the Metropolis-Hastings proposal API.}

      \item{Network initialization functions now return a pointer, and destroying functions free it.}

      \item{
	Miscellaneous improvements to low-level C code. In particular, more of the code now runs deterministically.
      }

      \item{
	Improvements to encapsulation of \code{Network} structures.
      }
    }
  }

  %% \subsection{UTILITIES}{
  %%   \itemize{
  %%     \item Likewise, "UTILITIES" is probably unnecessary.
  %%   }
  %% }

  %% \subsection{INSTALLATION}{
  %%   \itemize{
  %%     \item ... as is "INSTALLATION".
  %%   }
  %% }

  %% \subsection{PACKAGE INSTALLATION}{
  %%   \itemize{
  %%     \item ditto.
  %%   }
  %% }

}

\section{Changes in version 3.9.4}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
	A new API allows arbitrary combinations of dyad-level sample
	space constraints.
      }
      \item{
	Most terms that take categorical nodal attributes now have a \code{levels=} argument, which can be used to control the set and the ordering of the attribute levels.
      }
      \item{
	More informative error messages for term initialization.
      }
      \item{
	\code{\link{ergm.godfather}} has been reimplemented in the
	\code{ergm} package, along with its valued counterpart.
      }
      \item{
	New terms have been added or improved:
	\itemize{
	  \item{
	    \code{mm} for a more flexible specification
	    for mixing matrix and marginal effects;
	  }
	  \item{
	    valued versions of \code{sender}, \code{receiver},
	    \code{sociality} terms;
	  }
	  \item{
	    \code{ininterval} now allows interval openness to be
	    specified as a character string.
	  }
	}
      }
      \item{A new \code{term.options} control parameter allows additional
	arguments to be passed to term initializers.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
	A numerical instability in penalised pseudolikelihood
	calculation identified by Kevin Reuning has been fixed.
      }
      \item{
	An error in runtime traceplotting for models with offset terms
	reported by Kevin Reuning has been fixed.
      }
      \item{
	An error when trying to run the Robbins-Monro algorithm with parallel
	processing enabled, reported by George G. Vega Yon has been fixed.
      }
      \item{
	A bug in the \code{\link{meandeg}} term for directed
	networks has been fixed.
      }
      \item{
	Miscellaneous bugs in the \code{gof} family of functions have
	been fixed.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
	The coefficient names associated with terms enclosed in
	\code{offset()} are now also enclosed in \code{offset()}.
      }
      \item{
	\code{\link[=summary.ergm]{summary}} method for \code{ergm} now
	uses uses \eqn{z}-tests rather than \eqn{t}-tests, since there
	is no \dQuote{dispersion} parameter. In addition, it no longer
	takes a \code{digits=} argument.
      }
      \item{
	Handling of missing relations has been improved.
      }
      \item{
	Estimation progress when \code{verbose=TRUE} is now a bit less
	verbose; use \code{verbose=2} and up to get more information.
      }
      \item{
	The default value of the \code{eval.loglik=} argument to \code{\link{ergm}} is now
	controlled by \code{\link{options}} option \code{ergm.eval.loglik=}.
      }
      \item{
	\code{\link{plot.network.ergm}} has been deprecated.
      }
      \item{
	Some of the unofficial APIs have been changed for consistency
	and ease of use:
	\itemize{
	  \item{
	    A new set of helper functions has been added to
	    improve handling of nodal attributes. See
	    \link{nodal_attributes-API} for details.
	  }
	  \item{
	    \code{\link[=summary.ergm]{summary}} method for \code{ergm}
	    now exports a \code{$coefficients} element (like
	    \code{\link{summary.lm}}) and \code{$devtable} is a proper table.
	  }
	  \item{
	    New generics \code{\link{nparam}} and
	    \code{\link{param_names}} can be used to query the number
	    and the names of the parameters in a model.
	  }
	  \item{
	    \code{\link{ergm_model}} (replacing
	    \code{ergm.getmodel}) is now a full-fledged S3 class,
	    representing initialized \code{ergm} models. See its help for the
	    major methods. In particular, \code{ergm.getglobalstats} has
	    been repalced with \code{\link{summary.ergm_model}}.
	  }
	  \item{
	    \code{\link{simulate.formula}} functions now take an
	    \code{\link{ergm_model}} as an argument to \code{monitor=}.
	  }
	  \item{
	    Constraint initializers have been renamed from
	    \code{InitConstraint} to \code{InitErgmConstraint}. The
	    constraint implications table has been replaced with
	    \code{implies} and \code{impliedby} constraint lists.
	  }
	  \item{
	    Proposal initializers have been renamed from
	    \code{InitMHP} (and \code{InitWtMHP}) to
	    \code{InitErgmProposal} (and
	    \code{InitWtErgmProposal}). \code{ergm.MHP.table} has been renamed
	    to \code{\link{ergm_proposal_table}}.
	  }
	  \item{
	    Reference initializers have been renamed from
	    \code{InitReference} to \code{InitErgmReference}. Return
	    format has also been changed, replacing
	    \code{ergm.init.methods}.
	  }
	  \item{
	    \code{ergm.Clist.*} functions have been refined and
	    specialized. In particular, those that specifically return
	    numeric vectors (not lists of numeric vectors) for passing
	    to C code are now methods of the
	    \code{\link{to_ergm_Cdouble}} generic.
	  }
	  \item{
	    \code{ergm.getMCMCsample} and \code{ergm.mcmcslave} have
	    been renamed to \code{ergm_MCMC_sample} and
	    \code{ergm_MCMC_slave} respectively, and the arguments of
	    \code{ergm_MCMC_sample} have been made more user-friendly.
	  }
	}
      }
    }
  }

  %% \subsection{C-LEVEL FACILITIES}{
  %%   \itemize{
  %%   }
  %% }

  %% \subsection{UTILITIES}{
  %%   \itemize{
  %%     \item Likewise, "UTILITIES" is probably unnecessary.
  %%   }
  %% }

  %% \subsection{INSTALLATION}{
  %%   \itemize{
  %%     \item ... as is "INSTALLATION".
  %%   }
  %% }

  %% \subsection{PACKAGE INSTALLATION}{
  %%   \itemize{
  %%     \item ditto.
  %%   }
  %% }

}


\section{Changes in version 3.8.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
	Development of \code{ergm} and associated packages is now hosted
	on GitHub.
      }
      \item{
	Some terms, particularly \code{\link{dyadcov}} have been improved.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
	Bug fixes to contrastive divergence code.
      }

      \item{
	\code{\link{ddsp}} and \code{\link{dnsp}} terms now work
	correctly for 0 shared partners. (Thanks to Wouter Spekkink for
	reporting this bug.)
	}

      \item{
	Miscellaneous bug fixes, typo corrections, and robustifications.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
	Terms \code{degreepopularity}, \code{idegreepopularity}, and
	\code{odegreepopularity} have been renamed to \code{degree1.5},
	\code{idegree1.5}, and \code{odegree1.5} respectively.
      }
      \item{
	Progress reports for ERGM fits and simulations are now printed
	using \code{\link{message}}, so they can be suppressed with
	\code{\link{suppressMessages}}.
      }
      \item{
	A number of functions deprecated in 2015 have been removed. See
	\code{help('ergm-defunct')} for the list.
      }
    }
  }
}

\section{Changes in version 3.7.1}{


  \subsection{BUG FIXES}{
    \itemize{
      \item{
	Reverted a non-portable use of C function
	\code{alloca()}. (Thanks to Prof. Brian Ripley for pointing out
	the problem.)
      }
    }
  }


}


\section{Changes in version 3.7.0}{

  \subsection{NEW FEATURES}{
    \itemize{

      \item{
	New term \code{\link[ergm]{diff}}, a directed analogue to
	\code{\link{absdiff}} has been added, allowing for different
	ways to handle differences in nodal attributes.
      }

      \item{
	By default, \code{\link{gof}} now checks goodness-of-fit
	with respect to the estimated model in addition to the specified
	terms. This can be used to diagnose convergence.
      }


    }
  }


   \subsection{BUG FIXES}{
     \itemize{

       \item{Terms \code{\link{concurrent}}, \code{\link{degree}},
	 \code{\link{degrange}}, and \code{concurrentties}
	 \dQuote{\code{byarg}} argument has been corrected to
	 \dQuote{\code{by}}.}

       \item{The bug in the Contrastive Divergence code that caused ERGM
	 fits with \code{bd} constraint to hang has been fixed.}

       \item{Documentation for the curved ERGM terms has been clarified
	 in a number of places.}

       \item{Miscellaneous fixes and robustifications in the step-length
	 code.}

     }
   }


   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
       \item{
	 The geometric decay parameter in \code{gw*sp} and
	 \code{dgw*sp} terms has been renamed from \code{alpha} to
	 \code{decay}, for consistency with \code{gw*degree} terms.
       }

       \item{
	 It is now explicitly an error to pass a \code{control=}
	 parameter to the function directly. (Previously, such
	 parameters were silently ignored.)
       }
     }
   }
}


\section{Changes in version 3.6.1}{

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
	(\code{\link{ddsp}}, \code{\link{desp}}, and \code{\link{dnsp}})
	have been implemented. See \code{\link{ergm-terms}} for help.
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
      \item {MCMC control parameters changes in \code{\link{control.ergm}}
        \itemize{
          \item default for \code{MCMC.interval} increased from 100 to 1000
          \item default for \code{MCMC.burnin} changed from 10000 to mcmc.interval*16
          \item default for \code{MCMC.samplesize} decreased from 10000 to 1024
        }
      }
      \item{ MCMLE control parameter changes in \code{\link{control.ergm}}
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
      \code{\link[ergm.count:ergm.count-package]{ergm.count}})

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
      \code{\code{\link{control.ergm}}(drop=TRUE)}) can now detect both extremely high
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
