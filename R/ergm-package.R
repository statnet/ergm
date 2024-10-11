#  File R/ergm-package.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @details \insertNoCite{HuHa08e,KrHu23e}{ergm}
#'
#' For a complete list of the functions, use \code{library(help="ergm")} or
#' read the rest of the manual. For a simple demonstration, use
#' \code{demo(packages="ergm")}.
#' 
#' When publishing results obtained using this package, please cite the
#' original authors as described in \code{citation(package="ergm")}.
#' 
#' All programs derived from this package must cite it. Please see the
#' file `LICENSE` and [`https://statnet.org/attribution`](https://statnet.org/attribution).
#' 
#' Recent advances in the statistical modeling of random networks have had an
#' impact on the empirical study of social networks. Statistical exponential
#' family models (Strauss and Ikeda 1990) are a generalization of the Markov
#' random network models introduced by \insertCite{FrSt86m;textual}{ergm}, which in turn
#' derived from developments in spatial statistics \insertCite{Be74s}{ergm}. These models
#' recognize the complex dependencies within relational data structures.  To
#' date, the use of stochastic network models for networks has been limited by
#' three interrelated factors: the complexity of realistic models, the lack of
#' simulation tools for inference and validation, and a poor understanding of
#' the inferential properties of nontrivial models.
#' 
#' This manual introduces software tools for the representation, visualization,
#' and analysis of network data that address each of these previous
#' shortcomings.  The package relies on the [`network`]
#' package which allows networks to be represented in . The
#' \CRANpkg{ergm} package implements maximum likelihood
#' estimates of ERGMs to be calculated using Markov Chain Monte Carlo (via
#' [ergm()]). The package also provides tools for simulating networks
#' (via [simulate.ergm()]) and assessing model goodness-of-fit (see
#' [mcmc.diagnostics()] and [gof.ergm()]).
#'
#' A number of Statnet Project packages extend and enhance
#' \CRANpkg{ergm}. These include
#' \CRANpkg{tergm} (Temporal ERGM), which provides
#' extensions for modeling evolution of networks over time;
#' \CRANpkg{ergm.count}, which facilitates
#' exponential family modeling for networks whose dyadic measurements are
#' counts; and
#' \pkg{ergm.userterms}, available on GitHub at \url{https://github.com/statnet/ergm.userterms}, which
#' allows users to implement their own ERGM terms.
#' 
#' For detailed information on how to download and install the software, go to
#' the \CRANpkg{ergm} website: \url{https://statnet.org}. A
#' tutorial, support newsgroup, references and links to further resources are
#' provided there.
#'
#' @seealso [`ergmTerm`], [`ergmConstraint`], [`ergmReference`],
#'   [`ergmHint`], and [`ergmProposal`] for indices of model
#'   specification and estimation components visible to the \CRANpkg{ergm}'s API at any given time.
#' 
#' @references \insertAllCited{}
#'
#' Admiraal R, Handcock MS (2007).  \CRANpkg{networksis}: Simulate
#' bipartite graphs with fixed marginals through sequential importance
#' sampling.  Statnet Project, Seattle, WA.  Version 1,
#' \url{https://statnet.org}.
#' 
#' Bender-deMoll S, Morris M, Moody J (2008).  Prototype Packages for Managing
#' and Animating Longitudinal Network Data: \pkg{dynamicnetwork} and
#' \pkg{rSoNIA}.  \emph{Journal of Statistical Software}, 24(7).
#' \doi{10.18637/jss.v024.i07}
#' 
#' Boer P, Huisman M, Snijders T, Zeggelink E (2003).  StOCNET: an open
#' software system for the advanced statistical analysis of social networks.
#' Groningen: ProGAMMA / ICS, version 1.4 edition.
#' 
#' Butts CT (2007).  \CRANpkg{sna}: Tools for Social Network Analysis.  R package
#' version 2.3-2. \url{https://cran.r-project.org/package=sna}
#' 
#' Butts CT (2008).  \CRANpkg{network}: A Package for Managing Relational Data in .
#' \emph{Journal of Statistical Software}, 24(2).
#' \doi{10.18637/jss.v024.i02}
#' 
#' Butts C (2015). \CRANpkg{network}: Classes for Relational Data. The Statnet
#' Project (\url{https://statnet.org}). R package version 1.12.0,
#' \url{https://cran.r-project.org/package=network}.
#' 
#' Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \CRANpkg{statnet} Tutorial.  \emph{Journal of Statistical Software}, 24(8).
#' \doi{10.18637/jss.v024.i08}
#' 
#' Goodreau SM, Kitts J, Morris M (2008b).  Birds of a Feather, or Friend of a
#' Friend? Using Exponential Random Graph Models to Investigate Adolescent
#' Social Networks.  \emph{Demography}, 45, in press.
#' 
#' Handcock, M. S. (2003) Assessing Degeneracy in Statistical Models of Social
#' Networks, Working Paper #39, Center for Statistics and the Social Sciences,
#' University of Washington.
#' \url{https://csss.uw.edu/research/working-papers/assessing-degeneracy-statistical-models-social-networks}
#' 
#' Handcock MS (2003b).  \CRANpkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.0,
#' \url{https://statnet.org}.
#' 
#' Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
#' \CRANpkg{statnet}: Software Tools for the Statistical Modeling of Network Data.
#' Statnet Project, Seattle, WA.  Version 3, \url{https://statnet.org}.
#' 
#' Hunter, D. R. and Handcock, M. S. (2006) Inference in curved exponential
#' family models for networks, \emph{Journal of Computational and Graphical
#' Statistics}, 15: 565-583
#' 
#' Krivitsky PN, Handcock MS (2007).  \CRANpkg{latentnet}: Latent position and
#' cluster models for statistical networks.  Seattle, WA.  Version 2,
#' \url{https://statnet.org}.
#' 
#' Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. \emph{Electronic Journal of Statistics}, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#' 
#' Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' \emph{Journal of Statistical Software}, 24(4).
#' \doi{10.18637/jss.v024.i04}
#' 
#' Strauss, D., and Ikeda, M.(1990). Pseudolikelihood estimation for social
#' networks. \emph{Journal of the American Statistical Association}, 85,
#' 204-212.
#' @keywords models internal
"_PACKAGE"


#' Terms used in Exponential Family Random Graph Models
#'
#' @name ergmTerm
#' @aliases ergm-terms ergm.terms terms-ergm terms.ergm InitErgmTerm InitErgmWtTerm
#' @description This page explains how to specify the network statistics \eqn{g(y)} to functions in the [`ergm`][ergm-package] package and packages that extend it. It also provides an indexed list of the possible terms (and hence network statistics) visible to the \CRANpkg{ergm} API. Terms can also be searched via [`search.ergmTerms`], and help for an individual term can be obtained with `ergmTerm?<term>` or `help("<term>-ergmTerm")`.
#'
#' @section Specifying models:
#' \ERGMspec
#'
#' Network statistics \eqn{g(y)} and mappings \eqn{\eta(\theta)} are specified by a formula object, of the form `y ~ <term 1> + <term 2> ...`, where
#' `y` is a network object or a matrix that can be coerced to a network
#' object, and `<term 1>`, `<term 2>`, etc, are each terms chosen
#' from the list given below.  To create a network object in , use the
#' [`network`] function, then add nodal attributes to it
#' using the `%v%` operator if necessary.
#' 
#' ## Term operators
#' Operator terms like \ergmTerm{ergm}{B}{()} and \ergmTerm{ergm}{F}{()} take
#' formulas with other [`ergm`] terms as their arguments and transform them
#' by modifying their inputs (e.g., the network they evaluate) and/or their
#' outputs.
#' 
#' By convention, their names are capitalized and CamelCased.
#' 
#' ## Interactions
#' For binary ERGMs, interactions between [`ergm`] terms can be
#' specified in a manner similar to [`lm`] and others, as using the
#' `:` and `*` operators. However, they must be interpreted
#' carefully, especially for dyad-dependent terms. (Interactions involving
#' curved terms are not supported at this time.)
#' 
#' Generally, if term `a` has \eqn{p_a}{p[a]} statistics and `b` has
#' \eqn{p_b}{p[b]}, `a:b` will add \eqn{p_a \times p_b}{p[a]*p[b]}
#' statistics to the model, corresponding to each element of
#' \eqn{g_a(y)}{g[a](y)} interacted with each element of \eqn{g_b(y)}{g[b](y)}.
#' 
#' The interaction is defined as follows. Dyad-independent terms can be
#' expressed in the general form \eqn{g(y;x)=\sum_{i,j} }{sum[i,j]
#' x[i,j]*y[i,j]}\eqn{ x_{i,j}y_{i,j}}{sum[i,j] x[i,j]*y[i,j]} for some edge
#' covariate matrix \eqn{x}, \deqn{g_{a:b}(y)=\sum_{i,j}
#' x_{a,i,j}x_{b,i,j}y_{i,j}.}{g[a:b](y) = \sum[i,j] x[a,i,j]*x[b,i,j]*y[i,j].}
#' In other words, rather than being a product of their sufficient statistics
#' (\eqn{g_{a}(y)g_{b}(y)}{g[a](y)*g[b](y)}), it is a dyadwise product of their
#' dyad-level effects.
#' 
#' This means that an interaction between two dyad-independent terms can be
#' interpreted the same way as it would be in the corresponding logistic
#' regression for each potential edge. However, for undirected networks in
#' particular, this may lead to somewhat counterintuitive results. For example,
#' given two nodal covariates `"a"` and `"b"` (whose values for node
#' \eqn{i} are denoted \eqn{a_i}{a[i]} and \eqn{b_i}{b[i]}, respectively),
#' `nodecov("a")` adds one statistic of the form \eqn{\sum_{i,j}
#' (a_{i}+a_{j}) y_{i,j}}{sum[i,j] (a[i]+a[j])*y[i,j]} and analogously for
#' `nodecov("b")`, so `nodecov("a"):nodecov("b")` produces
#' \deqn{\sum_{i,j} (a_{i}+a_{j}) (b_{i}+b_{j}) y_{i,j}.}{sum[i,j]
#' (a[i]+a[j])*(b[i]+b[j])*y[i,j].}
#' 
#' ## Binary and valued ERGM terms
#' [`ergm`][ergm-package] functions such as [`ergm`] and
#' [`simulate`][simulate.formula] (for ERGMs) may operate in two
#' modes: binary and weighted/valued, with the latter activated by passing a
#' non-NULL value as the `response` argument, giving the edge attribute
#' name to be modeled/simulated.
#' 
#' ### Generalizations of binary terms
#' Binary ERGM statistics cannot be
#' used directly in valued mode and vice versa. However, a substantial number
#' of binary ERGM statistics --- particularly the ones with dyadic independence
#' --- have simple generalizations to valued ERGMs, and have been adapted in
#' [`ergm`][ergm-package]. They have the same form as their binary
#' ERGM counterparts, with an additional argument: `form`, which, at this
#' time, has two possible values: `"sum"` (the default) and
#' `"nonzero"`. The former creates a statistic of the form \eqn{\sum_{i,j}
#' x_{i,j} y_{i,j}}{sum[i,j] x[i,j]*y[i,j]}, where \eqn{y_{i,j}}{y[i,j]} is the
#' value of dyad \eqn{(i,j)} and \eqn{x_{i,j}}{x[i,j]} is the term's covariate
#' associated with it. The latter computes the binary version, with the edge
#' considered to be present if its value is not 0.  Valued version of some
#' binary ERGM terms have an argument `threshold`, which sets the value
#' above which a dyad is conidered to have a tie. (Value less than or equal to
#' `threshold` is considered a nontie.)
#' 
#' The \ergmTerm{ergm}{B}{()} operator term documented below can be used to pass other
#' binary terms to valued models, and is more flexible, at the cost of being
#' somewhat slower.
#' 
#' ## Nodal attribute levels and indices
#' Terms taking a categorical nodal covariate also take the `levels`
#' argument.  (There are analogous `b1levels` and `b2levels`
#' arguments for some terms that apply to bipartite networks, and the
#' `levels2` argument for mixing terms.)  The `levels` argument can
#' be used to control the set and the ordering of attribute levels.
#' 
#' Terms that allow the selection of nodes do so with the `nodes`
#' argument, which is interpreted in the same way as the `levels`
#' argument, where the categories are the relevant nodal indices themselves.
#' 
#' Both `levels` and `nodes` use the new level selection UI. (See
#' \link[=nodal_attributes]{Specifying Vertex attributes and Levels} (\verb{?
#' nodal_attributes}) for details.)
#' 
#' ### Legacy arguments
#' 
#' The legacy `base` and `keep` arguments are deprecated as of
#' version 3.10, and replaced by the `levels` UI. The `levels`
#' argument provides consistent and flexible mechanisms for specifying which
#' attribute levels to exclude (previously handled by `base`) and include
#' (previously handled by `keep`).  If `levels` or `nodes`
#' argument is given, then `base` and `keep` arguments are ignored.
#' The legacy arguments will most likely be removed in a future version.
#' 
#' Note that this exact behavior is new in version 3.10, and it differs
#' slightly from older versions: previously if both `levels` and
#' `base`/`keep` were given, `levels` argument was applied first
#' and then applied the `base`/`keep` argument. Since version 3.10,
#' `base`/`keep` would be ignored, even if old term behavior is
#' invoked (as described in the next section).
#' 
#' ## Term versioning
#' When a term's behavior has changed from prior version, it is often possible
#' to invoke the old behavior by setting and/or passing a `version` term
#' option, giving the verison (constructed by [`as.package_version`])
#' desired.
#' 
#' ## Custom `ergm` terms
#' Users and other packages may build custom terms, and package
#' \pkg{ergm.userterms} (\url{https://github.com/statnet/ergm.userterms}) provides
#' tools for implementing them.
#' 
#' The current recommendation for any package implementing additional terms is
#' to document the term with Roxygen comments and a name in the form
#' `termName-ergmTerm`. This ensures that \code{help("ergmTerm")} will list ERGM
#' terms available from all loaded packages.
#'
#' @section Terms included in the [`ergm`][ergm-package] package:
#' As noted above, a cross-referenced HTML version of the term documentation is
#' also available via `vignette('ergm-term-crossRef')` and terms
#' can also be searched via [`search.ergmTerms`].
#'
#' \ergmCSS
#'
#' ## Term index (plain)
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~!"operator"%in%.))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~!"operator"%in%.))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~!"operator"%in%.))}}
#'
#' ## Term index (operator)
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~"operator"%in%.))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~"operator"%in%.))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmTerm", keywords = ~"operator"%in%.))}}
#'
#' ## Frequently-used terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", keywords=~"frequently-used"%in%., display.keywords = subset(ergm::ergm_keyword(), popular)$name))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", keywords=~"frequently-used"%in%., display.keywords = subset(ergm::ergm_keyword(), popular)$name))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", keywords=~"frequently-used"%in%., display.keywords = subset(ergm::ergm_keyword(), popular)$name))}}
#'
#' ## Operator terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm", keywords=~"operator"%in%., display.keywords = subset(ergm::ergm_keyword(), popular & name!="operator")$name))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm", keywords=~"operator"%in%., display.keywords = subset(ergm::ergm_keyword(), popular & name!="operator")$name))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm", keywords=~"operator"%in%., display.keywords = subset(ergm::ergm_keyword(), popular & name!="operator")$name))}}
#' 
#' ## All terms
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmTerm"))}}
#' 
#' ## Terms by keywords
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmTerm"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmTerm"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmTerm"))}}
#'
#' @seealso [`ergm`][ergm-package] package, [`search.ergmTerms`], [`ergm`], [`network`], [`%v%`], [`%n%`]
#'
#' @references 
#' - Krivitsky P. N., Hunter D. R., Morris M., Klumb
#' C. (2021). "ergm 4.0: New features and improvements."
#' arXiv:2106.04997. \url{https://arxiv.org/abs/2106.04997}
#' 
#' - Bomiriya, R. P, Bansal, S., and Hunter, D. R. (2014).  Modeling
#' Homophily in ERGMs for Bipartite Networks.  Submitted.
#' 
#' - Butts, CT.  (2008).  "A Relational Event Framework for Social
#' Action." *Sociological Methodology,* 38(1).
#' 
#' - Davis, J.A. and Leinhardt, S.  (1972).  The Structure of Positive
#' Interpersonal Relations in Small Groups.  In J. Berger (Ed.),
#' *Sociological Theories in Progress, Volume 2*, 218--251.  Boston:
#' Houghton Mifflin.
#' 
#' - Holland, P. W. and S. Leinhardt (1981). An exponential family of
#' probability distributions for directed graphs.  *Journal of the
#' American Statistical Association*, 76: 33--50.
#' 
#' - Hunter, D. R. and M. S. Handcock (2006). Inference in curved
#' exponential family models for networks. *Journal of Computational and
#' Graphical Statistics*, 15: 565--583.
#' 
#' - Hunter, D. R. (2007). Curved exponential family models for social
#' networks. *Social Networks*, 29: 216--230.
#' 
#' - Krackhardt, D. and Handcock, M. S. (2007).  Heider versus Simmel:
#' Emergent Features in Dynamic Structures. *Lecture Notes in Computer
#' Science*, 4503, 14--27.
#' 
#' - Krivitsky P. N. (2012). Exponential-Family Random Graph Models for
#' Valued Networks. *Electronic Journal of Statistics*, 2012, 6,
#' 1100-1128. \doi{10.1214/12-EJS696}
#' 
#' - Robins, G; Pattison, P; and Wang, P.  (2009).  "Closure,
#' Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models
#' for Directed Social Networks." *Social Networks,* 31:105-117.
#' 
#' - Snijders T. A. B., G. G. van de Bunt, and C. E. G. Steglich.
#' Introduction to Stochastic Actor-Based Models for Network Dynamics.
#' *Social Networks*, 2010, 32(1), 44-60. \doi{10.1016/j.socnet.2009.02.004}
#' 
#' - Morris M, Handcock MS, and Hunter DR. Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 2008, 24(4), 1-24.
#' \doi{10.18637/jss.v024.i04}
#' 
#' - Snijders, T. A. B., P. E. Pattison, G. L. Robins, and M. S. Handcock
#' (2006). New specifications for exponential random graph models,
#' *Sociological Methodology*, 36(1): 99-153.
#' 
#' @keywords models
#' 
#' @examples
#' \dontrun{
#' ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#' 
#' ergm(molecule ~ edges + kstar(2:3) + triangle
#'                       + nodematch("atomic type",diff=TRUE)
#'                       + triangle + absdiff("atomic type"))
#' }
NULL
#TODO: Write a valued example.


#' Sample Space Constraints for Exponential-Family Random Graph Models
#'
#' @name ergmConstraint
#' @aliases ergm-constraints constraints-ergm ergm.constraints constraints.ergm
#' @description This page describes how to specify the constraints on the network sample space (the set of possible networks \eqn{Y}, the set of networks \eqn{y} for which \eqn{h(y)>0}) and sometimes the baseline weights \eqn{h(y)} to functions in the [`ergm`][ergm-package]
#' package. It also provides an indexed list of the constraints visible to the \CRANpkg{ergm}'s API. Constraints can also be searched via [`search.ergmConstraints`], and help for an individual constraint can be obtained with `ergmConstraint?<constraint>` or `help("<constraint>-ergmConstraint")`.
#'
#' @section Specifying constraints:
#' \ERGMspec
#' Constraints typically affect \eqn{Y}, or, equivalently, set \eqn{h(y)=0} for some \eqn{y}, but some (\dQuote{soft} constraints) set 
#' \eqn{h(y)} to values other than 0 and 1.
#' 
#' A constraints formula is a one- or two-sided formula whose left-hand side is
#' an optional direct selection of the `InitErgmProposal` function and
#' whose right-hand side is a series of one or more terms separated by
#' `"+"` and `"-"` operators, specifying the constraint.
#' 
#' The sample space (over and above the reference distribution) is determined
#' by iterating over the constraints terms from left to right, each term
#' updating it as follows: 
#' - If the constraint introduces complex
#' dependence structure (e.g., constrains degree or number of edges in the
#' network), then this constraint always restricts the sample space. It may
#' only have a `"+"` sign.
#' 
#' - If the constraint only restricts the set of dyads that may vary in the
#' sample space (e.g., block-diagonal structure or fixing specific dyads at
#' specific values) and has a `"+"` sign, the set of dyads that may
#' vary is restricted to those that may vary according to this constraint
#' *and* all the constraints to date.
#' 
#' - If the constraint only restricts the set of dyads that may vary in the
#' sample space but has a `"-"` sign, the set of dyads that may
#' vary is expanded to those that may vary according to this constraint
#' *or* all the constraints up to date.
#' 
#' For example, a constraints formula `~a-b+c-d` with all constraints
#' dyadic will allow dyads permitted by either `a` or `b` but only if they are
#' also permitted by `c`; as well as all dyads permitted by `d`. If `A`, `B`,
#' `C`, and `D` were logical matrices, the matrix of variable dyads would be
#' equal to `((A|B)&C)|D`.
#' 
#' Terms with a positive sign can be viewed as "adding" a constraint
#' while those with a negative sign can be viewed as "relaxing" a constraint.
#'
#' \subsection{Inheriting constraints from LHS [`network`]}{
#'
#' By default, [`%ergmlhs%`] attributes `constraints` or
#' `constraints.obs` (depending on which constraint) attached to the
#' LHS of the model formula or the `basis=` argument will be added in
#' front of the specified constraints formula. This is the desired
#' behaviour most of the time, since those constraints are usually
#' determined by how the network was constructed (e.g., structural
#' zeros in a block-diagonal network).
#'
#' For those situations in which this is not the desired behavior, a
#' `.` term (with a positive sign or no sign at all) can be used to
#' manually set the position of the inherited constraints in the
#' formula, and a `-.` (minus-dot) term anywhere in the constraints
#' formula will suppress the inherited formula altogether.
#'
#' }
#'
#' @section Constraints visible to the package:
#'
#' \ergmCSS
#'
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmConstraint"))}}
#'
#' ## All constraints
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmConstraint"))}}
#' 
#' ## Constraints by keywords
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmConstraint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmConstraint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmConstraint"))}}
#'
#' @references
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \CRANpkg{statnet} Tutorial. *Journal of Statistical Software*, 24(8).
#' \doi{10.18637/jss.v024.i08}
#'
#' - Hunter, D. R. and Handcock, M. S. (2006) *Inference in curved
#' exponential family models for networks*, Journal of Computational and
#' Graphical Statistics.
#'
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \CRANpkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  *Journal of Statistical Software*, 24(3).
#' \doi{10.18637/jss.v024.i03}
#'
#' - Karwa V, Krivitsky PN, and Slavkovi\'c AB (2016). Sharing Social Network
#' Data: Differentially Private Estimation of Exponential-Family Random Graph
#' Models. *Journal of the Royal Statistical Society, Series C*, 66(3):
#' 481-500. \doi{10.1111/rssc.12185}
#'
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' - Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 24(4). \doi{10.18637/jss.v024.i04}
#' @keywords models
NULL

#' MCMC Hints for Exponential-Family Random Graph Models
#'
#' @name ergmHint
#' @aliases ergm-hints hints-ergm ergm.hints hints.ergm hints
#' @description This page describes how to provide to the
#'   \CRANpkg{ergm}'s MCMC algorithms information about the sample space. Hints can also be searched via [`search.ergmHints`], and help for an individual hint can be obtained with `ergmHint?<hint>` or `help("<hint>-ergmHint")`.
#'
#' @section \dQuote{Hints} for MCMC:
#' \ERGMspec
#'
#' \newcommand{\Hint}{\dQuote{Hint}}
#' \newcommand{\hint}{\dQuote{hint}}
#' \newcommand{\Hints}{\dQuote{Hints}}
#' \newcommand{\hints}{\dQuote{hints}}
#'
#' It is often the case that there is additional information available
#' about the distribution of networks being modelled. For example, you
#' may be aware that the network is sparse or that there are strata
#' among the dyads. \Hints, typically passed on the right-hand side of `MCMC.prop`
#' and `obs.MCMC.prop` arguments to [control.ergm()],
#' [control.simulate.ergm()], and others, allow this information to be
#' provided. By default, hint [`sparse`][sparse-ergmHint] is in
#' effect.
#'
#' Unlike constraints, model terms, and reference distributions,
#' \hints{} do not affect the specification of the model. That is,
#' regardless of what \hints{} may or may not be in effect, the sample
#' space and the probabilities within it are the same. However,
#' \hints{} may affect the MCMC proposal distribution used by the
#' samplers.
#'
#' Note that not all proposals support all \hints: and if the most
#' suitable proposal available cannot incorporate a particular \hint,
#' a warning message will be printed.
#'
#' \Hints{} use the same underlying API as constraints, and, if present,
#' [`%ergmlhs%`] attributes `constraints` and `constraints.obs` will
#' be substituted in its place.
#'
#' @section Hints available to the package:
#'
#' The following hints are known to \CRANpkg{ergm} at this time:
#'
#' \ergmCSS
#'
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmHint"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmHint"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmHint"))}}
#' 
#' @references
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A
#' \CRANpkg{statnet} Tutorial. *Journal of Statistical Software*, 24(8).
#' \doi{10.18637/jss.v024.i08}
#' 
#' - Hunter, D. R. and Handcock, M. S. (2006) *Inference in curved
#' exponential family models for networks*, Journal of Computational and
#' Graphical Statistics.
#' 
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
#' \CRANpkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  *Journal of Statistical Software*, 24(3).
#' \doi{10.18637/jss.v024.i03}
#' 
#' - Karwa V, Krivitsky PN, and Slavkovi\'c AB (2016). Sharing Social Network
#' Data: Differentially Private Estimation of Exponential-Family Random Graph
#' Models. *Journal of the Royal Statistical Society, Series C*, 66(3):
#' 481-500. \doi{10.1111/rssc.12185}
#' 
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#' 
#' - Morris M, Handcock MS, Hunter DR (2008).  Specification of
#' Exponential-Family Random Graph Models: Terms and Computational Aspects.
#' *Journal of Statistical Software*, 24(4). \doi{10.18637/jss.v024.i04}
#' @keywords models
NULL


#' Reference Measures for Exponential-Family Random Graph Models
#'
#' @name ergmReference
#' @aliases ergm-references references-ergm ergm.references references.ergm
#' @description This page describes how to specify the reference measures (baseline distributions)
#' (the set of possible networks \eqn{Y} and the baseline weights \eqn{h(y)} to functions in the [`ergm`][ergm-package]
#' package. It also provides an indexed list of the references visible to the \CRANpkg{ergm}'s API. References can also be searched via [search.ergmReferences()], and help for an individual reference can be obtained with `ergmReference?<reference>` or `help("<reference>-ergmReference")`.
#'
#' @section Specifying reference measures:
#' \ERGMspec
#'
#' The reference measure \eqn{(Y,h(y))} is specified on the right-hand side of a one-sided formula passed
#' typically as the `reference` argument.
#'
#' @section Reference measures visible to the package:
#'
#' \ergmCSS
#'
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexLatex(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexText(ergm:::.buildTermsDataframe("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatIndexHtml(ergm:::.buildTermsDataframe("ergmReference"))}}
#'
#' ## All references
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixLatex(ergm:::.termMatrix("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixText(ergm:::.termMatrix("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatMatrixHtml(ergm:::.termMatrix("ergmReference"))}}
#' 
#' ## References by keywords
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocLatex(ergm:::.termToc("ergmReference"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocText(ergm:::.termToc("ergmReference"))}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatTocHtml(ergm:::.termToc("ergmReference"))}}
#'
#' @seealso [`ergm`][ergm-package], [`network`], \CRANpkg{sna}, [`summary.ergm`], [`print.ergm`], `\%v\%`, `\%n\%`
#' 
#' @references
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b). \CRANpkg{ergm}:
#' A Package to Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks. *Journal of Statistical Software*, 24(3).
#' \doi{10.18637/jss.v024.i03}
#' 
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' @keywords models
NULL

#' Metropolis-Hastings Proposal Methods for ERGM MCMC
#'
#' @name ergmProposal
#' @aliases ergm-proposals proposals-ergm ergm.proposals
#'   proposals.ergm InitErgmProposal InitWtErgmProposal
#' @description This page describes the low-level Metropolis--Hastings
#'   (MH) proposal algorithms. They are rarely invoked directly by the
#'   user but are rather selected based on the provided [sample space
#'   constraints][ergmConstraint] and [hints about the network
#'   process][ergmHint].  They can also be searched via
#'   [`search.ergmProposals`], and help for an individual proposal can
#'   be obtained with `ergmProposal?<proposal>` or
#'   `help("<proposal>-ergmProposal")`.
#'
#' @details [`ergm`] uses a Metropolis-Hastings (MH) algorithm to
#'   control the behavior of the Markov Chain Monte Carlo (MCMC) for
#'   sampling networks.  The MCMC chain is intended to step around the
#'   sample space of possible networks, generating a network at
#'   regular intervals to evaluate the statistics in the model.  For
#'   each MCMC step, one or more toggles are proposed to change the
#'   dyads to the opposite value. The probability of accepting the
#'   proposed change is determined by the MH acceptance ratio.  The
#'   role of the different MH methods implemented in
#'   [ergm()] is to vary how the sets of dyads are selected
#'   for toggle proposals.  This is used in some cases to improve the
#'   performance (speed and mixing) of the algorithm, and in other
#'   cases to constrain the sample space.
#'
#' @section Proposals available to the package:
#'
#' \ergmCSS
#'
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsHtml(ergm:::.buildProposalsList(), keepProposal=TRUE)}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsText(ergm:::.buildProposalsList(), keepProposal=TRUE)}}
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsLatex(ergm:::.buildProposalsList(), keepProposal=TRUE)}}
#'
#' Note that [`.dyads`][.dyads-ergmConstraint] is a meta-constraint, indicating that the proposal supports an arbitrary dyad-level constraint combination.
#'
#' @seealso [`ergm`][ergm-package] package, [`ergm`], [`ergmConstraint`], [`ergmHint`], [`ergm_proposal`]
#'
#' @references
#' - Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).  A \CRANpkg{statnet} Tutorial.
#' *Journal of Statistical Software*, 24(8). \doi{10.18637/jss.v024.i08}
#'
#' - Hunter, D. R. and Handcock, M. S. (2006) Inference in curved exponential family models for networks.
#' *Journal of Computational and Graphical Statistics*.
#'
#' - Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b). \CRANpkg{ergm}:
#' A Package to Fit, Simulate and Diagnose Exponential-Family Models for
#' Networks. *Journal of Statistical Software*, 24(3).
#' \doi{10.18637/jss.v024.i03}
#'
#' - Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
#' Networks. *Electronic Journal of Statistics*, 2012, 6, 1100-1128.
#' \doi{10.1214/12-EJS696}
#'
#' - Morris M, Handcock MS, Hunter DR (2008). Specification of Exponential-Family Random Graph Models:
#' Terms and Computational Aspects. *Journal of Statistical Software*, 24(4).
#' \doi{10.18637/jss.v024.i04}
#'
#' @keywords models
NULL
