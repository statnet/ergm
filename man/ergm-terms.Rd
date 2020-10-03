%  File man/ergm-terms.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2003-2020 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{InitErgmTerm}

\alias{absdiff}
\alias{absdiffcat}
\alias{altkstar}
\alias{asymmetric}
\alias{atleast}
\alias{atmost}
\alias{b1concurrent}
\alias{b1cov}
\alias{b1degrange}
\alias{b1degree}
\alias{b1dsp}
\alias{b1factor}
\alias{b1mindegree}
\alias{b1nodematch}
\alias{b1sociality}
\alias{b1star}
\alias{b1starmix}
\alias{b1twostar}
\alias{b2concurrent}
\alias{b2cov}
\alias{b2degrange}
\alias{b2degree}
\alias{b2dsp}
\alias{b2factor}
\alias{b2mindegree}
\alias{b2nodematch}
\alias{b2sociality}
\alias{b2star}
\alias{b2starmix}
\alias{b2twostar}
\alias{balance}
\alias{coincidence}
\alias{concurrent}
\alias{concurrentties}
\alias{ctriad}
\alias{ctriple}
\alias{cycle}
\alias{cyclicalties}
\alias{cyclicalweights}
\alias{ddsp}
\alias{degcor}
\alias{degcrossprod}
\alias{degrange}
\alias{degree}
\alias{degree1.5}
\alias{degreepopularity}
\alias{diff}
\alias{density}
\alias{desp}
\alias{dgwdsp}
\alias{dgwesp}
\alias{dgwnsp}
\alias{dnsp}
\alias{dsp}
\alias{dyadcov}
\alias{edgecov}
\alias{edges-term}
\alias{esp}
\alias{equalto}
\alias{greaterthan}
\alias{gwb1degree}
\alias{gwb1dsp}
\alias{gwb2degree}
\alias{gwb2dsp}
\alias{gwdegree}
\alias{gwdsp}
\alias{gwesp}
\alias{gwidegree}
\alias{gwnsp}
\alias{gwodegree}
\alias{hamming-term}
\alias{hammingmix}
\alias{idegrange}
\alias{idegree}
\alias{idegree1.5}
\alias{idegreepopularity}
\alias{ininterval}
\alias{intransitive}
\alias{isolatededges}
\alias{isolates}
\alias{istar}
\alias{kstar}
\alias{localtriangle}
\alias{m2star}
\alias{match}
\alias{mm}
\alias{meandeg}
\alias{mutual}
\alias{nearsimmelian}
\alias{nodecov}
\alias{nodecovar}
\alias{nodefactor}
\alias{nodeicov}
\alias{nodeicovar}
\alias{nodeifactor}
\alias{nodeisqrtcovar}
\alias{nodemain}
\alias{nodematch}
\alias{nodemix}
\alias{nodeocov}
\alias{nodeocovar}
\alias{nodeofactor}
\alias{nodeosqrtcovar}
\alias{nodesqrtcovar}
\alias{nonzero}
\alias{nsp}
\alias{odegrange}
\alias{odegree}
\alias{odegree1.5}
\alias{odegreepopularity}
\alias{opentriad}
\alias{ostar}
\alias{receiver}
\alias{sender}
\alias{simmelian}
\alias{simmelianties}
\alias{smalldiff}
\alias{smallerthan}
\alias{sociality}
\alias{sum}
\alias{threepath}
\alias{threetrail}
\alias{transitive}
\alias{transitiveties}
\alias{transitiveweights}
\alias{triadcensus}
\alias{triangle}
\alias{triangles}
\alias{tripercent}
\alias{ttriad}
\alias{ttriple}
\alias{twopath}
\docType{methods}

% The following macros define commonly repeated parts of the term description. It appears to be a limitation of the format that the whole macro must fit on one line.

% Versioning

\newcommand{\ThreeNineFour}{Note that \code{ergm} versions 3.9.4 and earlier used different arguments for this term. See the above section on versioning for invoking the old behavior.}

% Directedness

\newcommand{\Bipartite}{This term can only be used with undirected bipartite networks.}

\newcommand{\Undirected}{This term can only be used with undirected networks.}

\newcommand{\Directed}{This term can only be used with directed networks.}

% levels= arguments

\newcommand{\LevelsExplain}{The optional \code{levels} argument controls which levels of the attribute should be included and which should be excluded. (See \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details.) For example, if the \dQuote{fruit} attribute has levels \dQuote{orange}, \dQuote{apple}, \dQuote{banana}, and \dQuote{pear}, then to include just two levels, one for \dQuote{apple} and one for \dQuote{pear}, use any of \code{#1("fruit", levels=-(2:3))}, \code{#1("fruit", levels=c(1,4))}, and \code{#1("fruit", levels=c("apple", "pear"))}. Note: if you are using numeric values to specify the levels of a character variable, the levels will correspond to the alphabetically sorted character levels.}

\newcommand{\LevelsAll}{By default, all levels are included.}

\newcommand{\LevelsNotFirst}{To include all attribute values is usually not a good idea, because the sum of all such statistics equals the number of edges and hence a linear dependency would arise in any model also including \code{edges}. The default, \code{levels=-1}, is therefore to omit the first (in lexicographic order) attribute level. To include all levels, pass either \code{levels=TRUE} (i.e., keep all levels) or \code{levels=NULL} (i.e., do not filter levels).}

% Deprecation of base= and keep= arguments

\newcommand{\BaseDep}{The argument \code{base} is retained for backwards compatibility and may be removed in a future version. When both \code{base} and \code{levels} are passed, \code{levels} overrides \code{base}.}

\newcommand{\BaseDepNodes}{The argument \code{base} is retained for backwards compatibility and may be removed in a future version. When both \code{base} and \code{nodes} are passed, \code{nodes} overrides \code{base}.}

\newcommand{\BaseDepII}{The argument \code{base} is retained for backwards compatibility and may be removed in a future version. When both \code{base} and \code{levels2} are passed, \code{levels2} overrides \code{base}.}

\newcommand{\KeepDep}{The argument \code{keep} is retained for backwards compatibility and may be removed in a future version. When both \code{keep} and \code{levels} are passed, \code{levels} overrides \code{keep}.}

% Shared partner types

\newcommand{\SPTypes}{While there is only one shared partner configuration in the undirected case, nine distinct configurations are possible for directed graphs, selected using the \code{type} argument. Currently, terms may be defined with respect to five of these configurations; they are defined here as follows (using terminology from Butts (2008) and the \code{relevent} package): \describe{ \item{Outgoing Two-path (\code{"OTP"})}{vertex \eqn{k} is an OTP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k \to j}{i->k->j}. Also known as "transitive shared partner".}  \item{Incoming Two-path (\code{"ITP"})}{vertex \eqn{k} is an ITP shared partner of ordered pair \eqn{(i,j)} iff \eqn{j \to k \to i}{j->k->i}. Also known as "cyclical shared partner"} \item{Reciprocated Two-path (\code{"RTP"})}{vertex \eqn{k} is an RTP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \leftrightarrow k \leftrightarrow j}{i<->k<->j}.} \item{Outgoing Shared Partner (\code{"OSP"})}{vertex \eqn{k} is an OSP shared partner of ordered pair \eqn{(i,j)} iff \eqn{i \to k, j \to k}{i->k, j->k}.}  \item{Incoming Shared Partner (\code{"ISP"})}{vertex \eqn{k} is an ISP shared partner of ordered pair \eqn{(i,j)} iff \eqn{k \to i, k \to j}{k->i, k->j}.}} By default, outgoing two-paths (\code{"OTP"}) are calculated. Note that Robins et al. (2009) define closely related statistics to several of the above, using slightly different terminology.}

\newcommand{\GWCutoff}{The optional argument \code{cutoff} sets the number of underlying #1 terms to use in computing the statistics when \code{fixed=FALSE}, in order to reduce the computational burden.  Its default value can also be controlled by the \code{gw.cutoff} term option control parameter. (See \code{\link{control.ergm}}.)}

\newcommand{\SPtoDSP}{For directed networks, only outgoing two-path ("OTP") shared partners are counted.  In other words, for a #2 in a directed graph, the number of shared partners counted by \code{#1} is the number of nodes \code{k} that have edges \code{i -> k -> j}.  (These may also be called homogeneous shared partners.)  To count other types of shared partners instead, see \code{#3}.}

\title{Terms used in Exponential Family Random Graph Models}
\description{
  The function \code{\link{ergm}} is used to fit exponential
  random graph models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is
  \deqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}{h(y) exp\{eta(theta).g(y)\} / c(theta),} where
  \eqn{h(y)} is the reference measure (for valued network models),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  The network statistics \eqn{g(y)} are entered as terms in the
  function call to \code{\link{ergm}}.  This page describes
  the possible terms (and hence network statistics)
  included in \code{\link[=ergm-package]{ergm}} package.

  A cross-referenced HTML version of the term documentation is
  available via \code{vignette('ergm-term-crossRef')} and terms
  can also be searched via \code{\link{search.ergmTerms}}.
}

\section{Specifying models}{

  Terms to \code{\link{ergm}} are specified by a formula to represent the network and
  network statistics. This is done via a \code{formula}, that is,
  an \R formula object, of the form
  \code{y ~ <term 1> + <term 2> ...},
  where \code{y} is a network object or a matrix that can be coerced to a
  network
  object, and \code{<term 1>}, \code{<term 2>}, etc, are each terms chosen
  from the list given below.
  To create a network object in \R, use the \code{\link[network]{network}} function,
  then add nodal attributes to it using the \code{\%v\%} operator if necessary.

\subsection{Binary and valued ERGM terms}{

  \code{\link[=ergm-package]{ergm}} functions such as \code{\link{ergm}}
  and \code{\link[=simulate.formula]{simulate}} (for ERGMs) may operate
  in two modes: binary and weighted/valued, with the latter activated by
  passing a non-NULL value as the \code{response} argument, giving the
  edge attribute name to be modeled/simulated.

  Binary ERGM statistics cannot be used in valued mode and vice versa. However,
  a substantial number of binary ERGM statistics --- particularly the
  ones with dyadic indepenence --- have simple generalizations to valued
  ERGMs, and have been adapted in
  \code{\link[=ergm-package]{ergm}}. They have the same form as their
  binary ERGM counterparts, with an additional argument: \code{form}, which, at
  this time, has two possible values: \code{"sum"} (the default) and
  \code{"nonzero"}. The former creates a statistic of the form
  \eqn{\sum_{i,j} x_{i,j} y_{i,j}}, where \eqn{y_{i,j}} is the value of
  dyad \eqn{(i,j)} and \eqn{x_{i,j}} is the term's covariate associated
  with it. The latter computes the binary version, with the edge
  considered to be present if its value is not 0.

  Valued version of some binary ERGM terms have an argument
  \code{threshold}, which sets the value above which a dyad is conidered
  to have a tie. (Value less than or equal to \code{threshold} is
  considered a nontie.)
}

\subsection{Nodal attribute levels and indices}{

  Terms taking a categorical nodal covariate also take the \code{levels}
  argument.  (There are analogous \code{b1levels} and \code{b2levels} arguments
  for some terms that apply to bipartite networks, and the \code{levels2} argument
  for mixing terms.)  The \code{levels} argument can be used to control the set and the ordering of
  attribute levels.

  Terms that allow the selection of nodes do so with the \code{nodes}
  argument, which is interpreted in the same way as the \code{levels}
  argument, where the categories are the relevant nodal indices
  themselves.

  Both \code{levels} and \code{nodes} use the new level
  selection UI. (See \link[=nodal_attributes]{Specifying Vertex Attributes and
    Levels} for details.)

  \subsection{Legacy arguments}{

    The legacy \code{base} and \code{keep} arguments are
    deprecated as of version 3.10, and replaced by the
    \code{levels} UI. The \code{levels}
    argument provides consistent and flexible mechanisms for
    specifying which attribute levels to exclude (previously
    handled by \code{base}) and include (previously
    handled by \code{keep}).  If \code{levels} or
    \code{nodes} argument is given, then \code{base} and \code{keep}
    arguments are ignored. The legacy arguments
    will most likely be removed in a future version.

    Note that this exact behavior is new in version 3.10, and it differs
    slightly from older versions: previously if both \code{levels} and
    \code{base}/\code{keep} were given, \code{levels} argument was
    applied first and then applied the \code{base}/\code{keep}
    argument. Since version 3.10, \code{base}/\code{keep} would be
    ignored, even if old term behavior is invoked (as described in the
    next section).

  }

}

\subsection{Term versioning}{

  When a term's behavior has changed from prior version, it is often possible to invoke the old behavior by setting and/or passing a \code{version} term option, giving the verison (constructed by \code{\link{as.package_version}}) desired.
}

\subsection{Custom \code{ergm} terms}{

  Users and other packages may build custom terms, and package
  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}}
  provides tools for implementing them.

  The current recommendation for any package implementing additional
  terms is to create a help file with a name or
  alias \code{ergm-terms}, so that \code{help("ergm-terms")} will
  list ERGM terms available from all loaded packages.
}

}

\section{Terms included in the \code{\link[=ergm-package]{ergm}} package}{

  As noted above, a cross-referenced HTML version of the
  term documentation is
  available via \code{vignette('ergm-term-crossRef')} and terms
  can also be searched via \code{\link{search.ergmTerms}}.
  \describe{
    %%%%%%%%%%%%%%%%%%%%%%555
    % This section below is specially formatted so that it can be parsed to build additional
    % documentation about ergm terms. (parsing functions are in build_term_index.R)
    % The list of term tags must appear on the same line as the term name to be parsed correctly
    %%beginTerms%%
    \item{\code{absdiff(attr, pow=1)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (quantitative nodal attribute),
      \code{absdiff(attr, pow=1, form ="sum")}  (valued)  (dyad-independent) (directed)  (undirected) (quantitative nodal attribute)}{\emph{Absolute difference:}
      The \code{attr} argument specifies a quantitative attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds one network statistic to the model equaling the
      sum of \code{abs(attr[i]-attr[j])^pow} for all edges (i,j)
      in the network.

      \ThreeNineFour}

    \item{\code{absdiffcat(attr, base=NULL, levels=NULL)}  (binary) (dyad-independent) (directed) (undirected) (categorical nodal attribute),
      \code{absdiffcat(attr, base=NULL, levels=NULL, form="sum")}  (valued)  (dyad-independent)   (directed)  (undirected) (categorical nodal attribute)}{\emph{Categorical absolute
	difference:}
      The \code{attr} argument specifies a quantitative attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term
      adds one statistic for every possible nonzero distinct value of
      \code{abs(attr[i]-attr[j])} in the network; the value of each such
      statistic is the number of edges in the network with the corresponding
      absolute difference. The optional argument \code{levels} specifies which
      nonzero differences to include in or exclude from the model (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). For example, if the possible values of
      \code{abs(attr[i]-attr[j])} are 0, 0.5, 3, 3.5, and 10, then to omit
      0.5 and 10 one could set \code{levels=2:3} (we wish to retain the second and third nonzero difference,
      when differences are sorted in increasing order). Note that this term should
      generally be used only when the quantitative attribute has a limited number
      of possible values; an example is the \code{"Grade"} attribute of the
      \code{\link{faux.mesa.high}} or \code{\link{faux.magnolia.high}}
      datasets.

      \BaseDep}

    \item{\code{altkstar(lambda, fixed=FALSE)}  (binary) (undirected) (curved) (categorical nodal attribute)
    }{\emph{Alternating k-star:}
      This term adds one network statistic to the model equal to a weighted
      alternating sequence of k-star statistics with weight parameter
      \code{lambda}. This is the version given in Snijders et al. (2006). The
      \code{gwdegree} and \code{altkstar} produce mathematically equivalent
      models, as long as they are used together with the \code{edges} (or
      \code{kstar(1)}) term, yet the interpretation of the \code{gwdegree}
      parameters is slightly more straightforward than the interpretation of the
      \code{altkstar} parameters. For this reason, we recommend the use of the
      \code{gwdegree} instead of \code{altkstar}. See Section 3 and especially
      equation (13) of Hunter (2007) for details. The optional argument
      \code{fixed} indicates whether the \code{decay} parameter is fixed at the given value, or is to be
      fit as a curved exponential family model (see Hunter and Handcock, 2006).
      The default is \code{FALSE}, which means the scale parameter is not fixed
      and thus the model is a CEF model. \Undirected}

    \item{\code{asymmetric(attr=NULL, diff=FALSE, keep=NULL, levels=NULL)}  (binary)  (directed)  (dyad-independent) (triad-related)}{\emph{Asymmetric
	dyads:} This term adds one network statistic to the model equal to the
      number of pairs of actors for which exactly one of
      \eqn{(i{\rightarrow}j)}{(i,j)} or \eqn{(j{\rightarrow}i)}{(j,i)} exists.
      \Directed The optional argument
      \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  If \code{attr} is specified, only asymmetric pairs that match on the
      vertex attribute \code{attr} are counted. The optional modifiers
      \code{diff} and \code{levels} are used in the same way as for the \code{nodematch} term; refer
      to this term for details and an example.

      \KeepDep}

    \item{\code{atleast(threshold=0)}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Number of dyads with values
	greater than or equal to a threshold} Adds one statistic
      equaling to the number of dyads whose values equal or exceed
      \code{threshold}.}

    \item{\code{atmost(threshold=0)}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Number of dyads with values
	less than or equal to a threshold} Adds one statistic
      equaling to the number of dyads whose values equal or are exceeded by
      \code{threshold}.}

    \item{\code{b1concurrent(by=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{Concurrent node count for the
	first mode in a bipartite (aka two-mode) network:} This term adds one
      network statistic to the model, equal to the number of nodes in the first
      mode of the network with degree 2 or higher. The first mode of a bipartite
      network object is sometimes known as the "actor" mode. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details);
      it functions just like the \code{by} argument of the \code{b1degree} term.
      Without the optional argument, this statistic is equivalent to \code{b1mindegree(2)}.
      This term can only be
      used with undirected bipartite networks.}

    \item{\code{b1cov(attr)}  (binary) (undirected) (bipartite) (dyad-independent) (quantitative nodalattribute) (frequently-used),
      \code{b1cov(attr, form="sum")}  (valued)  (undirected) (bipartite) (dyad-independent) (quantitative nodal attribute) (frequently-used)
    }{\emph{Main effect of a covariate for
	the first mode in a bipartite (aka two-mode) network:}
      The \code{attr} argument specifies one or more quantitative attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
      value of \code{attr(i)} for all edges \eqn{(i,j)} in the network. This
      term may only be used with bipartite networks. For categorical attributes,
      see \code{b1factor}.

      \ThreeNineFour}

    \item{\code{b1degrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (bipartite)  (undirected)
    }{\emph{Degree range for the first mode in a bipartite
	(a.k.a. two-mode) network:}
      The \code{from} and \code{to} arguments are vectors of distinct
      integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals the number of nodes of the first mode
      ("actors") in the network of degree
      greater than or equal to
      \code{from[i]} but strictly less than \code{to[i]}, i.e. with
      edge count
      in semiopen interval \code{[from,to)}. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree range
      statistics are calculated for nodes having each separate
      value of the attribute.

      This term can only be used with bipartite networks; for directed networks
      see \code{idegrange} and \code{odegrange}. For undirected networks,
      see \code{degrange}, and see \code{b2degrange}
      for degrees of the second mode ("events").}

    \item{\code{b1degree(d, by=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute) (frequently-used)
    }{\emph{Degree for the first mode in a
	bipartite (aka two-mode) network:} The \code{d} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{d}; the \eqn{i}th such statistic equals the number of
      nodes of degree \code{d[i]} in the first mode of a bipartite network, i.e.
      with exactly \code{d[i]} edges. The first mode of a bipartite network object
      is sometimes known as the "actor" mode. The optional argument \code{by} specifies
      a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified
      then each node's degree is tabulated only with other nodes having the same
      value of the \code{by} attribute.

      \Bipartite}

    \item{\code{b1dsp(d)}  (binary)  (bipartite) (undirected)}{\emph{Dyadwise shared partners for dyads in the first bipartition:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the number of dyads in the first bipartition with exactly
      \code{d[i]} shared partners.  (Those shared partners, of course, must be members 
      of the second bipartition.)  This term can only be used with bipartite networks.}

    \item{\code{b1factor(attr, base=1, levels=-1)}  (binary)  (bipartite) (undirected)  (dyad-independent) (frequently-used) (categorical nodal attribute),
      \code{b1factor(attr, base=1, levels=-1, form="sum")}  (valued)  (bipartite)  (undirected)  (dyad-independent) (frequently-used) (categorical nodal attribute)
    }{\emph{Factor attribute effect for
	the first mode in a bipartite (aka two-mode) network:}
      The \code{attr} argument specifies a categorical vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds
      multiple network statistics to the model, one for each of (a subset of) the
      unique values of the \code{attr} attribute. Each of these statistics
      gives the number of times a node with that attribute in the first mode of
      the network appears in an edge. The first mode of a bipartite network object
      is sometimes known as the "actor" mode.

      \LevelsExplain{b1factor}

      \LevelsNotFirst

      \BaseDep

      \Bipartite}

    \item{\code{b1mindegree(d)}  (binary)  (bipartite)  (undirected)}{\emph{Minimum degree for the first mode in a
	bipartite (aka two-mode) network:} The \code{d} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{d}; the \eqn{i}th such statistic equals the number of
      nodes in the first mode of a bipartite network with at least degree \code{d[i]}.
      The first mode of a bipartite network object is sometimes known as the "actor" mode.

      \Bipartite}

    \item{\code{b1nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1,
	byb2attr=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (dyad-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Nodal attribute-based homophily effect for
	the first mode in a bipartite (aka two-mode) network:}  This term is introduced
      in Bomiriya et al (2014).
      The \code{attr} argument specifies a categorical vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). Out of the two arguments
      (discount parameters) \code{alpha} and \code{beta}, both of which take values from [0,1],
      only one should be set at a time. If none is set to a value other than 1, this term will
      simply be a homophily based two-star statistic. This term adds one statistic to the model
      unless \code{diff} is set to \code{TRUE}, in which case the term adds multiple network
      statistics to the model, one for each of (a subset of) the unique values of the \code{attr}
      attribute. To include only the attribute values you wish, use the \code{levels} arguments.

      \KeepDep

      If an \code{alpha}
      discount parameter is used, each of these statistics gives the sum of
      the number of common second-mode nodes raised to the power \code{alpha} for each pair of
      first-mode nodes with that attribute. If a \code{beta} discount parameter is used, each
      of these statistics gives half the sum of the number of two-paths with two first-mode nodes
      with that attribute as the two ends of the two path raised to the power \code{beta} for each
      edge in the network. The \code{byb2attr} argument specifies a
      second mode categorical attribute. Setting this argument
      will separate the orginal statistics based on the values of the set second mode attribute---
      i.e. for example, if \code{diff} is \code{FALSE}, then the sum of all the statistics for
      each level of this second-mode attribute will be equal to the original \code{b1nodematch}
      statistic where \code{byb2attr} set to \code{NULL}. \Bipartite}

    \item{\code{b1sociality(nodes=-1)}  (binary) (bipartite) (undirected) (dyad-independent)
      ,
      \code{b1sociality(nodes=-1, form="sum")}  (valued) (bipartite) (undirected) (dyad-independent)
    }{\emph{Degree:}
      This term adds one network statistic for each node in the first bipartition, equal to the number of
      ties of that node.  By default, \code{nodes=-1} means that the statistic for the
      first node will be omitted, but this argument may be changed to control
      which statistics are included.  The \code{nodes} argument is interpreted using the new UI for level specification
      (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details), where both the attribute and the sorted
      unique values are the vector of vertex indices \code{1:nb1}, where
      \code{nb1} is the size of the first bipartition.  This term can only be used with
      bipartite networks.  For directed networks, see \code{sender} and
      \code{receiver}.  For unipartite networks, see \code{sociality}.}

    \item{\code{b1star(k, attr=NULL, levels=NULL)}  (binary) (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars for the first mode in a
	bipartite (aka two-mode) network:} The \code{k} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{k}. The \eqn{i}th such statistic counts the number of
      distinct \code{k[i]}-stars whose center node is in the first mode of the
      network. The first mode of a bipartite network object is sometimes known as
      the "actor" mode. A \eqn{k}-star is defined to be a center node \eqn{N} and
      a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the
      ties \eqn{\{N, O_i\}}{\{N, O[i]\}} exist for \eqn{i=1, \dots, k}. The optional argument
      \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified then the count is over
      the number of \eqn{k}-stars (with center node in the first mode) where all
      nodes have the same value of the attribute. This term can only be used for
      undirected bipartite networks. Note that \code{b1star(1)} is equal to
      \code{b2star(1)} and to \code{edges}.}

    \item{\code{b1starmix(k, attr, base=NULL, diff=TRUE)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Mixing
	matrix for k-stars centered on the first mode of a bipartite network:} Only
      a single value of \eqn{k} is allowed. This term counts all k-stars in which
      the b2 nodes (called events in some contexts) are homophilous in the sense
      that they all share the same value of \code{attr}. However, the b1 node
      (in some contexts, the actor) at the center of the k-star does NOT have to
      have the same value as the b2 nodes; indeed, the values taken by the b1
      nodes may be completely distinct from those of the b2 nodes, which allows
      for the use of this term in cases where there are two separate nodal
      attributes, one for the b1 nodes and another for the b2 nodes (in this case,
      however, these two attributes should be combined to form a single nodal
      attribute, \code{attr}). A different statistic is created for each
      value of \code{attr} seen in a b1 node, even if no k-stars are observed
      with this value. Whether a different statistic is created for each value
      seen in a b2 node depends on the value of the \code{diff} argument: When
      \code{diff=TRUE}, the default, a different statistic is created for each
      value and thus the behavior of this term is reminiscent of the
      \code{nodemix} term, from which it takes its name; when \code{diff=FALSE},
      all homophilous k-stars are counted together, though these k-stars are still
      categorized according to the value of the central b1 node.}

    \item{\code{b1twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)}  (binary) (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Two-star
	census for central nodes centered on the first mode of a bipartite network:}
      This term takes two nodal attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details), one for b1 nodes (actors in some
      contexts) and one for b2 nodes (events in some contexts). Only
      \code{b1attr} is required; if \code{b2attr} is not passed, it is
      assumed to be the same as \code{b1attr}. Assuming that there are
      \eqn{n_1} values of \code{b1attr} among the b1 nodes and \eqn{n_2}
      values of \code{b2attr} among the b2 nodes, then the total number of
      distinct categories of two stars according to these two attributes is
      \eqn{n_1(n_2)(n_2+1)/2}.  By default, this model term creates a distinct statistic
      counting each of these categories.  The \code{b1levels}, \code{b2levels}, \code{base}, and \code{levels2} arguments
      may be used to leave some of these categories out (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).

      \BaseDep \BaseDepII}

    \item{\code{b2concurrent(by=NULL)}  (binary)  (bipartite)  (undirected) (frequently-used)
    }{\emph{Concurrent node count for the
	second mode in a bipartite (aka two-mode) network:} This term adds one
      network statistic to the model, equal to the number of nodes in the second
      mode of the network with degree 2 or higher. The second mode of a bipartite
      network object is sometimes known as the "event" mode. The optional argument
      \code{by} specifies	a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details);
      it functions just like the \code{by} argument of the \code{b2degree} term.
      Without the optional argument, this statistic is equivalent to \code{b2mindegree(2)}.

      \Bipartite}

    \item{\code{b2cov(attr)}  (binary) (undirected) (bipartite) (dyad-independent) (quantitative nodal attribute) (frequently-used),
      \code{b2cov(attr, form="sum")}  (valued) (undirected) (bipartite) (dyad-independent) (quantitative nodal attribute) (frequently-used)
    }{\emph{Main effect of a covariate for
	the second mode in a bipartite (aka two-mode) network:}
      The \code{attr} argument specifies one or more quantitative attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
      value of \code{attr(j)} for all edges \eqn{(i,j)} in the network. This
      term may only be used with bipartite networks. For categorical attributes,
      see \code{b2factor}.

      \ThreeNineFour}

    \item{\code{b2degrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (bipartite)  (undirected)
    }{\emph{Degree range for the second mode in a bipartite
	(a.k.a. two-mode) network:}
      The \code{from} and \code{to} arguments are vectors of distinct
      integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals the number of nodes of the second mode
      ("events") in the network of degree
      greater than or equal to
      \code{from[i]} but strictly less than \code{to[i]}, i.e. with
      edge count
      in semiopen interval \code{[from,to)}. The optional argument
      \code{by} specifies
      a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree range
      statistics are calculated for nodes having each separate
      value of the attribute.

      This term can only be used with bipartite networks; for directed networks
      see \code{idegrange} and \code{odegrange}. For undirected networks,
      see \code{degrange}, and see \code{b1degrange}
      for degrees of the first mode ("actors").}

    \item{\code{b2degree(d, by=NULL)}  (binary)  (bipartite)  (undirected)  (categorical nodal attribute) (frequently-used)
    }{\emph{Degree for the second mode in a
	bipartite (aka two-mode) network:} The \code{d} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{d}; the \eqn{i}th such statistic equals the number of
      nodes of degree \code{d[i]} in the second mode of a bipartite network, i.e.
      with exactly \code{d[i]} edges. The second mode of a bipartite network
      object is sometimes known as the "event" mode. The optional term
      \code{by} specifies
      a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified
      then each node's degree is tabulated only with other nodes having the same
      value of the \code{by} attribute.

      \Bipartite}

    \item{\code{b2dsp(d)}  (binary)  (bipartite) (undirected)}{\emph{Dyadwise shared partners for dyads in the second bipartition:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the number of dyads in the second bipartition with exactly
      \code{d[i]} shared partners.  (Those shared partners, of course, must be members 
      of the first bipartition.)  This term can only be used with bipartite networks.}

    \item{\code{b2factor(attr, base=1, levels=-1)}  (binary)  (bipartite) (undirected)  (dyad-independent) (categorical nodal attribute) (frequently-used),
      \code{b2factor(attr, base=1, levels=-1, form="sum")}  (valued)  (bipartite)  (undirected)  (dyad-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Factor attribute effect for
	the second mode in a bipartite (aka two-mode) network :}
      The \code{attr} argument specifies
      a categorical vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds
      multiple network statistics to the model, one for each of (a subset of) the
      unique values of the \code{attr} attribute. Each of these statistics
      gives the number of times a node with that attribute in the second mode of
      the network appears in an edge. The second mode of a bipartite network
      object is sometimes known as the "event" mode.

      \LevelsExplain{b2factor}

      \LevelsNotFirst

      \BaseDep

      \Bipartite}

    \item{\code{b2mindegree(d)}  (binary) (bipartite)  (undirected)
    }{\emph{Minimum degree for the second mode in a
	bipartite (aka two-mode) network:} The \code{d} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{d}; the \eqn{i}th such statistic equals the number of
      nodes  in the second mode of a bipartite network with at least degree \code{d[i]}.
      The second mode of a bipartite network object is sometimes known as the "event" mode.

      \Bipartite}

    \item{\code{b2nodematch(attr, diff=FALSE, keep=NULL, alpha=1, beta=1,
	byb1attr=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (dyad-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Nodal attribute-based homophily effect for
	the second mode in a bipartite (aka two-mode) network:} This term is introduced
      in Bomiriya et al (2014).
      The \code{attr} argument specifies a categorical vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  Out of the two arguments
      (discount parameters) \code{alpha} and \code{beta}, both which takes values from [0,1],
      only one should be set at a time. If none is set to a value other than 1, this term will
      simply be a homophily based two-star statistic. This term adds one statistic to the model
      unless \code{diff} is set to \code{TRUE}, in which case the term adds multiple network
      statistics to the model, one for each of (a subset of) the unique values of the \code{attr}
      attribute. To include only the attribute values you wish, use the \code{levels}
      argument.

      \KeepDep

      If an \code{alpha}
      discount parameter is used, each of these statistics gives the sum of
      the number of common first-mode nodes raised to the power \code{alpha} for each pair of
      second-mode nodes with that attribute. If a \code{beta} discount parameter is used, each
      of these statistics gives half the sum of the number of two-paths with two second-mode nodes
      with that attribute as the two ends of the two path raised to the power \code{beta} for each
      edge in the network. The \code{byb1attr} argument specifies a
      first mode categorical attribute. Setting this argument
      will separate the orginal statistics based on the values of the set first mode attribute---
      i.e. for example, if \code{diff} is \code{FALSE}, then the sum of all the statistics for
      each level of this first-mode attribute will be equal to the original \code{b2nodematch}
      statistic where \code{byb1attr} set to \code{NULL}.

      \Bipartite}

    \item{\code{b2sociality(nodes=-1)}  (binary) (bipartite) (undirected) (dyad-independent)
      ,
      \code{b2sociality(nodes=-1, form="sum")}  (valued) (bipartite) (undirected) (dyad-independent)
    }{\emph{Degree:}
      This term adds one network statistic for each node in the second bipartition, equal to the number of
      ties of that node.  By default, \code{nodes=-1} means that the statistic for the
      first node (in the second bipartition) will be omitted, but this argument may be changed to control
      which statistics are included.  The \code{nodes} argument is interpreted using the new UI for level specification
      (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details), where both the attribute and the sorted
      unique values are the vector of vertex indices \code{(nb1 + 1):n}, where
      \code{nb1} is the size of the first bipartition and \code{n} is the total number of nodes in the network.  Thus \code{nodes=120} will include only the statistic
      for the 120th node in the second biparition, while \code{nodes=I(120)} will include only the statistic for the 120th node in the entire network. \Bipartite For directed networks, see \code{sender} and
      \code{receiver}.  For unipartite networks, see \code{sociality}.}

    \item{\code{b2star(k, attr=NULL, levels=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars for the second mode in a
	bipartite (aka two-mode) network:} The \code{k} argument is a vector of
      distinct integers. This term adds one network statistic to the model for
      each element in \code{k}. The \eqn{i}th such statistic counts the number of
      distinct \code{k[i]}-stars whose center node is in the second mode of the
      network. The second mode of a bipartite network object is sometimes known as
      the "event" mode. A \eqn{k}-star is defined to be a center node \eqn{N} and
      a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the
      ties \eqn{\{N, O_i\}}{\{N, O[i]\}} exist for \eqn{i=1, \dots, k}. The optional argument
      \code{attr} specifies
      a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified then the count is over
      the number of \eqn{k}-stars (with center node in the second mode) where all
      nodes have the same value of the attribute. This term can only be used for
      undirected bipartite networks. Note that \code{b2star(1)} is equal to
      \code{b1star(1)} and to \code{edges}.}

    \item{\code{b2starmix(k, attr, base=NULL, diff=TRUE)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Mixing
	matrix for k-stars centered on the second mode of a bipartite network:}
      This term is exactly the same as \code{b1starmix} except that the roles of
      b1 and b2 are reversed.}

    \item{\code{b2twostar(b1attr, b2attr, base=NULL, b1levels=NULL, b2levels=NULL, levels2=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Two-star
	census for central nodes centered on the second mode of a bipartite
	network:} This term is exactly the same as \code{b1twostar} except that the
      roles of b1 and b2 are reversed.}

    \item{\code{balance}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Balanced triads:}
      This term adds one network statistic to the model equal to the number of
      triads in the network that are balanced. The balanced triads are those of
      type \code{102} or \code{300} in the categorization of Davis and Leinhardt (1972). For details on the 16 possible triad types, see
      \code{?triad.classify} in the \code{{sna}} package. For an undirected
      network, the balanced triads are those with an odd number of ties (i.e., 1
      and 3).}

    \item{\code{coincidence(levels=NULL,active=0)}  (binary)  (bipartite)  (undirected)
    }{\emph{Coincident node count for the
	second mode in a bipartite (aka two-mode) network:} By default this term adds one
      network statistic to the model for each pair of nodes of mode two. It is
      equal to the number of (first mode) mutual partners of that pair.
      The first mode of a bipartite
      network object is sometimes known as the "actor" mode and the seconds as the "event" mode. So this is the number of actors going to both events in the pair. The optional argument
      \code{levels} specifies which pairs of nodes in mode two to include (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). The second optional argument, \code{active}, selects pairs for which the observed count is at least \code{active}.  If both \code{levels} and \code{active} are
      specified, then \code{active} is ignored.  (Thus, indices passed as \code{levels} should correspond to indices when \code{levels} = NULL and \code{active} = 0.) This term can only be
      used with undirected bipartite networks.

      \ThreeNineFour}

    \item{\code{concurrent(by=NULL, levels=NULL)}  (binary) (undirected) (categorical nodal attribute)
    }{\emph{Concurrent node count:}
      This term adds one network statistic to the model, equal to the number of
      nodes in the network with degree 2 or higher. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details);
      it functions just like the \code{by} argument of the \code{degree} term.
      This term can only be used with undirected
      networks.}

    \item{\code{concurrentties(by=NULL, levels=NULL)}  (binary) (undirected) (categorical nodal attribute)
    }{\emph{Concurrent tie count:}
      This term adds one network statistic to the model, equal to the number of
      ties incident on each actor beyond the first. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details);
      it functions just like the \code{by} argument of the \code{degree} term.
      This term can only be used with undirected
      networks.}

    \item{\code{ctriple(attr=NULL, diff=FALSE, levels=NULL)}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
      ,
      a.k.a. \code{ctriad}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    }{\emph{Cyclic triples:} By default, this term adds one
      statistic to the model, equal to the number of cyclic triples in the
      network, defined as a set of edges of the form \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k), (k{\rightarrow}i)\}}{\{(i,j), (j,k), (k,i)\}}. Note that
      for all directed networks, \code{triangle} is equal to
      \code{ttriple+ctriple}, so at most two of these three terms can be in a
      model. The optional argument \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  If \code{attr} is
      specified and \code{diff} is \code{FALSE}, then the statistic is the number of cyclic triples where all
      three nodes have the same value of the attribute.  If \code{attr} is specified and \code{diff} is \code{TRUE}, then one statistic is added to the model for each value of \code{attr} (or each
      value of \code{attr} specified by \code{levels} if that argument is passed), equal to the number of cyclic triples where all
      three nodes have that value of the attribute.  \Directed}

    \item{\code{cycle(k, semi=FALSE)} (binary) (directed) (undirected)}{\emph{k-Cycle Census:}
      The \code{k} argument must be a vector of integers giving the cycle lengths to count.
      Directed cycle lengths may range from \code{2} to \code{N} (the network size); undirected
      cycle lengths and semicycle lengths may range from \code{3} to \code{N}; length 2 semicycles
      are not currently supported.  Note that directed 2-cycles are equivalent to mutual 
      dyads.
      
      This term adds one network statistic to the model for each value of \code{k},
      corresponding to the number of \code{k}-cycles (or, alternately, semicycles)
      in the graph.
      
      The optional argument \code{semi} is a logical indicating whether semicycles 
      (rather than directed cycles) should be counted; this is ignored in the 
      undirected case.
      
      This term can be used with either directed or undirected networks.
    }

    \item{\code{cyclicalties(attr=NULL, levels=NULL)}   (binary)  (directed), \code{cyclicalties(threshold=0)}   (valued)  (directed)  (undirected)}{\emph{Cyclical ties:}
      This term adds one statistic, equal to the number of ties
      \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
      \eqn{j} to \eqn{i}.  (Related to the \code{ttriple} term.)
      The binary version takes a nodal attribute \code{attr}, and, if
      given, all three nodes
      involved (\eqn{i}, \eqn{j}, and the node on the two-path) must match
      on this attribute in order for \eqn{i\rightarrow j}{i-->j} to be counted.
      The binary version of this term can only be used with directed
      networks. The valued version can be used with both directed and undirected.}

    \item{\code{cyclicalweights(twopath="min",combine="max",affect="min")} (valued)  (directed)  (undirected)}{\emph{Cyclical
	weights:} This statistic implements the cyclical weights
      statistic, like that defined by Krivitsky (2012), Equation 13,
      but with the focus dyad being \eqn{y_{j,i}} rather than
      \eqn{y_{i,j}}. The
      currently implemented options for \code{twopath} is the minimum
      of the constituent dyads (\code{"min"}) or their geometric mean
      (\code{"geomean"}); for \code{combine}, the maximum of the
      2-path strengths (\code{"max"}) or their sum (\code{"sum"}); and
      for \code{affect}, the minimum of the focus dyad and the
      combined strength of the two paths (\code{"min"}) or their
      geometric mean (\code{"geomean"}). For each of these options,
      the first (and the default) is more stable but also more
      conservative, while the second is more sensitive but more likely
      to induce a multimodal distribution of networks.}

    \item{\code{ddsp(d, type="OTP")} (binary) (directed)}{  \emph{Directed dyadwise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{dyads} in the network with exactly \code{d[i]} shared partners. \Directed

      \SPTypes}

    \item{\code{degrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)} (binary)  (undirected) (categorical nodal attribute)
    }{\emph{Degree range:}
      The \code{from} and \code{to} arguments are vectors of distinct
      integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals the number of nodes in the network of degree
      greater than or equal to
      \code{from[i]} but strictly less than \code{to[i]}, i.e. with edges
      in semiopen interval \code{[from,to)}. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree range
      statistics are calculated for nodes having each separate
      value of the attribute.

      This term can only be used with undirected networks; for directed networks
      see \code{idegrange} and \code{odegrange}. This term can be used
      with bipartite networks, and will count nodes of both first and second mode in
      the specified degree range. To count only nodes of the first mode ("actors"), use \code{b1degrange}
      and to count only those fo the second mode ("events"), use \code{b2degrange}.}

    \item{\code{degree(d, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (undirected)  (categorical nodal attribute) (frequently-used)}{\emph{Degree:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the number of nodes in the network of degree
      \code{d[i]}, i.e. with exactly \code{d[i]} edges. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.
      This term can only be used with undirected networks; for directed networks
      see \code{idegree} and \code{odegree}.}

    \item{\code{degree1.5}  (binary)  (undirected)
    }{\emph{Degree to the 3/2 power:} This
      term adds one network statistic to the model equaling the sum over
      the actors of each actor's degree taken to the 3/2 power (or,
      equivalently, multiplied by its square root). This term is an
      undirected analog to the terms of Snijders et al. (2010), equations
      (11) and (12). This term can only be used with undirected networks.}

    \item{\code{degreepopularity}  (binary)  (undirected) (deprecated)
    }{\emph{Degree popularity (deprecated):} see \code{degree1.5}.}

    \item{\code{degcrossprod}  (binary)  (undirected)}{\emph{Degree Cross-Product:}
      This term adds one network statistic equal to the mean of the cross-products
      of the degrees of all pairs of nodes in the network which are tied.
      Only coded for undirected networks.}

    \item{\code{degcor}  (binary)  (undirected)}{\emph{Degree Correlation:}
      This term adds one network statistic equal to the correlation
      of the degrees of all pairs of nodes in the network which are tied.
      Only coded for undirected networks.}

    \item{\code{density}  (binary)  (dyad-independent)  (directed)  (undirected)}{\emph{Density:}
      This term adds one network statistic equal to the density of the network.
      For undirected networks, \code{density} equals \code{kstar(1)} or
      \code{edges} divided by \eqn{n(n-1)/2}; for directed networks,
      \code{density} equals \code{edges} or \code{istar(1)} or \code{ostar(1)}
      divided by \eqn{n(n-1)}.}

    \item{\code{diff(attr, pow=1, dir="t-h", sign.action="identity")} (binary) (dyad-independent) (frequently-used) (directed) (undirected) (quantitative nodal attribute), \code{diff(attr, pow=1, dir="t-h", sign.action="identity", form ="sum")} (valued) (dyad-independent) (directed) (undirected) (bipartite) (quantitative nodal attribute)}{\emph{Difference:}
      The \code{attr} argument specifies a quantitative vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). For values of \code{pow} other than
      \code{0}, this term adds one network statistic to the model,
      equaling the sum, over directed edges \eqn{(i,j)}, of
      \code{sign.action(attr[i]-attr[j])^pow} if \code{dir} is
      \code{"t-h"} (the default), \code{"tail-head"}, or \code{"b1-b2"}
      and of \code{sign.action(attr[j]-attr[i])^pow} if
      \code{"h-t"}, \code{"head-tail"}, or \code{"b2-b1"}. That is, the
      argument \code{dir} determines which vertex's attribute is
      subtracted from which, with tail being the origin of a directed edge
      and head being its destination, and bipartite networks' edges being
      treated as going from the first part (b1) to the second (b2).

      If \code{pow==0}, the exponentiation is replaced by the signum
      function: \code{+1} if the difference is positive, \code{0} if there
      is no difference, and \code{-1} if the difference is negative. Note
      that this function is applied \emph{after} the
      \code{sign.action}. The comparison is exact, so when using
      calculated values of \code{attr}, ensure that values that you
      want to be considered equal are, in fact, equal.

      The following \code{sign.actions} are possible:
      \describe{
	\item{\code{"identity"} (the default)}{no transformation of the
	  difference regardless of sign}
	\item{\code{"abs"}}{absolute value of the difference: equivalent
	  to the \code{absdiff} term}
	\item{\code{"posonly"}}{positive differences are kept, negative
	  differences are replaced by 0}
	\item{\code{"negonly"}}{negative differences are kept, positive
	  differences are replaced by 0}
      }

      Note that this term may not be meaningful for unipartite undirected
      networks unless \code{sign.action=="abs"}. When used on such a
      network, it behaves as if all edges were directed, going from the
      lower-indexed vertex to the higher-indexed vertex.}

    \item{\code{desp(d, type="OTP")} (binary) (directed)}{  \emph{Directed edgewise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{edges} in the network with exactly \code{d[i]} shared partners. \Directed

      \SPTypes}

    \item{\code{dgwdsp(decay=0, fixed=FALSE, cutoff=30, type="OTP")} (binary) (directed)}{  \emph{Geometrically weighted dyadwise shared partner distribution:} This term adds one network statistic to the model equal to the geometrically weighted dyadwise shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).  Note that the GWDSP statistic is equal to the sum of GWNSP plus GWESP.

      \SPTypes

      \GWCutoff{DSP}}

    \item{\code{dgwesp(decay=0, fixed=FALSE, cutoff=30, type="OTP")} (binary) (directed)}{ \emph{Geometrically weighted edgewise shared partner distribution:} This term adds a statistic equal to the geometrically weighted \emph{edgewise} (not dyadwise) shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).

      \SPTypes

      \GWCutoff{ESP}}

    \item{\code{dgwnsp(decay=0, fixed=FALSE, cutoff=30, type="OTP")} (binary) (directed)}{ \emph{Geometrically weighted non-edgewise shared partner distribution:} This term is just like gwesp and gwdsp except it adds a statistic equal to the geometrically weighted nonedgewise (that is, over dyads that do not have an edge) shared partner distribution with decay parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for this parameter may be fixed (if \code{fixed=TRUE}), or it may be used instead as the starting value for the estimation of \code{decay} in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).

      \SPTypes

      \GWCutoff{NSP}}

    \item{\code{dnsp(d, type="OTP")} (binary) (directed)}{  \emph{Directed non-edgewise shared partners:} This term adds one network statistic to the model for each element in \code{d} where the \eqn{i}th such statistic equals the number of \emph{non-edges} in the network with exactly \code{d[i]} shared partners. \Directed

      \SPTypes
    }

    \item{\code{dsp(d)}  (binary)  (directed)  (undirected)}{\emph{Dyadwise shared partners:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the number of dyads in the network with exactly
      \code{d[i]} shared partners.  This term can be used with directed and
      undirected networks.
      
      \SPtoDSP{dsp}{(directed) dyad \code{(i,j)}}{ddsp}}

    \item{\code{dyadcov(x, attrname=NULL)}  (binary)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Dyadic covariate:} The \code{x}
      argument is either a square matrix of covariates, one for each possible
      edge in the network, the name of a network attribute of covariates,
      or a network; if the latter, optional argument \code{attrname} provides
      the name of the quantitative edge attribute to use for covariate values
      (in this case, missing edges in \code{x} are assigned a covariate value
      of zero).
      This term adds three statistics to the model, each equal to the sum of the
      covariate values for all dyads occupying one of the three possible non-empty
      dyad states (mutual, upper-triangular asymmetric, and lower-triangular
      asymmetric dyads, respectively), with the empty or null state serving as a
      reference category. If the network is undirected, \code{x} is either a
      matrix of edgewise covariates, or a network; if the latter, optional
      argument \code{attrname} provides the name of the edge attribute to use for
      edge values. This term adds one statistic to the model, equal to the sum of
      the covariate values for each edge appearing in the network. The
      \code{edgecov} and \code{dyadcov} terms are equivalent for undirected
      networks.}

    \item{\code{edgecov(x, attrname=NULL)}  (binary)  (dyad-independent)  (directed)  (undirected) (frequently-used)
      , \code{edgecov(x,
	attrname=NULL, form="sum")}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Edge covariate:} The \code{x}
      argument is either a square matrix of covariates, one for each possible
      edge in the network, the name of a network attribute of covariates,
      or a network; if the latter, optional argument \code{attrname} provides
      the name of the quantitative edge attribute to use for covariate values
      (in this case, missing edges in \code{x} are assigned a covariate value
      of zero). This term adds one statistic to the model, equal to the sum
      of the covariate values for each edge appearing in the network. The
      \code{edgecov} term applies to both directed and undirected networks. For
      undirected networks the covariates are also assumed to be undirected. The
      \code{edgecov} and \code{dyadcov} terms are equivalent for undirected
      networks.}

    \item{\code{edges}  (binary)  (valued)  (dyad-independent)  (directed)  (undirected) (frequently-used)
      , a.k.a \code{nonzero}  (valued)  (directed)  (undirected) (dyad-independent)
    }{\emph{Edges:} This term adds one network statistic equal
      to the number of edges (i.e. nonzero values) in the network. For undirected networks, \code{edges}
      is equal to \code{kstar(1)}; for directed networks, \code{edges} is equal to
      both \code{ostar(1)} and \code{istar(1)}.}

    \item{\code{esp(d)}  (binary)  (directed)  (undirected)}{\emph{Edgewise shared partners:}
      This is just like the \code{dsp} term, except this term adds one network
      statistic to the model for each element in \code{d} where the \eqn{i}th such
      statistic equals the number of \emph{edges} (rather than dyads) in the
      network with exactly \code{d[i]} shared partners. This term can be used with
      directed and undirected networks.
      
      \SPtoDSP{esp}{(directed) edge \code{i -> j}}{desp}}

    \item{\code{equalto(value=0, tolerance=0)}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Number of dyads with
	values equal to a specific value (within tolerance):}
      Adds one statistic equal to the number of dyads whose values
      are within \code{tolerance} of \code{value}, i.e., between
      \code{value-tolerance} and \code{value+tolerance}, inclusive.}

    \item{\code{greaterthan(threshold=0)}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Number of dyads with
	values strictly greater than a threshold:}
      Adds one statistic equal to the number of dyads whose values
      exceed \code{threshold}.}

    \item{\code{gwb1degree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
	degree distribution for the first mode in a bipartite (aka two-mode)
	network:}
      This term adds one network statistic to the model equal to the weighted
      degree distribution with decay controlled by the \code{decay} parameter, which should be non-negative,
      for nodes in the
      first mode of a bipartite network. The first mode of a bipartite network
      object is sometimes known as the "actor" mode.
      The \code{decay} parameter is the same as theta_s in
      equation (14) in Hunter (2007).  The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used as merely the starting value for the estimation
      in a curved exponential family model (the default).

      \GWCutoff{degree}

      If \code{attr} is specified (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details) then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.
      This term can only be used with undirected bipartite
      networks.}

    \item{\code{gwb1dsp(decay=0, fixed=FALSE, cutoff=30)}  (binary)  (bipartite) (undirected) (curved)}{\emph{Geometrically weighted
	dyadwise shared partner distribution for dyads in the first bipartition:}
      This term adds one network statistic to the model equal to the geometrically
      weighted dyadwise shared partner distribution for dyads in the first bipartition,
      with decay parameter \code{decay} parameter, which should be non-negative. The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).
      This term can only be used with bipartite networks.

      \GWCutoff{b1dsp}}

    \item{\code{gwb2degree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
	degree distribution for the second mode in a bipartite (aka two-mode)
	network:}
      This term adds one network statistic to the model equal to the weighted
      degree distribution with decay controlled by the which should be non-negative,
      for nodes in the
      second mode of a bipartite network. The second mode of a bipartite network
      object is sometimes known as the "event" mode.
      The \code{decay} parameter is the same as theta_s in
      equation (14) in Hunter (2007).  The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used as merely the starting value for the estimation
      in a curved exponential family model (the default).

      \GWCutoff{degree}

      If \code{attr} is specified (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details) then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.
      This term can only be used with undirected bipartite
      networks.}

    \item{\code{gwb2dsp(decay=0, fixed=FALSE, cutoff=30)}  (binary)  (bipartite) (undirected) (curved)}{\emph{Geometrically weighted
	dyadwise shared partner distribution for dyads in the second bipartition:}
      This term adds one network statistic to the model equal to the geometrically
      weighted dyadwise shared partner distribution for dyads in the second bipartition,
      with decay parameter \code{decay} parameter, which should be non-negative. The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).
      This term can only be used with bipartite networks.

      \GWCutoff{b2dsp}}

    \item{\code{gwdegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)}  (binary)  (undirected)  (curved) (frequently-used)
    }{\emph{Geometrically weighted
	degree distribution:}
      This term adds one network statistic to the model equal to the weighted
      degree distribution with decay controlled by the \code{decay} parameter.
      The \code{decay} parameter is the same as theta_s in
      equation (14) in Hunter (2007).  The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).

      \GWCutoff{degree}

      If \code{attr} is specified (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details) then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.  This term can only be used with undirected networks.}

    \item{\code{gwdsp(decay=0, fixed=FALSE, cutoff=30)}  (binary)  (directed)  (undirected)  (curved)}{\emph{Geometrically weighted
	dyadwise shared partner distribution:}
      This term adds one network statistic to the model equal to the geometrically
      weighted dyadwise shared partner distribution with decay parameter
      \code{decay} parameter, which should be non-negative. The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006).
      This term can be used with directed and
      undirected networks.
      
      \SPtoDSP{gwdsp}{(directed) dyad \code{(i,j)}}{dgwdsp}

      \GWCutoff{DSP}
    }

    \item{\code{gwesp(decay=0, fixed=FALSE, cutoff=30)}  (binary)  (frequently-used)   (directed)  (undirected)  (curved)}{\emph{Geometrically weighted
	edgewise shared partner distribution:}
      This term is just like \code{gwdsp} except it adds a statistic equal to the
      geometrically weighted \emph{edgewise} (not dyadwise) shared partner
      distribution with decay parameter
      \code{decay} parameter, which should be non-negative. The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This term can be used with directed and
      undirected networks.
      
      \SPtoDSP{gwesp}{(directed) edge \code{i -> j}}{dgwesp}

      \GWCutoff{ESP}
    }

    \item{\code{gwidegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted
	in-degree distribution:} This term adds one network statistic to the model
      equal to the weighted in-degree distribution with decay parameter
      \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This
      term can only be used with directed networks.

      \GWCutoff{degree}

      If \code{attr} is specified (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details) then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.}

    \item{\code{gwnsp(decay=0, fixed=FALSE, cutoff=30)}  (binary)  (directed)  (undirected) (curved)}{\emph{Geometrically weighted
	nonedgewise shared partner distribution:} This term is just like
      \code{gwesp} and \code{gwdsp} except it adds a statistic equal to
      the geometrically weighted \emph{nonedgewise} (that is, over dyads
      that do not have an edge) shared partner distribution with weight
      parameter \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The optional argument \code{fixed} indicates
      whether the \code{decay} parameter is fixed at the given value, or is to be fit as a curved
      exponential-family model (see Hunter and Handcock, 2006).  The
      default is \code{FALSE}, which means the scale parameter is not
      fixed and thus the model is a CEF model. This term can be used with
      directed and undirected networks.
      
      \SPtoDSP{gwnsp}{(directed) non-edge \code{(i,j)}}{dgwnsp}


      \GWCutoff{NSP}
    }

    \item{\code{gwodegree(decay, fixed=FALSE, attr=NULL, cutoff=30, levels=NULL)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted
	out-degree distribution:} This term adds one network statistic to the model
      equal to the weighted out-degree distribution with decay parameter
      \code{decay} parameter, which should be non-negative. (this parameter was called \code{alpha} prior to \code{ergm 3.7}). The value supplied for
      this parameter may be fixed (if \code{fixed=TRUE}),
      or it may be used instead as the starting value for the estimation of \code{decay}
      in a curved exponential family model (when \code{fixed=FALSE}, the default) (see Hunter and Handcock, 2006). This
      term can only be used with directed networks.

      \GWCutoff{degree}

      If \code{attr} is specified (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details) then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.}

    \item{\code{hamming(x, cov, attrname=NULL)}  (binary)  (dyad-independent)  (directed)  (undirected)}{\emph{Hamming distance:}
      This term adds one statistic to the model equal to the weighted or
      unweighted Hamming distance of the network from the network specified by
      \code{x}. (If no argument is given, \code{x} is taken to be the observed
      network, i.e., the network on the left side of the \eqn{\sim}{~} in the formula
      that defines the ERGM.) Unweighted Hamming distance is defined as the total
      number of pairs \eqn{(i,j)} (ordered or unordered, depending on whether the
      network is directed or undirected) on which the two networks differ. If the
      optional argument \code{cov} is specified, then the weighted Hamming
      distance is computed instead, where each pair \eqn{(i,j)} contributes a
      pre-specified weight toward the distance when the two networks differ on
      that pair. The argument \code{cov} is either a matrix of edgewise weights or
      a network; if the latter, the optional argument \code{attrname} provides the
      name of the edge attribute to use for weight values.}

    \item{\code{hammingmix(attr, x, base=NULL, levels=NULL, levels2=NULL)}  (binary)  (directed)  (dyad-independent)}{\emph{Hamming
	distance within mixing:}
      By default, this term adds one statistic to the model for every possible pairing of
      attribute values of the network for the vertex attribute specified by \code{attr} (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels}).  Each such statistic is the Hamming distance (i.e., the number of differences)
      between the appropriate subset of dyads in
      the network and the corresponding subset in \code{x}.  The ordering of the
      attribute values is alphabetical by default, but this may be overridden using the \code{levels} arguments.
      The optional arguments \code{levels} and \code{levels2} allow the user to control what pairings are included (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels}), and the order in which they appear. For example \code{levels2=-2} will
      omit the second statistic, making it the de facto reference
      category.

      \BaseDepII

      \Directed}

    \item{\code{idegrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (directed) (categorical nodal attribute)
    }{\emph{In-degree range:}
      The \code{from} and \code{to} arguments are vectors of distinct
      integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals the number of nodes in the network of in-degree
      greater than or equal to
      \code{from[i]} but strictly less than \code{to[i]}, i.e. with
      in-edge count
      in semiopen interval \code{[from,to)}. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree range
      statistics are calculated for nodes having each separate
      value of the attribute.

      This term can only be used with directed networks; for undirected
      networks (bipartite and not)
      see \code{degrange}. For degrees of specific modes of bipartite
      networks, see \code{b1degrange} and \code{b2degrange}. For
      in-degrees, see \code{idegrange}.}

    \item{\code{idegree(d, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{In-degree:} The \code{d} argument
      is a vector of distinct integers. This term adds one network statistic to
      the model for each element in \code{d}; the \eqn{i}th such statistic equals
      the number of nodes in the network of in-degree \code{d[i]}, i.e. the number
      of nodes with exactly \code{d[i]} in-edges. The optional term
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.
      This term can only be used with directed networks; for undirected networks
      see \code{degree}.}

    \item{\code{idegree1.5}  (binary)  (directed)
    }{\emph{In-degree to the 3/2 power:}
      This term adds one network statistic to the model equaling the sum
      over the actors of each actor's indegree taken to the 3/2 power
      (or, equivalently, multiplied by its square root). This term is
      analogous to the term of Snijders et al. (2010), equation (12). This
      term can only be used with directed networks.}

    \item{\code{idegreepopularity}  (binary)  (directed) (deprecated)
    }{\emph{In-degree popularity (deprecated):} see \code{idegree1.5}.}

    \item{\code{ininterval(lower=-Inf, upper=+Inf, open=c(TRUE,TRUE))}  (valued)  (directed)  (undirected) (dyad-independent)
    }{\emph{Number of dyads
	whose values are in an interval}
      Adds one statistic equaling to the number of dyads whose values
      are between \code{lower} and
      \code{upper}. Argument \code{open}
      is a \code{logical} vector of length 2 that controls whether
      the interval is open (exclusive) on the lower and on the upper
      end, respectively. \code{open} can also be specified as one of
      \code{"[]"}, \code{"(]"}, \code{"[)"}, and \code{"()"}.}

    \item{\code{intransitive}  (binary)  (directed)  (triad-related)}{\emph{Intransitive triads:}
      This term adds one statistic to the model, equal to the number of triads in
      the network that are intransitive. The intransitive triads are those of type
      \code{111D}, \code{201}, \code{111U}, \code{021C}, or \code{030C} in the
      categorization of Davis and Leinhardt (1972). For details on the 16 possible
      triad types, see \code{\link[sna]{triad.classify}} in the
      \code{\link[sna]{sna}} package. Note the distinction from the \code{ctriple}
      term. \Directed}

    \item{\code{isolatededges}  (binary)  (undirected) (bipartite)
    }{\emph{Isolated edges:} This term adds one statistic to the
      model equal to the number of isolated edges in the network, i.e., the number
      of edges each of whose endpoints has degree 1.  This term can only be used
      with undirected networks.}      
      
    \item{\code{isolates}  (binary)  (directed)  (undirected) (frequently-used)
    }{\emph{Isolates:} This term adds one statistic to the
      model equal to the number of isolates in the network. For an undirected
      network, an isolate is defined to be any node with degree zero. For a
      directed network, an isolate is any node with both in-degree and out-degree
      equal to zero.}

    \item{\code{istar(k, attr=NULL, levels=NULL)}  (binary)  (directed) (categorical nodal attribute)
    }{\emph{In-stars:} The \code{k} argument is a
      vector of distinct integers. This term adds one network statistic to the
      model for each element in \code{k}. The \eqn{i}th such statistic counts the
      number of distinct \code{k[i]}-instars in the network, where a
      \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
      different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the ties
      \eqn{(O_j{\rightarrow}N)}{(O_j, N)} exist for \eqn{j=1, \dots, k}. The
      optional argument \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified
      then the count is over the number of \eqn{k}-instars where all nodes have
      the same value of the attribute. This term can only be used for directed
      networks; for undirected networks see \code{kstar}. Note that
      \code{istar(1)} is equal to both \code{ostar(1)} and \code{edges}.}

    \item{\code{kstar(k, attr=NULL, levels=NULL)}  (binary)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars:}
      The \code{k} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{k}. The \eqn{i}th
      such statistic counts the number of distinct \code{k[i]}-stars in the
      network, where a \eqn{k}-star is defined to be a node \eqn{N} and a set of
      \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the ties
      \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}. The optional argument
      \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified then the count is over
      the number of \eqn{k}-stars where all nodes have the same value of the
      attribute. This term can only be used for undirected networks; for directed
      networks, see \code{istar}, \code{ostar}, \code{twopath} and \code{m2star}.
      Note that \code{kstar(1)} is equal to \code{edges}.}

    \item{\code{smallerthan(threshold=0)}  (valued)  (directed)  (undirected) (dyad-independent)}{\emph{Number of dyads with
	values strictly smaller than a threshold:}
      Adds one statistic equaling to the number of dyads whose values
      exceeded by \code{threshold}.}

    \item{\code{localtriangle(x)}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Triangles within neighborhoods:}
      This term adds one statistic to the model equal to the number of triangles
      in the network between nodes \dQuote{close to} each other. For an undirected
      network, a local triangle is defined to be any set of three edges between
      nodal pairs \eqn{\{(i,j), (j,k), (k,i)\}} that are in the same neighborhood.
      For a directed network, a triangle is defined as any set of three edges
      \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)} and either
      \eqn{(k{\rightarrow}i)} or \eqn{(k{\leftarrow}i)} where again all nodes are
      within the same neighborhood. The argument \code{x} is an undirected
      network or an symmetric adjacency matrix that specifies whether the two nodes
      are in the same neighborhood. Note that \code{triangle}, with or without an argument, is a
      special case of \code{localtriangle}.}

    \item{\code{m2star}  (binary)  (directed)}{\emph{Mixed 2-stars, a.k.a 2-paths:}
      This term adds one statistic to the model, equal to the number of mixed
      2-stars in the network, where a mixed 2-star is a pair of distinct edges
      \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)}. A mixed 2-star is
      sometimes called a 2-path because it is a directed path of length 2 from
      \eqn{i} to \eqn{k} via \eqn{j}. However, in the case of a 2-path the focus
      is usually on the endpoints \eqn{i} and \eqn{k}, whereas for a mixed 2-star
      the focus is usually on the midpoint \eqn{j}. This term can only be used
      with directed networks; for undirected networks see \code{kstar(2)}. See
      also \code{twopath}.}

    \item{\code{meandeg}  (binary)  (dyad-independent)  (directed)  (undirected)}{\emph{Mean vertex degree:}
      This term adds one network statistic to the model equal to the
      average degree of a node.  Note that this term is a constant multiple of
      both \code{edges} and \code{density}.}

    \item{\code{mm(attrs, levels=NULL, levels2=NULL)} (binary)  (dyad-independent)  (frequently-used) (directed)  (undirected) (categorical nodal attribute), \code{mm(attrs, levels=NULL, levels2=NULL, form="sum")} (valued)  (dyad-independent)  (frequently-used) (directed)  (undirected) (categorical nodal attribute)}{\emph{Mixing matrix cells and margins:}
      \code{attrs} is a two-sided formula whose LHS gives the attribute or
      attribute function (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels}) for the rows of the mixing matrix and whose RHS gives
      that for its columns. A one-sided formula (e.g.,
      \code{~A}) is symmetrized (e.g., \code{A~A}). \code{levels}
      similarly specifies the subset of rows and columns to be
      used. \code{levels2} can then be used to filter which specific cells
      of the matrix to include. A two-sided formula with a dot on one side
      calculates the margins of the mixing matrix, analogously to \code{nodefactor}, with
      \code{A~.} calculating the row/sender/b1 margins and \code{.~A}
      calculating the column/receiver/b2 margins.}

    \item{\code{mutual(same=NULL, by=NULL, diff=FALSE, keep=NULL, levels=NULL)} (binary)  (directed)  (frequently-used), \code{mutual(form="min",threshold=0)}  (valued)  (directed)
    }{\emph{Mutuality:}
      In binary ERGMs, equal to the number of
      pairs of actors \eqn{i} and \eqn{j} for which \eqn{(i{\rightarrow}j)}{(i,j)}
      and \eqn{(j{\rightarrow}i)}{(j,i)} both exist. For valued ERGMs, equal to \eqn{\sum_{i<j} m(y_{i,j},y_{j,i})},
      where \eqn{m} is determined by \code{form} argument: \code{"min"}
      for \eqn{\min(y_{i,j},y_{j,i})}, \code{"nabsdiff"} for
      \eqn{-|y_{i,j},y_{j,i}|}, \code{"product"} for
      \eqn{y_{i,j}y_{j,i}}, and \code{"geometric"} for
      \eqn{\sqrt{y_{i,j}}\sqrt{y_{j,i}}}. See Krivitsky (2012) for a
      discussion of these statistics. \code{form="threshold"} simply
      computes the binary \code{mutuality} after
      thresholding at \code{threshold}.

      This term can only be used
      with directed networks. The binary version also has the following capabilities: if the optional \code{same} argument is passed
      (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details),
      only mutual pairs that match on the attribute are counted;
      separate counts for each unique matching value can be obtained by using
      \code{diff=TRUE} with \code{same}; and if \code{by} is passed (again, see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels}),
      then each node is counted separately for each mutual pair in which it
      occurs and the counts are tabulated by unique values of the attribute.
      This means that the sum of the mutual statistics when \code{by} is used
      will equal twice the standard mutual statistic.  Only one of \code{same}
      or \code{by} may be used, and only the former is affected by \code{diff};
      if both \code{same} and \code{by} are passed, \code{by} is ignored.
      Finally, if \code{levels} is passed, this tells which statistics should be kept
      whenever the \code{mutual} term would ordinarily result in multiple statistics (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels}).

      \KeepDep}

    \item{\code{nearsimmelian}  (binary)  (directed)  (triad-related)}{\emph{Near simmelian triads:}
      This term adds one statistic to the model equal to the number of near
      Simmelian triads, as defined by Krackhardt and Handcock (2007). This is a
      sub-graph of size three which is exactly one tie short of being complete.
      \Directed}

    \item{\code{nodecov(attr)}  (binary)  (dyad-independent)  (frequently-used)   (directed)  (undirected) (quantitative nodal attribute), \code{nodecov(attr, form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (quantitative nodal attribute), a.k.a. \code{nodemain}  (binary)   (directed)  (undirected)}{\emph{Main effect of a covariate:}
      The \code{attr} argument specifies one or more quantitative attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the sum of
      \code{attr(i)} and \code{attr(j)} for all edges \eqn{(i,j)} in the
      network. For categorical attributes, see \code{nodefactor}. Note that for
      directed networks, \code{nodecov} equals \code{nodeicov} plus
      \code{nodeocov}.

      \ThreeNineFour}

    \item{\code{nodecovar}  (valued)  (directed)  (undirected) (quantitative nodal attribute)
    }{\emph{Uncentered covariance of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} (y_{i,j}y_{i,k}+y_{k,j}y_{k,j})}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic.}

    \item{\code{nodefactor(attr, base=1, levels=-1)}  (binary)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute) (frequently-used)
      ,
      \code{nodefactor(attr, base=1, levels=-1, form="sum")}  (dyad-independent)  (valued)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Factor attribute
	effect:} The \code{attr} argument specifies one or more categorical attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds multiple network statistics to the
      model, one for each of (a subset of) the unique values of the
      \code{attr} attribute (or each combination of the attributes
      given). Each of these statistics gives the number of times a node
      with that attribute or those attributes appears in an edge in the
      network.

      \LevelsExplain{nodefactor}

      \LevelsNotFirst

      \BaseDep}

    \item{\code{nodeicov(attr)}  (binary)  (directed) (quantitative nodal attribute) (frequently-used)
      , \code{nodeicov(attr, form="sum")}  (valued)  (directed) (quantitative nodal attribute)
    }{\emph{Main effect of a covariate for
	in-edges:}
      The \code{attr} argument specifies one or more quantitative attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
      value of \code{attr(j)} for all edges \eqn{(i,j)} in the network. This
      term may only be used with directed networks. For categorical attributes,
      see \code{nodeifactor}.

      \ThreeNineFour}

    \item{\code{nodeicovar}  (valued)  (directed) (quantitative nodal attribute)
    }{\emph{Uncentered covariance of in-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{k,j}y_{k,j}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.}

    \item{\code{nodeifactor(attr, base=1, levels=-1)}  (binary)  (dyad-independent)  (directed) (categorical nodal attribute) (frequently-used)
      ,
      \code{nodeifactor(attr, base=1, levels=-1, form="sum")}  (valued)  (dyad-independent)  (directed) (categorical nodal attribute)
    }{\emph{Factor attribute
	effect for in-edges:} The \code{attr} argument specifies one or more categorical attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds multiple network
      statistics to the model, one for each of (a subset of) the unique
      values of the \code{attr} attribute (or each combination of the
      attributes given). Each of these statistics gives the number of
      times a node with that attribute or those attributes appears as the
      terminal node of a directed tie.

      \LevelsExplain{nodeifactor}

      \LevelsNotFirst

      \BaseDep

      For an analogous term for quantitative vertex attributes, see \code{nodeicov}.}

    \item{\code{nodeisqrtcovar}  (valued)  (directed)  (non-negative)
    }{\emph{Uncentered covariance of square roots of in-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{k,j}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.}

    \item{\code{nodematch(attr, diff=FALSE, keep=NULL, levels=NULL)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (categorical nodal attribute)
      , \code{nodematch(attr, diff=FALSE, keep=NULL, levels=NULL, form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
      a.k.a. \code{match}  (binary)  (directed)  (dyad-independent)  (undirected) (categorical nodal attribute)
    }{\emph{Uniform homophily and differential homophily:}
      The \code{attr} argument specifies one or more attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  When \code{diff=FALSE}, this term adds one network statistic
      to the model, which counts the number of edges \eqn{(i,j)} for which
      \code{attr(i)==attr(j)}. This is also called ''uniform homophily,'' because each group is assumed to have the same propensity for within-group ties.  When multiple attribute names are given, the
      statistic counts only ties for which all of the attributes
      match. When \code{diff=TRUE}, \eqn{p} network statistics are added
      to the model, where \eqn{p} is the number of unique values of the
      \code{attr} attribute. The \eqn{k}th such statistic counts the
      number of edges \eqn{(i,j)} for which \code{attr(i) ==
	attr(j) == value(k)}, where \code{value(k)} is the \eqn{k}th
      smallest unique value of the \code{attr} attribute. This is also called ''differential homophily,'' because each group is allowed to have a unique propensity for within-group ties.  Note that a statistical test of uniform vs. differential homophily should be conducted using the ANOVA function.

      By default, matches on all levels \eqn{k} are
      counted. \LevelsExplain{nodematch} This works for both
      \code{diff=TRUE} and \code{diff=FALSE}.

      \KeepDep}

    \item{\code{nodemix(attr, base=NULL, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=NULL)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (categorical nodal attribute)
      ,
      \code{nodemix(attr, base=NULL, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=NULL, form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Nodal attribute
	mixing:} The \code{attr} argument specifies one or more categorical vertex attributes
      (see \link[=nodal_attributes]{Specifying Vertex Attributes and Levels} for details). By default, this term adds one network statistic to
      the model for each possible pairing of attribute values. The
      statistic equals the number of edges in the network in which the
      nodes have that pairing of values. (When multiple attributes are specified, a
      statistic is added for each combination of attribute values for
      those attributes.) In other words, this term produces one statistic for
      every entry in the mixing matrix for the attribute(s). By default, the ordering of
      the attribute values is lexicographic: alphabetical (for nominal categories) or
      numerical (for ordered categories), but this can be overridden using
      the \code{levels} arguments.
      The optional arguments \code{levels}, \code{levels2},
      \code{b1levels}, and \code{b2levels} control what statistics are included in the model, and the order in which they appear.
      \code{levels2} apply to all networks; \code{levels} applies to unipartite networks; \code{b1levels} and \code{b2levels} apply to
      bipartite networks (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels}).

      \BaseDepII}

    \item{\code{nodeocov(attr)}  (binary)  (directed) (dyad-independent)(quantitative nodal attribute)
      , \code{nodeocov(attr, form="sum")}  (valued)  (directed) (dyad-independent)
      (quantitative nodal attribute)
    }{\emph{Main effect of a covariate for
	out-edges:}
      The \code{attr} argument specifies one or more quantitative attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  This term adds a single network statistic for each quantitative attribute or matrix column to the model equaling the total
      value of \code{attr(i)} for all edges \eqn{(i,j)} in the network. This
      term may only be used with directed networks. For categorical attributes,
      see \code{nodeofactor}.

      \ThreeNineFour}

    \item{\code{nodeocovar}  (valued)  (directed) (quantitative nodal attribute)
    }{\emph{Uncentered covariance of out-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{i,j}y_{i,k}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.}

    \item{\code{nodeofactor(attr, base=1, levels=-1)}  (binary)  (dyad-independent)  (directed) (categorical nodal attribute)
      ,
      \code{nodeofactor(attr, base=1, levels=-1, form="sum")}  (valued)  (dyad-independent)  (categorical nodal attribute) (directed)}{\emph{Factor attribute
	effect for out-edges:} The \code{attr} argument specifies one or more categorical attributes (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds multiple network
      statistics to the model, one for each of (a subset of) the unique
      values of the \code{attr} attribute (or each combination of the
      attributes given). Each of these statistics gives the number of
      times a node with that attribute or those attributes appears as the
      node of origin of a directed tie.

      \LevelsExplain{nodeofactor}

      \LevelsNotFirst

      \BaseDep

      \Directed}

    \item{\code{nodeosqrtcovar}  (valued)  (directed)  (non-negative)
    }{\emph{Uncentered covariance of square roots of out-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{i,k}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.}

    \item{\code{nodesqrtcovar(center=TRUE)}  (valued)  (non-negative)  (directed)  (undirected)
    }{\emph{Covariance of square roots of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k}
	(\sqrt{y_{i,j}}\sqrt{y_{i,k}}+\sqrt{y_{k,j}}\sqrt{y_{k,j}})} if
      \code{center=FALSE}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic. If \code{center=TRUE} (the default), the statistic is
      instead \eqn{\sum_{i,j,k}
	((\sqrt{y_{i,j}}-\bar{\sqrt{y}})(\sqrt{y_{i,k}}-\bar{\sqrt{y}})+(\sqrt{y_{k,j}}-\bar{\sqrt{y}})(\sqrt{y_{k,j}}-\bar{\sqrt{y}}))},
      where \eqn{\bar{\sqrt{y}}} is the mean of the square root of
      dyad values.}

    \item{\code{nsp(d)}  (binary)  (directed)  (undirected)}{\emph{Nonedgewise shared partners:} This is
      just like the \code{dsp} and \code{esp} terms, except this term adds
      one network statistic to the model for each element in \code{d}
      where the \eqn{i}th such statistic equals the number of
      \emph{non-edges} (that is, dyads that do not have an edge) in the network
      with exactly \code{d[i]} shared partners. This term can be used with
      directed and undirected networks. 
      
      \SPtoDSP{nsp}{(directed) non-edge \code{(i,j)}}{dnsp}}

    \item{\code{odegrange(from, to=+Inf, by=NULL, homophily=FALSE, levels=NULL)} (binary)  (directed) (categorical nodal attribute)
    }{\emph{Out-degree range:}
      The \code{from} and \code{to} arguments are vectors of distinct
      integers (or \code{+Inf}, for \code{to} (its default)). If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals the number of nodes in the network of out-degree
      greater than or equal to
      \code{from[i]} but strictly less than \code{to[i]}, i.e. with
      out-edge count
      in semiopen interval \code{[from,to)}. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree range
      statistics are calculated for nodes having each separate
      value of the attribute.

      This term can only be used with directed networks; for undirected
      networks (bipartite and not)
      see \code{degrange}. For degrees of specific modes of bipartite
      networks, see \code{b1degrange} and \code{b2degrange}. For
      in-degrees, see \code{idegrange}.}

    \item{\code{odegree(d, by=NULL, homophily=FALSE, levels=NULL)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{Out-degree:} The \code{d} argument
      is a vector of distinct integers. This term adds one network statistic to
      the model for each element in \code{d}; the \eqn{i}th such statistic equals
      the number of nodes in the network of out-degree \code{d[i]}, i.e. the
      number of nodes with exactly \code{d[i]} out-edges. The optional argument
      \code{by} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).
      If this is specified and \code{homophily} is \code{TRUE},
      then degrees are calculated using the subnetwork consisting of only
      edges whose endpoints have the same value of the \code{by} attribute.
      If \code{by} is specified and
      \code{homophily} is \code{FALSE} (the default), then separate degree
      statistics are calculated for nodes having each separate
      value of the attribute.
      This term can only be used with directed networks; for undirected networks
      see \code{degree}.}

    \item{\code{odegree1.5}  (binary)  (directed)
    }{\emph{Out-degree to the 3/2 power:}
      This term adds one network statistic to the model equaling the sum
      over the actors of each actor's outdegree taken to the 3/2 power
      (or, equivalently, multiplied by its square root). This term is
      analogous to the term of Snijders et al. (2010), equation (12). This
      term can only be used with directed networks.}

    \item{\code{odegreepopularity}  (binary)  (directed) (deprecated)
    }{\emph{Out-degree popularity (deprecated):} see \code{odegree1.5}.}

    \item{\code{opentriad}  (binary)  (undirected) (triad-related)}{\emph{Open triads:} This term
      adds one statistic to the model equal to the number of 2-stars minus
      three times the number of triangles in the network. It is currently
      only implemented for undirected networks.}

    \item{\code{ostar(k, attr=NULL, levels=NULL)}  (binary)  (directed) (categorical nodal attribute)
    }{\emph{k-Outstars:} The \code{k} argument is
      a vector of distinct integers. This term adds one network statistic to the
      model for each element in \code{k}. The \eqn{i}th such statistic counts the
      number of distinct \code{k[i]}-outstars in the network, where a
      \eqn{k}-outstar is defined to be a node \eqn{N} and a set of \eqn{k}
      different nodes \eqn{\{O_1, \dots, O_k\}}{\{O[1], ..., O[k]\}} such that the ties
      \eqn{(N{\rightarrow}O_j)}{(N,O_j)} exist for \eqn{j=1, \dots, k}. The
      optional argument \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If this is specified
      then the count is the number of \eqn{k}-outstars where all nodes have the
      same value of the attribute. This term can only be used with directed
      networks; for undirected networks see \code{kstar}. Note that
      \code{ostar(1)} is equal to both \code{istar(1)} and \code{edges}.}

    \item{\code{receiver(base=1, nodes=-1)}  (binary)  (directed)  (dyad-independent)
      ,
      \code{receiver(base=1, nodes=-1, form="sum")}  (valued)  (directed)  (dyad-independent)}{\emph{Receiver effect:}
      This term adds one network statistic for each node equal to the number of
      in-ties for that node. This measures the popularity of the node. The term
      for the first node is omitted by default because of linear dependence that
      arises if this term is used together with \code{edges}, but its coefficient
      can be computed as the negative of the sum of the coefficients of all the
      other actors. That is, the average coefficient is zero, following the
      Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
      1981). The \code{base} and \code{nodes} arguments allow the user to determine which nodes'
      statistics should be included or excluded (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  The argument \code{nodes} is preferred to \code{base}, although \code{base}
      carries a default value of 1 for backwards compatibility.  (If both \code{base} and \code{nodes} are supplied, then \code{nodes} overrides \code{base}.)  This
      term can only be used with directed networks. For undirected networks, see
      \code{sociality}.}

    \item{\code{sender(base=1, nodes=-1)}  (binary)  (directed)  (dyad-independent)
      ,
      \code{sender(base=1, nodes=-1, form="sum")}  (valued)  (directed)  (dyad-independent)}{\emph{Sender effect:}
      This term adds one network statistic for each node equal to the number of
      out-ties for that node. This measures the activity of the node. The term for
      the first node is omitted by default because of linear dependence that
      arises if this term is used together with \code{edges}, but its coefficient
      can be computed as the negative of the sum of the coefficients of all the
      other actors. That is, the average coefficient is zero, following the
      Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
      1981). The \code{nodes} arguments allow the user to determine which nodes'
      statistics should be included or excluded (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).

      \BaseDepNodes

      \Directed For undirected networks, see \code{sociality}.}

    \item{\code{simmelian}  (binary)  (directed)  (triad-related)}{\emph{Simmelian triads:} This term adds one
      statistic to the model equal to the number of Simmelian triads, as defined
      by Krackhardt and Handcock (2007). This is a complete sub-graph of size
      three. \Directed}

    \item{\code{simmelianties}  (binary)  (triad-related)  (directed)}{\emph{Ties in simmelian triads:} This term adds
      one statistic to the model equal to the number of ties in the network that
      are associated with Simmelian triads, as defined by Krackhardt and Handcock
      (2007). Each Simmelian has six ties in it but, because Simmelians can
      overlap in terms of nodes (and associated ties), the total number of ties in
      these Simmelians is less than six times the number of Simmelians. Hence this
      is a measure of the clustering of Simmelians (given the number of
      Simmelians). \Directed}

    \item{\code{smalldiff(attr, cutoff)}  (binary)  (dyad-independent)  (directed)  (undirected) (quantitative nodal attribute)
    }{\emph{Number of
	ties between actors with similar (but not necessarily identical)
	attribute values:} The \code{attr} argument specifies a quantitative vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). This term adds one statistic, having as its
      value the number of edges in the network for which the incident
      actors' attribute values differ less than \code{cutoff}; that is,
      number of edges between \code{i} to \code{j} such that
      \code{abs(attr[i]-attr[j])<cutoff}.}

    \item{\code{sociality(attr=NULL, base=1, levels=NULL, nodes=-1)}  (binary) (undirected) (dyad-independent) (categorical nodal attribute)
      ,
      \code{sociality(attr=NULL, base=1, levels=NULL, nodes=-1, form="sum")}  (valued) (undirected) (dyad-independent) (categorical nodal attribute)
    }{\emph{Undirected degree:}
      This term adds one network statistic for each node equal to the number of
      ties of that node. The optional \code{attr} argument specifies a categorical vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details). If provided, this term only counts ties
      between nodes with the same value of the attribute (an actor-specific
      version of the \code{nodematch} term), restricted to be one of the values specified by \code{levels} if \code{levels} is not \code{NULL}.  \Undirected  For directed networks, see \code{sender} and
      \code{receiver}. By default, \code{nodes=-1} means that the statistic for the
      first node will be omitted, but this argument may be changed to control
      which statistics are included just as for the \code{nodes} argument of \code{sender} and
      \code{receiver} terms.

      \BaseDepNodes}

    \item{\code{sum(pow=1)}  (valued)  (directed)  (undirected)}{\emph{Sum of dyad values (optionally taken
	to a power):} This term adds one statistic equal to the sum of
      dyad values taken to the power \code{pow}, which defaults to 1.}

    \item{\code{threetrail(keep=NULL, levels=NULL)} (binary)  (directed)  (undirected)
      (triad-related),
    }{\emph{ Three-trails:}
      a.k.a. \code{threepath}.
      For an undirected network, this term adds one statistic equal to the number
      of 3-trails, where a 3-trail is defined as a \dQuote{trail} of length three that
      traverses three distinct edges.
      Note that a 3-trail need not
      include four distinct nodes; in particular, a triangle counts as three
      3-trails.  For a directed network, this term adds four statistics
      (or some subset of these four specified by the \code{levels} argument),
      one for each of the four distinct types of directed three-paths. If the
      nodes of the path are written from left to right such that the middle edge
      points to the right (R), then the four types are RRR, RRL, LRR, and LRL.
      That is, an RRR 3-trail is of the form
      \eqn{i\rightarrow j\rightarrow k\rightarrow l}{i-->j-->k-->l}, and RRL
      3-trail is of the form
      \eqn{i\rightarrow j\rightarrow k\leftarrow l}{i-->j-->k<--l}, etc.
      Like in the undirected case, there is no requirement that the nodes be
      distinct in a directed 3-trail.  However, the three edges must all be
      distinct.  Thus, a mutual tie \eqn{i\leftrightarrow j}{i<-->j} does not
      count as a 3-trail of the form
      \eqn{i\rightarrow j\rightarrow i\leftarrow j}{i-->j-->i<--j}; however,
      in the subnetwork \eqn{i\leftrightarrow j \rightarrow k}{i<-->j-->k},
      there are two directed 3-trails, one LRR
      (\eqn{k\leftarrow j\rightarrow i\leftarrow j}{k<--j-->i-->j})
      and one RRR
      (\eqn{j\rightarrow i\rightarrow j\leftarrow k}{k<--j-->i-->j}).

      \KeepDep This term used to be (inaccurately) called \code{threepath}. That
      name has been deprecated and may be removed in a future version.}

    \item{\code{transitive}  (binary)  (directed)  (triad-related)}{\emph{Transitive triads:}
      This term adds one statistic to the model, equal to the number of triads in
      the network that are transitive. The transitive triads are those of type
      \code{120D}, \code{030T}, \code{120U}, or \code{300} in the categorization
      of Davis and Leinhardt (1972). For details on the 16 possible triad types,
      see \code{\link[sna]{triad.classify}} in the \code{\link[sna]{sna}} package.
      Note the distinction from the \code{ttriple} term. This term can only be
      used with directed networks.}

    \item{\code{transitiveties(attr=NULL, levels=NULL)}   (binary)  (directed) (triad-related) (categorical nodal attribute)
      , \code{transitiveties(threshold=0)}   (valued)  (directed)  (undirected) (triad-related)}{\emph{Transitive ties:}
      This term adds one statistic, equal to the number of ties
      \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
      \eqn{i} to \eqn{j}.  (Related to the \code{ttriple} term.)
      The binary version takes a nodal attribute \code{attr}, and, if
      given, all three nodes
      involved (\eqn{i}, \eqn{j}, and the node on the two-path) must match
      on this attribute in order for \eqn{i\rightarrow j}{i-->j} to be counted.
      The binary version of this term can only be used with directed
      networks. The valued version can be used with both directed and
      undirected.}

    \item{\code{transitiveweights(twopath="min",combine="max",affect="min")}  (valued)  (directed)  (undirected)  (non-negative) (triad-related)
    }{\emph{Transitive
	weights:} This statistic implements the transitive weights
      statistic defined by Krivitsky (2012), Equation 13. The
      currently implemented options for \code{twopath} is the minimum
      of the constituent dyads (\code{"min"}) or their geometric mean
      (\code{"geomean"}); for \code{combine}, the maximum of the
      2-path strengths (\code{"max"}) or their sum (\code{"sum"}); and
      for \code{affect}, the minimum of the focus dyad and the
      combined strength of the two paths (\code{"min"}) or their
      geometric mean (\code{"geomean"}). For each of these options,
      the first (and the default) is more stable but also more
      conservative, while the second is more sensitive but more likely
      to induce a multimodal distribution of networks.}

    \item{\code{triadcensus(levels)}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Triad census:}
      For a directed network, this term adds one network statistic for each of
      an arbitrary subset of the 16 possible types of triads categorized by
      Davis and Leinhardt (1972) as \code{003, 012, 102, 021D, 021U, 021C, 111D,
	111U, 030T, 030C, 201, 120D, 120U, 120C, 210,} and \code{300}. Note that at
      least one category should be dropped; otherwise a linear dependency will
      exist among the 16 statistics, since they must sum to the total number of
      three-node sets. By default, the category \code{003}, which is the category
      of completely empty three-node sets, is dropped. This is considered category
      zero, and the others are numbered 1 through 15 in the order given above. By
      using the \code{levels} argument (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details), the user may specify a set of terms to add other than the default
      value of \code{1:15}. Each statistic is the count of the corresponding triad
      type in the network. For details on the 16 types, see \code{?triad.classify}
      in the \code{{sna}} package, on which this code is based. For an undirected
      network, the triad census is over the four types defined by the number of
      ties (i.e., 0, 1, 2, and 3), and the default is to add \code{1:3}, which is
      to say that the 0 is dropped; however, this too may be controlled by
      changing the \code{levels} argument.}

    \item{\code{triangle(attr=NULL, diff=FALSE, levels=NULL)}  (binary)  (frequently-used)  (triad-related)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Triangles:}
      By default, this term adds one statistic to the model equal to the number of triangles
      in the network. For an undirected network, a triangle is defined to be any
      set \eqn{\{(i,j), (j,k), (k,i)\}} of three edges. For a directed network, a
      triangle is defined as any set of three edges \eqn{(i{\rightarrow}j)}{(i,j)}
      and \eqn{(j{\rightarrow}k)}{(j,k)} and either \eqn{(k{\rightarrow}i)}{(k,i)}
      or \eqn{(k{\leftarrow}i)}{(i,k)}. The former case is called a ``transitive
      triple'' and the latter is called a ``cyclic triple'', so in the case of a
      directed network, \code{triangle} equals \code{ttriple} plus \code{ctriple}
      --- thus at most two of these three terms can be in a model. The optional
      argument \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  If \code{attr} is specified and \code{diff} is \code{FALSE},	
      then the count is restricted to those triples of nodes with
      equal values of the vertex attribute specified by \code{attr}. If \code{attr} is specified and \code{diff} is \code{TRUE},
      then one statistic is added for each value of \code{attr} (or each value specified by \code{levels} if that argument is passed),
      equal to the number of triangles where all three nodes have that value of the attribute.}

    \item{\code{tripercent(attr=NULL, diff=FALSE, levels=NULL)}  (binary)  (undirected)  (triad-related) (categorical nodal attribute)
    }{\emph{Triangle percentage:}
      By default, this term adds one statistic to the model equal to 100 times the ratio of
      the number of triangles in the network to the sum of the number of triangles
      and the number of 2-stars not in triangles (the latter is considered a
      potential but incomplete triangle).  In case the denominator equals zero,
      the statistic is defined to be zero. For the definition of triangle, see
      \code{triangle}. The optional argument \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  If \code{attr} is specified and \code{diff} is \code{FALSE}, the counts
      (both numerator and denominator) are restricted to those triples of nodes with equal values
      of the vertex attribute specified by \code{attr}.  If \code{attr} is specified and \code{diff} is \code{TRUE}, then one statistic is added
      for each value of \code{attr} (or each value specified by \code{levels} if that argument is passed), where the counts
      (both numerator and denominator) are restricted to those triples of nodes with that value
      of the vertex attribute specified by \code{attr}.  This is often called
      the mean correlation coefficient. This term can only be
      used with undirected networks; for directed networks, it is difficult to
      define the numerator and denominator in a consistent and meaningful way.}

    \item{\code{ttriple(attr=NULL, diff=FALSE, levels=NULL)}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
      , a.k.a. \code{ttriad}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    }{\emph{Transitive triples:}
      By default, this term adds one statistic to the model, equal to the number of transitive
      triples in the network, defined as a set of edges \eqn{\{(i{\rightarrow}j),
	(j{\rightarrow}k), (i{\rightarrow}k)\}}{\{(i,j), (j,k), (i,k)\}}. Note that
      \code{triangle} equals \code{ttriple+ctriple} for a directed network, so at
      most two of the three terms can be in a model. The optional argument
      \code{attr} specifies a vertex attribute (see \link[=nodal_attributes]{Specifying Vertex
	Attributes and Levels} for details).  If \code{attr} is specified and \code{diff} is \code{FALSE}, then the count is over
      the number of transitive triples where all three nodes have the same value
      of the attribute.  If \code{attr} is specified and \code{diff} is \code{TRUE}, then one statistic is added for each value of \code{attr}
      (or each value of \code{attr} specified by \code{levels} if that argument is passed), equal to the number of transitive triples where all three nodes
      have that value of \code{attr}.  \Directed}

    \item{\code{twopath}  (binary)  (directed)  (undirected)}{\emph{2-Paths:}
      This term adds one statistic to the model, equal to the number of 2-paths in
      the network. For a directed network this is defined as a pair of edges
      \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)}, where \eqn{i} and
      \eqn{j} must be distinct. That is, it is a directed path of length 2 from
      \eqn{i} to \eqn{k} via \eqn{j}. For directed networks a 2-path is also a
      mixed 2-star but the interpretation is usually different; see \code{m2star}.
      For undirected networks a twopath is defined as a pair of edges
      \eqn{\{i,j\}, \{j,k\}}. That is, it is an undirected path of length 2 from
      \eqn{i} to \eqn{k} via \eqn{j}, also known as a 2-star.}
  }
}
\references{
  \itemize{
    \item Bomiriya, R. P, Bansal, S., and Hunter, D. R. (2014).  Modeling
    Homophily in ERGMs for Bipartite Networks.  Submitted.

    \item Butts, CT.  (2008).  \dQuote{A Relational Event Framework for Social Action.} \emph{Sociological Methodology,} 38(1).

    \item Davis, J.A. and Leinhardt, S.  (1972).  The Structure of Positive Interpersonal
    Relations in Small Groups.  In J. Berger (Ed.), \emph{Sociological Theories in Progress,
      Volume 2}, 218--251.  Boston: Houghton Mifflin.

    \item Holland, P. W. and S. Leinhardt (1981). An exponential family
    of probability distributions for directed graphs.
    \emph{Journal of the American Statistical Association}, 76:
    33--50.

    \item Hunter, D. R. and M. S. Handcock (2006). Inference in curved
    exponential family models for networks. \emph{Journal of Computational
      and Graphical Statistics}, 15: 565--583.

    \item Hunter, D. R. (2007). Curved exponential family models for social
    networks. \emph{Social Networks}, 29: 216--230.

    \item Krackhardt, D. and Handcock, M. S. (2007).  Heider versus Simmel: Emergent Features
    in Dynamic Structures. \emph{Lecture Notes in Computer Science}, 4503, 14--27.

    \item Krivitsky P. N. (2012). Exponential-Family Random Graph Models for Valued
    Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
    1100-1128. \doi{10.1214/12-EJS696}

    \item Robins, G; Pattison, P; and Wang, P.  (2009).  \dQuote{Closure, Connectivity, and Degree Distributions: Exponential Random Graph (p*) Models for Directed Social Networks.}  \emph{Social Networks,} 31:105-117.

    \item Snijders T. A. B., G. G. van de Bunt, and
    C. E. G. Steglich. Introduction to Stochastic Actor-Based Models for
    Network Dynamics. \emph{Social Networks}, 2010, 32(1),
    44-60. \doi{10.1016/j.socnet.2009.02.004}

    \item Morris M, Handcock MS, and Hunter DR. Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects. \emph{Journal of Statistical Software}, 2008, 24(4), 1-24. \url{https://www.jstatsoft.org/v24/i04}

    \item Snijders, T. A. B., P. E. Pattison, G. L. Robins, and
    M. S. Handcock (2006). New specifications for exponential random graph
    models, \emph{Sociological Methodology}, 36(1): 99-153.
  }
}
\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{search.ergmTerms}}, \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\examples{
\dontrun{
ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)

ergm(molecule ~ edges + kstar(2:3) + triangle
                      + nodematch("atomic type",diff=TRUE)
                      + triangle + absdiff("atomic type"))
}
% TODO: Write a valued example.
}
\keyword{models}
