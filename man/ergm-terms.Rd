%  File man/ergm-terms.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2015 Statnet Commons
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
\alias{b1concurrent}
\alias{b1degrange}
\alias{b1degree}
\alias{b1factor}
\alias{b1mindegree}
\alias{b1nodematch}
\alias{b1star}
\alias{b1starmix}
\alias{b1twostar}
\alias{b2concurrent}
\alias{b2degrange}
\alias{b2degree}
\alias{b2factor}
\alias{b2mindegree}
\alias{b2nodematch}
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
\alias{degcor}
\alias{degcrossprod}
\alias{degrange}
\alias{degree}
\alias{degreepopularity}
\alias{density}
\alias{dsp}
\alias{dyadcov}
\alias{edgecov}
\alias{edges}
\alias{esp}
\alias{greaterthan}
\alias{gwb1degree}
\alias{gwb2degree}
\alias{gwdegree}
\alias{gwdsp}
\alias{gwesp}
\alias{gwidegree}
\alias{gwnsp}
\alias{gwodegree}
\alias{hamming}
\alias{hammingmix}
\alias{idegrange}
\alias{idegree}
\alias{idegreepopularity}
\alias{ininterval}
\alias{intransitive}
\alias{isolates}
\alias{istar}
\alias{kstar}
\alias{localtriangle}
\alias{m2star}
\alias{match}
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
\alias{odegreepopularity}
\alias{opentriad}
\alias{ostar}
\alias{receiver}
\alias{sender}
\alias{simmelian}
\alias{simmelianties}
\alias{smalldiff}
\alias{sociality}
\alias{sum}
\alias{threepath}
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
\alias{InitErgmTerm.absdiff}
\alias{InitErgmTerm.absdiffcat}
\alias{InitErgmTerm.altkstar}
\alias{InitErgmTerm.asymmetric}
\alias{InitErgmTerm.b1concurrent}
\alias{InitErgmTerm.b1degrange}
\alias{InitErgmTerm.b1degree}
\alias{InitErgmTerm.b1factor}
\alias{InitErgmTerm.b1mindegree}
\alias{InitErgmTerm.b1nodematch}
\alias{InitErgmTerm.b1star}
\alias{InitErgmTerm.b1starmix}
\alias{InitErgmTerm.b1twostar}
\alias{InitErgmTerm.b2concurrent}
\alias{InitErgmTerm.b2degrange}
\alias{InitErgmTerm.b2degree}
\alias{InitErgmTerm.b2factor}
\alias{InitErgmTerm.b2nodematch}
\alias{InitErgmTerm.b2mindegree}
\alias{InitErgmTerm.b2star}
\alias{InitErgmTerm.b2starmix}
\alias{InitErgmTerm.b2twostar}
\alias{InitErgmTerm.balance}
\alias{InitErgmTerm.coincidence}
\alias{InitErgmTerm.concurrent}
\alias{InitErgmTerm.concurrentties}
\alias{InitErgmTerm.ctriple}
\alias{InitErgmTerm.ctriad}
\alias{InitErgmTerm.cycle}
\alias{InitErgmTerm.cyclicalties}
\alias{InitErgmTerm.degcor}
\alias{InitErgmTerm.degcrossprod}
\alias{InitErgmTerm.degrange}
\alias{InitErgmTerm.degree}
\alias{InitErgmTerm.degreepopularity}
\alias{InitErgmTerm.density}
\alias{InitErgmTerm.dsp}
\alias{InitErgmTerm.dyadcov}
\alias{InitErgmTerm.edgecov}
\alias{InitErgmTerm.edges}
\alias{InitErgmTerm.esp}
\alias{InitErgmTerm.gwb1degree}
\alias{InitErgmTerm.gwb2degree}
\alias{InitErgmTerm.gwdegree}
\alias{InitErgmTerm.gwdsp}
\alias{InitErgmTerm.gwesp}
\alias{InitErgmTerm.gwidegree}
\alias{InitErgmTerm.gwnsp}
\alias{InitErgmTerm.gwodegree}
\alias{InitErgmTerm.hamming}
\alias{InitErgmTerm.hammingmix}
\alias{InitErgmTerm.idegrange}
\alias{InitErgmTerm.idegree}
\alias{InitErgmTerm.idegreepopularity}
\alias{InitErgmTerm.intransitive}
\alias{InitErgmTerm.isolates}
\alias{InitErgmTerm.istar}
\alias{InitErgmTerm.kstar}
\alias{InitErgmTerm.localtriangle}
\alias{InitErgmTerm.m2star}
\alias{InitErgmTerm.match}
\alias{InitErgmTerm.meandeg}
\alias{InitErgmTerm.mutual}
\alias{InitErgmTerm.nearsimmelian}
\alias{InitErgmTerm.nodecov}
\alias{InitErgmTerm.nodefactor}
\alias{InitErgmTerm.nodeicov}
\alias{InitErgmTerm.nodeifactor}
\alias{InitErgmTerm.nodemain}
\alias{InitErgmTerm.nodematch}
\alias{InitErgmTerm.nodemix}
\alias{InitErgmTerm.nodeocov}
\alias{InitErgmTerm.nodeofactor}
\alias{InitErgmTerm.nsp}
\alias{InitErgmTerm.odegrange}
\alias{InitErgmTerm.odegree}
\alias{InitErgmTerm.odegreepopularity}
\alias{InitErgmTerm.opentriad}
\alias{InitErgmTerm.ostar}
\alias{InitErgmTerm.receiver}
\alias{InitErgmTerm.sender}
\alias{InitErgmTerm.simmelian}
\alias{InitErgmTerm.simmelianties}
\alias{InitErgmTerm.smalldiff}
\alias{InitErgmTerm.sociality}
\alias{InitErgmTerm.threepath}
\alias{InitErgmTerm.transitive}
\alias{InitErgmTerm.transitiveties}
\alias{InitErgmTerm.triadcensus}
\alias{InitErgmTerm.triangle}
\alias{InitErgmTerm.triangles}
\alias{InitErgmTerm.tripercent}
\alias{InitErgmTerm.ttriple}
\alias{InitErgmTerm.ttriad}
\alias{InitErgmTerm.twopath}
\alias{InitWtErgmTerm.absdiff}
\alias{InitWtErgmTerm.absdiffcat}
\alias{InitWtErgmTerm.atleast}
\alias{InitWtErgmTerm.cyclicalties}
\alias{InitWtErgmTerm.cyclicalweights}
\alias{InitWtErgmTerm.edgecov}
\alias{InitWtErgmTerm.edges}
\alias{InitWtErgmTerm.greaterthan}
\alias{InitWtErgmTerm.ininterval}
\alias{InitWtErgmTerm.match}
\alias{InitWtErgmTerm.mutual}
\alias{InitWtErgmTerm.nodecov}
\alias{InitWtErgmTerm.nodecovar}
\alias{InitWtErgmTerm.nodefactor}
\alias{InitWtErgmTerm.nodematch}
\alias{InitWtErgmTerm.nodemix}
\alias{InitWtErgmTerm.nodeicov}
\alias{InitWtErgmTerm.nodeicovar}
\alias{InitWtErgmTerm.nodeifactor}
\alias{InitWtErgmTerm.nodeisqrtcovar}
\alias{InitWtErgmTerm.nodemain}
\alias{InitWtErgmTerm.nodeocov}
\alias{InitWtErgmTerm.nodeocovar}
\alias{InitWtErgmTerm.nodeofactor}
\alias{InitWtErgmTerm.nodeosqrtcovar}
\alias{InitWtErgmTerm.nodesqrtcovar}
\alias{InitWtErgmTerm.nonzero}
\alias{InitWtErgmTerm.sum}
\alias{InitWtErgmTerm.transitiveties}
\alias{InitWtErgmTerm.transitiveweights}
\alias{InitErgmTerm.conddegmple}
\alias{InitErgmTerm.indices}
\docType{methods}
\title{Terms used in Exponential Family Random Graph Models}
\description{
  The function \code{\link{ergm}} is used to fit exponential 
  random graph models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}, where
  \eqn{h(y)} is the reference measure (for valued network models),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same 
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  The network statistics \eqn{g(y)} are entered as terms in the 
  function call to \code{\link{ergm}}.

  This page describes the possible terms (and hence network statistics)
  included in \code{\link[=ergm-package]{ergm}} package. Other packages
  may add their own terms, and package
  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} provides tools
  for implementing them.

  The current recommendation for any package implementing additional
  terms is to create a help file with a name or alias \code{ergm-terms},
  so that \code{help("ergm-terms")} will list ERGM terms available from
  all loaded packages.
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
  then add nodal attributes to it using the \code{\%v\%} operator if necessary.}

\section{Binary and valued ERGM terms}{

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

\section{Terms to represent network statistics included in the
  \code{\link[=ergm-package]{ergm}} package}{
  A cross-referenced html version of the term documentation is is available via \code{vignette('ergm-term-crossRef')} and terms can also be searched via \code{\link{search.ergmTerms}}.
  \describe{
    %%%%%%%%%%%%%%%%%%%%%%555
    % This section below is specially formatted so that it can be parsed to build additional 
    % documentation about ergm terms. (parsing functions are in build_term_index.R)
    % The list of term tags must appear on the same line as the term name to be parsed correctly
    %%beginTerms%%
    \item{\code{absdiff(attrname, pow=1)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (quantitative nodal attribute),
    \code{absdiff(attrname, pow=1, form ="sum")}  (valued)  (dyad-independent) (directed)  (undirected) (quantitative nodal attribute)}{\emph{Absolute difference:}
    The \code{attrname} argument is a character string giving the name
    of a quantitative attribute in the network's vertex attribute
    list. This term adds one network statistic to the model equaling the
    sum of \code{abs(attrname[i]-attrname[j])^pow} for all edges (i,j)
    in the network. }
    
    \item{\code{absdiffcat(attrname, base=NULL)}  (binary) (dyad-independent) (directed) (undirected) (categorical nodal attribute),
    \code{absdiffcat(attrname, base=NULL, form="sum")}  (valued)  (dyad-independent)   (directed)  (undirected) (categorical nodal attribute)}{\emph{Categorical absolute
    difference:} 
    The \code{attrname} argument is a character string giving the name of a
    quantitative attribute in the network's vertex attribute list. This term
    adds one statistic for every possible nonzero distinct value of
    \code{abs(attrname[i]-attrname[j])} in the network; the value of each such
    statistic is the number of edges in the network with the corresponding
    absolute difference. The optional \code{base} argument is a vector
    indicating which nonzero differences, in order from smallest to largest,
    should be omitted from the model (i.e., treated like the zero-difference
    category). The \code{base} argument, if used, should contain indices, not
    differences themselves. For instance, if the possible values of
    \code{abs(attrname[i]-attrname[j])} are 0, 0.5, 3, 3.5, and 10, then to omit
    0.5 and 10 one should set \code{base=c(1, 4)}. Note that this term should
    generally be used only when the quantitative attribute has a limited number
    of possible values; an example is the \code{"Grade"} attribute of the
    \code{\link{faux.mesa.high}} or \code{\link{faux.magnolia.high}} datasets.}
    
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
    \code{fixed} indicates whether the scale parameter \code{lambda} is to be
    fit as a curved exponential family model (see Hunter and Handcock, 2006).
    The default is \code{FALSE}, which means the scale parameter is not fixed
    and thus the model is a CEF model. This term can only be used with
    undirected networks.}
    
    \item{\code{asymmetric(attrname=NULL, diff=FALSE, keep=NULL)}  (binary)  (directed)  (dyad-independent) (triad-related)}{\emph{Asymmetric
    dyads:} This term adds one network statistic to the model equal to the
    number of pairs of actors for which exactly one of
    \eqn{(i{\rightarrow}j)}{(i,j)} or \eqn{(j{\rightarrow}i)}{(j,i)} exists.
    This term can only be used with directed networks. If the optional
    \code{attrname} argument is used, only asymmetric pairs that match on the
    named vertex attribute are counted. The optional modifiers \code{diff} and
    \code{keep} are used in the same way as for the \code{nodematch} term; refer
    to this term for details and an example.}

  \item{\code{atleast(threshold=0)}  (valued)  (directed)  (undirected)}{\emph{Number of ties with values
      greater than or equal to a threshold} Adds one statistic
    equaling to the number of ties whose values equal or exceed
    \code{threshold}.  }
  
  \item{\code{b1concurrent(by=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
  }{\emph{Concurrent node count for the
    first mode in a bipartite (aka two-mode) network:} This term adds one
    network statistic to the model, equal to the number of nodes in the first
    mode of the network with degree 2 or higher. The first mode of a bipartite
    network object is sometimes known as the "actor" mode. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list;
    it functions just like the \code{by} argument of the \code{b1degree} term.
    Without the optional argument, this statistic is equivalent to \code{b1mindegree(2)}.
    This term can only be
    used with undirected bipartite networks. }

  \item{\code{b1degrange(from, to=+Inf, by=NULL, homophily=FALSE)}  (binary)  (bipartite)  (undirected) 
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
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

    \item{\code{b1degree(d, by=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute) (frequently-used)
    }{\emph{Degree for the first mode in a
    bipartite (aka two-mode) network:} The \code{d} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{d}; the \eqn{i}th such statistic equals the number of
    nodes of degree \code{d[i]} in the first mode of a bipartite network, i.e.
    with exactly \code{d[i]} edges. The first mode of a bipartite network object
    is sometimes known as the "actor" mode. The optional argument \code{by} is
    a character string giving the name of an attribute in the network's vertex
    attribute list. If this is specified 
    then each node's degree is tabulated only with other nodes having the same
    value of the \code{by} attribute.
    This term can only be used with undirected bipartite networks.}
    
    \item{\code{b1factor(attrname, base=1)}  (binary)  (bipartite)  (undirected)  (dyad-independent) (frequently-used) (categorical nodal attribute)
    }{\emph{Factor attribute effect for
    the first mode in a bipartite (aka two-mode) network :} 
    The \code{attrname} argument is a character string giving the name of a
    categorical attribute in the network's vertex attribute list. This term adds
    multiple network statistics to the model, one for each of (a subset of) the
    unique values of the \code{attrname} attribute. Each of these statistics
    gives the number of times a node with that attribute in the first mode of
    the network appears in an edge. The first mode of a bipartite network object
    is sometimes known as the "actor" mode. To include all attribute values is
    usually not a good idea, because the sum of all such statistics equals the
    number of edges and hence a linear dependency would arise in any model also
    including \code{edges}. Thus, the \code{base} argument tells which value(s)
    (numbered in order according to the \code{sort} function) should be omitted.
    The default value, \code{base=1}, means that the smallest (i.e., first in
    sorted order) attribute value is omitted. For example, if the \dQuote{fruit}
    factor has levels \dQuote{orange}, \dQuote{apple}, \dQuote{banana}, and
    \dQuote{pear}, then to add just two terms, one for \dQuote{apple} and one
    for \dQuote{pear}, then set \dQuote{banana} and \dQuote{orange} to the base
    (remember to sort the values first) by using \code{nodefactor("fruit",
    base=2:3)}. This term can only be used with undirected bipartite networks.}
    
    \item{\code{b1mindegree(d)}  (binary)  (bipartite)  (undirected)}{\emph{Minimum degree for the first mode in a
    bipartite (aka two-mode) network:} The \code{d} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{d}; the \eqn{i}th such statistic equals the number of
    nodes in the first mode of a bipartite network with at least degree \code{d[i]}. 
    The first mode of a bipartite network object is sometimes known as the "actor" mode.
    This term can only be used with undirected bipartite networks.}
    
    \item{\code{b1nodematch(attrname, diff=FALSE, keep=NULL, by=NULL, alpha=1, beta=1, 
    byb2attr=NULL)}  (binary)  (bipartite)  (undirected) (dyadic-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Nodal attribute-based homophily effect for
    the first mode in a bipartite (aka two-mode) network:}  This term is introduced
    in Bomiriya et al (2014).
    The \code{attrname} argument is a character string giving the name of a
    categorical attribute in the network's vertex attribute list. Out of the two arguments
    (discount parameters) \code{alpha} and \code{beta}, both which takes values from [0,1], 
    only one should be set at a time. If none is set to a value other than 1, this term will 
    simply be a homophily based two-star statistic. This term adds one statistic to the model 
    unless \code{diff} is set to \code{TRUE}, in which case the term adds multiple network 
    statistics to the model, one for each of (a subset of) the unique values of the \code{attrname} 
    attribute. To include only the attribute values you wish, use the \code{keep} argument. 
    If an \code{alpha} discount parameter is used, each of these statistics gives the sum of 
    the number of common second-mode nodes raised to the power \code{alpha} for each pair of 
    first-mode nodes with that attribute. If a \code{beta} discount parameter is used, each 
    of these statistics gives half the sum of the number of two-paths with two first-mode nodes 
    with that attribute as the two ends of the two path raised to the power \code{beta} for each 
    edge in the network. The \code{byb2attr} argument is a character string giving the name of a
   	second mode categorical attribute in the network's attribute list. Setting this argument
   	will separate the orginal statistics based on the values of the set second mode attribute---
   	i.e. for example, if \code{diff} is \code{FALSE}, then the sum of all the statistics for 
   	each level of this second-mode attribute will be equal to the original \code{b1nodematch}
   	statistic where \code{byb2attr} set to \code{NULL}. This term can only be used with 
   	undirected bipartite networks.}
    
    \item{\code{b1star(k, attrname=NULL)}  (binary) (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars for the first mode in a
    bipartite (aka two-mode) network: } The \code{k} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{k}. The \eqn{i}th such statistic counts the number of
    distinct \code{k[i]}-stars whose center node is in the first mode of the
    network. The first mode of a bipartite network object is sometimes known as
    the "actor" mode. A \eqn{k}-star is defined to be a center node \eqn{N} and
    a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}} such that the
    ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}. The optional argument
    \code{attrname} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified then the count is over
    the number of \eqn{k}-stars (with center node in the first mode) where all
    nodes have the same value of the attribute. This term can only be used for
    undirected bipartite networks. Note that \code{b1star(1)} is equal to
    \code{b2star(1)} and to \code{edges}. }
    
    \item{\code{b1starmix(k, attrname, base=NULL, diff=TRUE)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Mixing
    matrix for k-stars centered on the first mode of a bipartite network:} Only
    a single value of \eqn{k} is allowed. This term counts all k-stars in which
    the b2 nodes (called events in some contexts) are homophilous in the sense
    that they all share the same value of \code{attrname}. However, the b1 node
    (in some contexts, the actor) at the center of the k-star does NOT have to
    have the same value as the b2 nodes; indeed, the values taken by the b1
    nodes may be completely distinct from those of the b2 nodes, which allows
    for the use of this term in cases where there are two separate nodal
    attributes, one for the b1 nodes and another for the b2 nodes (in this case,
    however, these two attributes should be combined to form a single nodal
    attribute called \code{attrname}. A different statistic is created for each
    value of \code{attrname} seen in a b1 node, even if no k-stars are observed
    with this value. Whether a different statistic is created for each value
    seen in a b2 node depends on the value of the \code{diff} argument: When
    \code{diff=TRUE}, the default, a different statistic is created for each
    value and thus the behavior of this term is reminiscent of the
    \code{nodemix} term, from which it takes its name; when \code{diff=FALSE},
    all homophilous k-stars are counted together, though these k-stars are still
    categorized according to the value of the central b1 node. The \code{base}
    term may be used to control which of the possible terms are left out of the
    model: By default, all terms are included, but if \code{base} is set to a
    vector of indices then the corresponding terms (in the order they would be
    created when \code{base=NULL}) are left out.}

    \item{\code{b1twostar(b1attrname, b2attrname, base=NULL)}  (binary) (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Two-star
    census for central nodes centered on the first mode of a bipartite network:}
    This term takes two nodal attribute names, one for b1 nodes (actors in some
    contexts) and one for b2 nodes (events in some contexts). Only
    \code{b1attrname} is required; if \code{b2attrname} is not passed, it is
    assumed to be the same as \code{b1attrname}. Assuming that there are
    \eqn{n_1} values of \code{b1attrname} among the b1 nodes and \eqn{n_2}
    values of \code{b2attrname} among the b2 nodes, then the total number of
    distinct categories of two stars according to these two attributes is
    \eqn{n_1(n_2)(n_2+1)/2}. This model term creates a distinct statistic
    counting each of these categories. The \code{base} term may be used to leave
    some of these categories out; when passed as a vector of integer indices (in
    the order the statistics would be created when \code{base=NULL}), the
    corresponding terms will be left out.}
    
    \item{\code{b2concurrent(by=NULL)}  (binary)  (bipartite)  (undirected) (frequently-used)
    }{\emph{Concurrent node count for the
    second mode in a bipartite (aka two-mode) network:} This term adds one
    network statistic to the model, equal to the number of nodes in the second
    mode of the network with degree 2 or higher. The second mode of a bipartite
    network object is sometimes known as the "event" mode. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list;
    it functions just like the \code{by} argument of the \code{b2degree} term.
    Without the optional argument, this statistic is equivalent to \code{b2mindegree(2)}.
    This term can only be used with undirected bipartite networks. }

    \item{\code{b2degrange(from, to=+Inf, by=NULL, homophily=FALSE)}  (binary)  (bipartite)  (undirected)
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified 
    then each node's degree is tabulated only with other nodes having the same
    value of the \code{by} attribute.
    This term can only be used with undirected bipartite networks.}
    
    \item{\code{b2factor(attrname, base=1)}  (binary)  (bipartite)  (undirected)  (dyad-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Factor attribute effect for
    the second mode in a bipartite (aka two-mode) network :} 
    The \code{attrname} argument is a character string giving the name of a
    categorical attribute in the network's vertex attribute list. This term adds
    multiple network statistics to the model, one for each of (a subset of) the
    unique values of the \code{attrname} attribute. Each of these statistics
    gives the number of times a node with that attribute in the second mode of
    the network appears in an edge. The second mode of a bipartite network
    object is sometimes known as the "event" mode. To include all attribute
    values is usually not a good idea, because the sum of all such statistics
    equals the number of edges and hence a linear dependency would arise in any
    model also including \code{edges}. Thus, the \code{base} argument tells
    which value(s) (numbered in order according to the \code{sort} function)
    should be omitted. The default value, \code{base=1}, means that the smallest
    (i.e., first in sorted order) attribute value is omitted. For example, if
    the \dQuote{fruit} factor has levels \dQuote{orange}, \dQuote{apple},
    \dQuote{banana}, and \dQuote{pear}, then to add just two terms, one for
    \dQuote{apple} and one for \dQuote{pear}, then set \dQuote{banana} and
    \dQuote{orange} to the base (remember to sort the values first) by using
    \code{nodefactor("fruit", base=2:3)}. This term can only be used with
    undirected bipartite networks.}
    
    \item{\code{b2mindegree(d)}  (binary) (bipartite)  (undirected)
    }{\emph{Minimum degree for the second mode in a
    bipartite (aka two-mode) network:} The \code{d} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{d}; the \eqn{i}th such statistic equals the number of
    nodes  in the second mode of a bipartite network with at least degree \code{d[i]}. 
    The second mode of a bipartite network object is sometimes known as the "event" mode.
    This term can only be used with undirected bipartite networks.}
    
    \item{\code{b2nodematch(attrname, diff=FALSE, keep=NULL, by=NULL, alpha=1, beta=1, 
    byb1attr=NULL)}  (binary)  (bipartite)  (undirected) (dyadic-independent) (categorical nodal attribute) (frequently-used)
    }{\emph{Nodal attribute-based homophily effect for
    the second mode in a bipartite (aka two-mode) network:} This term is introduced
    in Bomiriya et al (2014).
    The \code{attrname} argument is a character string giving the name of a
    categorical attribute in the network's vertex attribute list. Out of the two arguments
    (discount parameters) \code{alpha} and \code{beta}, both which takes values from [0,1], 
    only one should be set at a time. If none is set to a value other than 1, this term will 
    simply be a homophily based two-star statistic. This term adds one statistic to the model 
    unless \code{diff} is set to \code{TRUE}, in which case the term adds multiple network 
    statistics to the model, one for each of (a subset of) the unique values of the \code{attrname} 
    attribute. To include only the attribute values you wish, use the \code{keep} argument. 
    If an \code{alpha} discount parameter is used, each of these statistics gives the sum of 
    the number of common first-mode nodes raised to the power \code{alpha} for each pair of 
    second-mode nodes with that attribute. If a \code{beta} discount parameter is used, each 
    of these statistics gives half the sum of the number of two-paths with two second-mode nodes 
    with that attribute as the two ends of the two path raised to the power \code{beta} for each 
    edge in the network. The \code{byb1attr} argument is a character string giving the name of a
   	first mode categorical attribute in the network's attribute list. Setting this argument
   	will separate the orginal statistics based on the values of the set first mode attribute---
   	i.e. for example, if \code{diff} is \code{FALSE}, then the sum of all the statistics for 
   	each level of this first-mode attribute will be equal to the original \code{b2nodematch}
   	statistic where \code{byb1attr} set to \code{NULL}. This term can only be used with 
   	undirected bipartite networks.}

    \item{\code{b2star(k, attrname=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars for the second mode in a
    bipartite (aka two-mode) network: } The \code{k} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{k}. The \eqn{i}th such statistic counts the number of
    distinct \code{k[i]}-stars whose center node is in the second mode of the
    network. The second mode of a bipartite network object is sometimes known as
    the "event" mode. A \eqn{k}-star is defined to be a center node \eqn{N} and
    a set of \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}} such that the
    ties \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}. The optional argument
    \code{attrname} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified then the count is over
    the number of \eqn{k}-stars (with center node in the second mode) where all
    nodes have the same value of the attribute. This term can only be used for
    undirected bipartite networks. Note that \code{b2star(1)} is equal to
    \code{b1star(1)} and to \code{edges}. }
    
    \item{\code{b2starmix(k, attrname, base=NULL, diff=TRUE)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Mixing
    matrix for k-stars centered on the second mode of a bipartite network:} 
    This term is exactly the same as \code{b1starmix} except that the roles of
    b1 and b2 are reversed.}
    
    \item{\code{b2twostar(b1attrname, b2attrname, base=NULL)}  (binary)  (bipartite)  (undirected) (categorical nodal attribute)
    }{\emph{ Two-star
    census for central nodes centered on the second mode of a bipartite
    network:} This term is exactly the same as \code{b1twostar} except that the
    roles of b1 and b2 are reversed.}    

    \item{\code{balance}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Balanced triads:} 
    This term adds one network statistic to the model equal to the number of
    triads in the network that are balanced. The balanced triads are those of
    type \code{102} or \code{300} in the categorization of Davis and Leinhardt (1972). For details on the 16 possible triad types, see
    \code{?triad.classify} in the \code{{sna}} package. For an undirected
    network, the balanced triads are those with an even number of ties (i.e., 0
    and 2).}
    
    \item{\code{coincidence(d=NULL,active=0)}  (binary)  (bipartite)  (undirected) 
    }{\emph{Coincident node count for the
    second mode in a bipartite (aka two-mode) network:} By default this term adds one
    network statistic to the model for each pair of nodes of mode two. It is
equal to the number of (first mode) mutual partners of that pair.
    The first mode of a bipartite
    network object is sometimes known as the "actor" mode and the seconds as the "event" mode. So this is the number of actors going to both events in the pair. The optional argument
    \code{d} is a two-column matrix of (row-wise) pairs indices where the first row is less than the second row. The second optional argument, \code{active}, selects pairs for which the observed count is at least \code{active}.
    This term can only be
    used with undirected bipartite networks. }
    
    \item{\code{concurrent(by=NULL)}  (binary) (undirected) (categorical nodal attribute)
    }{\emph{Concurrent node count:} 
    This term adds one network statistic to the model, equal to the number of
    nodes in the network with degree 2 or higher. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list;
    it functions just like the \code{by} argument of the \code{degree} term.
    This term can only be used with undirected
    networks. }

    \item{\code{concurrentties(by=NULL)}  (binary) (undirected) (categorical nodal attribute)
    }{\emph{Concurrent tie count:} 
    This term adds one network statistic to the model, equal to the number of
    ties incident on each actor beyond the first. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list;
    it functions just like the \code{by} argument of the \code{degree} term.
    This term can only be used with undirected
    networks. }

    
    \item{\code{ctriple(attrname=NULL)}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    , 
    a.k.a. \code{ctriad}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    }{\emph{Cyclic triples:} This term adds one
    statistic to the model, equal to the number of cyclic triples in the
    network, defined as a set of edges of the form \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k), (k{\rightarrow}i)\}}{\{(i,j), (j,k), (k,i)\}}. Note that
    for all directed networks, \code{triangle} is equal to
    \code{ttriple+ctriple}, so at most two of these three terms can be in a
    model. The optional argument \code{attrname} is a character string giving
    the name of an attribute in the network's vertex attribute list. If this is
    specified then the count is over the number of cyclic triples where all
    three nodes have the same value of the attribute. This term can only be used
    with directed networks.}
    
    \item{\code{cycle(k)}  (binary)  (directed)  (undirected)}{\emph{Cycles:} 
    The \code{k} argument is a vector of distinct integers. This term adds one
    network statistic to the model for each element in \code{k}; the \eqn{i}th
    such statistic equals the number of cycles in the network with length
    exactly \code{k[i]}. The cycle statistic applies to both directed and
    undirected networks. For directed networks, it counts directed cycles of
    length \eqn{k}, as opposed to undirected cycles in the undirected case. The
    directed cycle terms of lengths 2 and 3 are equivalent to \code{mutual} and
    \code{ctriple} (respectively). The undirected cycle term of length 3 is
    equivalent to \code{triangle}, and there is no undirected cycle term of
    length 2.}

  \item{\code{cyclicalties(attrname=NULL)}   (binary)  (directed), \code{cyclicalties(threshold=0)}   (valued)  (directed)  (undirected)}{\emph{Cyclical ties:} 
    This term adds one statistic, equal to the number of ties 
    \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
    \eqn{i} to \eqn{j}.  (Related to the \code{ttriple} term.)
    The binary version takes a nodal attribute \code{attrname}, and, if
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
	to induce a multimodal distribution of networks.
  }

  \item{\code{degrange(from, to=+Inf, by=NULL, homophily=FALSE)} (binary)  (undirected) (categorical nodal attribute)
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
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
  
  \item{\code{degree(d, by=NULL, homophily=FALSE)}  (binary)  (undirected)  (categorical nodal attribute) (frequently-used)}{\emph{Degree:} 
    The \code{d} argument is a vector of distinct integers. This term adds one
    network statistic to the model for each element in \code{d}; the \eqn{i}th
    such statistic equals the number of nodes in the network of degree
    \code{d[i]}, i.e. with exactly \code{d[i]} edges. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with undirected networks; for directed networks
    see \code{idegree} and \code{odegree}. }

  \item{\code{degreepopularity}  (binary)  (undirected)
  }{\emph{Degree popularity:} This
    term adds one network statistic to the model equaling the sum over
    the actors of each actor's degree taken to the 3/2 power (or,
    equivalently, multiplied by its square root). This term is an
    undirected analog to the terms of Snijders et al. (2010), equations
    (11) and (12). This term can only be used with undirected networks.}
  
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
    divided by \eqn{n(n-1)}. }
    
    \item{\code{dsp(d)}  (binary)  (directed)  (undirected)}{\emph{Dyadwise shared partners:} 
    The \code{d} argument is a vector of distinct integers. This term adds one
    network statistic to the model for each element in \code{d}; the \eqn{i}th
    such statistic equals the number of dyads in the network with exactly
    \code{d[i]} shared partners. This term can be used with directed and
    undirected networks. For directed networks the count is over homogeneous
    shared partners only (i.e., only partners on a directed two-path connecting
    the nodes in the dyad).}
    
    \item{\code{dyadcov(x, attrname=NULL)}  (binary)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Dyadic covariate:}
    If the network is directed, \code{x} is either a (symmetric) matrix of
    covariates, one for each possible dyad \eqn{(i,j)}, or an undirected
    network; if the latter, optional argument \code{attrname} provides the name
    of the quantitative edge attribute to use for covariate values (in this
    case, missing edges in \code{x} are assigned a covariate value of zero).
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
    both \code{ostar(1)} and \code{istar(1)}. }
    
    \item{\code{esp(d)}  (binary)  (directed)  (undirected)}{\emph{Edgewise shared partners:} 
    This is just like the \code{dsp} term, except this term adds one network
    statistic to the model for each element in \code{d} where the \eqn{i}th such
    statistic equals the number of \emph{edges} (rather than dyads) in the
    network with exactly \code{d[i]} shared partners. This term can be used with
    directed and undirected networks. For directed networks the count is over
    homogeneous shared partners only (i.e., only partners on a directed two-path
    connecting the nodes in the edge and in the same direction).}

       
    \item{\code{greaterthan(threshold=0)}  (valued)  (directed)  (undirected) (dyadic-independent)}{\emph{Number of dyads with
	values strictly greater than a threshold:}
      Adds one statistic equaling to the number of ties whose values
      exceed \code{threshold}.
    }
  
  
    \item{\code{gwb1degree(decay, fixed=FALSE, cutoff=30)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
    degree distribution for the first mode in a bipartite (aka two-mode)
    network:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the \code{decay} parameter,
    for nodes in the
    first mode of a bipartite network. The first mode of a bipartite network
    object is sometimes known as the "actor" mode. 
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used as merely the starting value for the estimation
    in a curved exponential family model (the default).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.
    This term can only be used with undirected bipartite
    networks.}
    
    \item{\code{gwb2degree(decay, fixed=FALSE, cutoff=30)}  (binary)  (bipartite)  (undirected) (curved)}{\emph{Geometrically weighted
    degree distribution for the second mode in a bipartite (aka two-mode)
    network:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the \code{decay} parameter,
    for nodes in the
    second mode of a bipartite network. The second mode of a bipartite network
    object is sometimes known as the "event" mode. 
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used as merely the starting value for the estimation
    in a curved exponential family model (the default).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.
    This term can only be used with undirected bipartite
    networks.}
    
    \item{\code{gwdegree(decay, fixed=FALSE, cutoff=30)}  (binary)  (undirected)  (curved) (frequently-used)
    }{\emph{Geometrically weighted
    degree distribution:} 
    This term adds one network statistic to the model equal to the weighted
    degree distribution with decay controlled by the \code{decay} parameter.
    The \code{decay} parameter is the same as theta_s in
    equation (14) in Hunter (2007).  The value supplied for
    this parameter may be fixed (if \code{fixed=TRUE}),
    or it may be used as merely the starting value for the estimation
    in a curved exponential family model (the default).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.    This term can only be used with undirected networks.}

    \item{\code{gwdsp(alpha=0, fixed=FALSE, cutoff=30)}  (binary)  (directed)  (undirected)  (curved)}{\emph{Geometrically weighted
    dyadwise shared partner distribution:} 
    This term adds one network statistic to the model equal to the geometrically
    weighted dyadwise shared partner distribution with weight parameter
    \code{alpha} \eqn{> 0}. The optional argument \code{fixed} indicates whether
    the scale parameter \code{lambda} is to be fit as a curved exponential
    family model (see Hunter and Handcock, 2006). The default is \code{FALSE},
    which means the scale parameter is not fixed and thus the model is a CEF
    model. This term can be used with directed and undirected networks. For
    directed networks the count is over homogeneous shared partners only (i.e.,
    only partners on a directed two-path connecting the nodes in the dyad).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.}
    
    \item{\code{gwesp(alpha=0, fixed=FALSE, cutoff=30)}  (binary)  (frequently-used)   (directed)  (undirected)  (curved)}{\emph{Geometrically weighted
    edgewise shared partner distribution:} 
    This term is just like \code{gwdsp} except it adds a statistic equal to the
    geometrically weighted \emph{edgewise} (not dyadwise) shared partner
    distribution with weight parameter \code{alpha}. The optional argument
    \code{fixed} indicates whether the scale parameter \code{lambda} is to be
    fit as a curved exponential-family model (see Hunter and Handcock, 2006).
    The default is \code{FALSE}, which means the scale parameter is not fixed
    and thus the model is a CEF model. This term can be used with directed and
    undirected networks. For directed networks the geometric weighting is over
    homogeneous shared partners only (i.e., only partners on a directed two-path
    connecting the nodes in the edge and in the same direction).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.}
 
    \item{\code{gwidegree(decay, fixed=FALSE, cutoff=30)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted
    in-degree distribution:} This term adds one network statistic to the model
    equal to the weighted in-degree distribution with weight parameter
    \code{decay}. The optional argument \code{fixed} indicates whether the scale
    parameter \code{lambda} is to be fit as a curved exponential family model
    (see Hunter and Handcock, 2006). The default is \code{FALSE}, which means
    the scale parameter is not fixed and thus the model is a CEF model. This
    term can only be used with directed networks.
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.}

    \item{\code{gwnsp(alpha=0, fixed=FALSE, cutoff=30)}  (binary)  (directed)  (undirected) (curved)}{\emph{Geometrically weighted
    nonedgewise shared partner distribution:} This term is just like
    \code{gwesp} and \code{gwdsp} except it adds a statistic equal to
    the geometrically weighted \emph{nonedgewise} (that is, over dyads
    that do not have an edge) shared partner distribution with weight
    parameter \code{alpha}. The optional argument \code{fixed} indicates
    whether the scale parameter \code{lambda} is to be fit as a curved
    exponential-family model (see Hunter and Handcock, 2006).  The
    default is \code{FALSE}, which means the scale parameter is not
    fixed and thus the model is a CEF model. This term can be used with
    directed and undirected networks. For directed networks the
    geometric weighting is over homogeneous shared partners only (i.e.,
    only partners on a directed two-path connecting the nodes in the
    non-edge and in the same direction).
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.}
  
    \item{\code{gwodegree(decay, fixed=FALSE, cutoff=30)}  (binary)  (directed) (curved)}{\emph{Geometrically weighted
    out-degree distribution:} This term adds one network statistic to the model
    equal to the weighted out-degree distribution with weight parameter
    \code{decay}. The optional argument \code{fixed} indicates whether the scale
    parameter \code{lambda} is to be fit as a curved exponential family model
    (see Hunter and Handcock, 2006). The default is \code{FALSE}, which means
    the scale parameter is not fixed and thus the model is a CEF model. This
    term can only be used with directed networks.
    The optional argument \code{cutoff}
    is only relevant if \code{fixed=FALSE}. In that case it only uses this 
    number of terms in computing the statistics to reduce the computational
    burden.}
    
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
    
    \item{\code{hammingmix(attrname, x, base=0)}  (binary)  (directed)  (dyad-independent)}{\emph{Hamming 
    distance within mixing:} 
    This term adds one statistic to the model for every possible pairing of  
    attribute values of the network for the vertex attribute named \code{attrname}. 
    Each such statistic is the Hamming distance (i.e., the number of differences) 
    between the appropriate subset of dyads in
    the network and the corresponding subset in \code{x}. The ordering of the
    attribute values is alphabetical. 
    The option \code{base} gives the index of
    statistics to be omitted from the tabulation. For example \code{base=2} will
    omit the second statistic, making it the de facto reference category.
    This term can only be used with directed networks.}

  \item{\code{idegrange(from, to=+Inf, by=NULL, homophily=FALSE)}  (binary)  (directed) (categorical nodal attribute)
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
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

  
  \item{\code{idegree(d, by=NULL, homophily=FALSE)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{In-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of nodes in the network of in-degree \code{d[i]}, i.e. the number
    of nodes with exactly \code{d[i]} in-edges. The optional term
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with directed networks; for undirected networks
    see \code{degree}. }

  \item{\code{idegreepopularity}  (binary)  (directed)}{\emph{In-degree popularity:}
    This term adds one network statistic to the model equaling the sum
    over the actors of each actor's in-degree taken to the 3/2 power
    (or, equivalently, multiplied by its square root). This term is
    analogous to the term of Snijders et al. (2010), equation (11). This
    term can only be used with directed networks.}

  \item{\code{ininterval(lower=-Inf, upper=+Inf, open=c(TRUE,TRUE))}  (valued)  (directed)  (undirected) (dyadic-independent)
  }{\emph{Number of ties
	whose values are in an interval}
      Adds one statistic equaling to the number of ties whose values
      are between \code{lower} and
      \code{upper}. Argument \code{open}
      is a \code{logical} vector of length 2 that controls whether
      the interval is open (exclusive) on the lower and on the upper
      end, respectively.
    }

  
    \item{\code{intransitive}  (binary)  (directed)  (triad-related)}{\emph{Intransitive triads:} 
    This term adds one statistic to the model, equal to the number of triads in
    the network that are intransitive. The intransitive triads are those of type
    \code{111D}, \code{201}, \code{111U}, \code{021C}, or \code{030C} in the
    categorization of Davis and Leinhardt (1972). For details on the 16 possible
    triad types, see \code{\link[sna]{triad.classify}} in the
    \code{\link[sna]{sna}} package. Note the distinction from the \code{ctriple}
    term. This term can only be used with directed networks.}
    
    \item{\code{isolates}  (binary)  (directed)  (undirected) (frequently-used)
    }{\emph{Isolates:} This term adds one statistic to the
    model equal to the number of isolates in the network. For an undirected
    network, an isolate is defined to be any node with degree zero. For a 
    directed network, an isolate is any node with both in-degree and out-degree
    equal to zero.}
    
    \item{\code{istar(k, attrname=NULL)}  (binary)  (directed) (categorical nodal attribute)
    }{\emph{In-stars: } The \code{k} argument is a
    vector of distinct integers. This term adds one network statistic to the
    model for each element in \code{k}. The \eqn{i}th such statistic counts the
    number of distinct \code{k[i]}-instars in the network, where a
    \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}} such that the ties
    \eqn{(O_j{\rightarrow}N)}{(O_j, N)} exist for \eqn{j=1, \dots, k}. The
    optional argument \code{attrname} is a character string giving the name of
    an attribute in the network's vertex attribute list. If this is specified
    then the count is over the number of \eqn{k}-instars where all nodes have
    the same value of the attribute. This term can only be used for directed
    networks; for undirected networks see \code{kstar}. Note that
    \code{istar(1)} is equal to both \code{ostar(1)} and \code{edges}. }
    
    \item{\code{kstar(k, attrname=NULL)}  (binary)  (undirected) (categorical nodal attribute)
    }{\emph{k-Stars: } 
    The \code{k} argument is a vector of distinct integers. This term adds one
    network statistic to the model for each element in \code{k}. The \eqn{i}th
    such statistic counts the number of distinct \code{k[i]}-stars in the
    network, where a \eqn{k}-star is defined to be a node \eqn{N} and a set of
    \eqn{k} different nodes \eqn{\{O_1, \dots, O_k\}} such that the ties
    \eqn{\{N, O_i\}} exist for \eqn{i=1, \dots, k}. The optional argument
    \code{attrname} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified then the count is over
    the number of \eqn{k}-stars where all nodes have the same value of the
    attribute. This term can only be used for undirected networks; for directed
    networks, see \code{istar}, \code{ostar}, \code{twopath} and \code{m2star}.
    Note that \code{kstar(1)} is equal to \code{edges}. }
    
    \item{\code{localtriangle(x)}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Triangles within neighborhoods: } 
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
    special case of \code{localtriangle}. }
    
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
  
    \item{\code{mutual(same=NULL, diff=FALSE, by=NULL, keep=NULL)} (binary)  (directed)  (dyad-independent)  (frequently-used), \code{mutual(form="min",threshold=0)}  (valued)  (directed) (dyad-independent) 
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
    the name of a vertex attribute,
    only mutual pairs that match on the attribute are counted;
    separate counts for each unique matching value can be obtained by using
    \code{diff=TRUE} with \code{same}; and if \code{by} is passed the name of a vertex attribute,
    then each node is counted separately for each mutual pair in which it 
    occurs and the counts are tabulated by unique values of the attribute.
    This means that the sum of the mutual statistics when \code{by} is used
    will equal twice the standard mutual statistic.  Only one of \code{same}
    or \code{by} may be used, and only the former is affected by \code{diff};
    if both \code{same} and \code{by} are passed, \code{by} is ignored.
    Finally, if \code{keep} is passed a numerical vector, this vector of integers
    tells which statistics should be kept whenever the \code{mutual} term would
    ordinarily result in multiple statistics.
  }
    
    \item{\code{nearsimmelian}  (binary)  (directed)  (triad-related)}{\emph{Near simmelian triads:} 
    This term adds one statistic to the model equal to the number of near
    Simmelian triads, as defined by Krackhardt and Handcock (2007). This is a
    sub-graph of size three which is exactly one tie short of being complete.
    This term can only be used with directed networks.}
    
    \item{\code{nodecov(attrname)}  (binary)  (dyad-independent)  (frequently-used)   (directed)  (undirected) (quantitative nodal attribute)
    , \code{nodecov(attrname,
    form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (quantitative nodal attribute)
    , a.k.a. \code{nodemain}  (binary)   (directed)  (undirected)}{\emph{Main effect of a covariate:} 
    The \code{attrname} argument is a character string giving the name of a
    numeric (not categorical) attribute in the network's vertex attribute list.
    This term adds a single network statistic to the model equaling the sum of
    \code{attrname(i)} and \code{attrname(j)} for all edges \eqn{(i,j)} in the
    network. For categorical attributes, see \code{nodefactor}. Note that for
    directed networks, \code{nodecov} equals \code{nodeicov} plus
    \code{nodeocov}.}

    
    \item{\code{nodecovar}  (valued)  (directed)  (undirected) (quantitative nodal attribute)
    }{\emph{Uncentered covariance of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} (y_{i,j}y_{i,k}+y_{k,j}y_{k,j})}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic.
    }
  
    
    \item{\code{nodefactor(attrname, base=1)}  (binary)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute) (frequently-used) 
    ,
    \code{nodefactor(attrname, base=1, form="sum")}  (dyad-independent)  (valued)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Factor attribute
    effect:} The \code{attrname} argument is a character vector giving
    one or more names of categorical attributes in the network's vertex
    attribute list. This term adds multiple network statistics to the
    model, one for each of (a subset of) the unique values of the
    \code{attrname} attribute (or each combination of the attributes
    given). Each of these statistics gives the number of times a node
    with that attribute or those attributes appears in an edge in the
    network. In particular, for edges whose endpoints both have the same
    attribute values, this value is counted twice. To include all
    attribute values is usually not a good idea -- though this may be 
    accomplished if desired by setting \code{base=0} -- because the sum
    of all such statistics equals twice the number of edges and hence a 
    linear dependency would arise in any model also including
    \code{edges}. Thus, the \code{base} argument tells which value(s) 
    (numbered in order according to the \code{sort} function) should be
    omitted. The default value, \code{base=1}, means that the smallest 
    (i.e., first in sorted order) attribute value is omitted. For
    example, if the \dQuote{fruit} factor has levels \dQuote{orange},
    \dQuote{apple}, \dQuote{banana}, and \dQuote{pear}, then to add just
    two terms, one for \dQuote{apple} and one for \dQuote{pear}, then
    set \dQuote{banana} and \dQuote{orange} to the base (remember to
    sort the values first) by using \code{nodefactor("fruit",
    base=2:3)}. For an analogous term for quantitative vertex
    attributes, see \code{nodecov}. }
    
    \item{\code{nodeicov(attrname)}  (binary)  (directed) (quantitative nodal attribute) (frequently-used)
    , \code{nodeicov(attrname, form="sum")}  (valued)  (directed) (quantitative nodal attribute)
    }{\emph{Main effect of a covariate for
    in-edges:} 
    The \code{attrname} argument is a character string giving the name of a
    numeric (not categorical) attribute in the network's vertex attribute list.
    This term adds a single network statistic to the model equaling the total
    value of \code{attrname(j)} for all edges \eqn{(i,j)} in the network. This
    term may only be used with directed networks. For categorical attributes,
    see \code{nodeifactor}.}

      \item{\code{nodeicovar}  (valued)  (directed) (quantitative nodal attribute)
      }{\emph{Uncentered covariance of in-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{k,j}y_{k,j}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.
    }

  
    \item{\code{nodeifactor(attrname, base=1)}  (binary)  (dyad-independent)  (directed) (categorical nodal attribute) (frequently-used)
    ,
    \code{nodeifactor(attrname, base=1, form="sum")}  (valued)  (dyad-independent)  (directed) (categorical nodal attribute)
    }{\emph{Factor attribute
    effect for in-edges:} The \code{attrname} argument is a character
    vector giving one or more names of a categorical attribute in the
    network's vertex attribute list. This term adds multiple network
    statistics to the model, one for each of (a subset of) the unique
    values of the \code{attrname} attribute (or each combination of the
    attributes given). Each of these statistics gives the number of
    times a node with that attribute or those attributes appears as the
    terminal node of a directed tie. To include all attribute values is
    usually not a good idea -- though this may be accomplished if desired
    by setting \code{base=0} --
    because the sum of all such statistics
    equals the number of edges and hence a linear dependency would arise
    in any model also including \code{edges}. Thus, the \code{base}
    argument tells which value(s) (numbered in order according to the
    \code{sort} function) should be omitted. The default value,
    \code{base=1}, means that the smallest (i.e., first in sorted order)
    attribute value is omitted. For example, if the \dQuote{fruit}
    factor has levels \dQuote{orange}, \dQuote{apple}, \dQuote{banana},
    and \dQuote{pear}, then to add just two terms, one for
    \dQuote{apple} and one for \dQuote{pear}, then set \dQuote{banana}
    and \dQuote{orange} to the base (remember to sort the values first)
    by using \code{nodefactor("fruit", base=2:3)}. For an analogous term
    for quantitative vertex attributes, see \code{nodeicov}. }

      \item{\code{nodeisqrtcovar}  (valued)  (directed)  (non-negative) (quantitative nodal attribute)
      }{\emph{Uncentered covariance of square roots of in-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{k,j}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{istar(2)}}
      statistic.
    }

  
    \item{\code{nodematch(attrname, diff=FALSE,
    keep=NULL)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (categorical nodal attribute)
    , \code{nodematch(attrname, diff=FALSE,
    keep=NULL, form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
     a.k.a. \code{match}  (binary)  (directed)  (dyad-independent)  (undirected) (categorical nodal attribute)
     }{\emph{Uniform homophily and differential homophily:}
    The \code{attrname} argument is a character vector giving one or
    more names of attributes in the network's vertex attribute
    list. When \code{diff=FALSE}, this term adds one network statistic
    to the model, which counts the number of edges \eqn{(i,j)} for which
    \code{attrname(i)==attrname(j)}. (When multiple names are given, the
    statistic counts only those on which all the named attributes
    match.) When \code{diff=TRUE}, \eqn{p} network statistics are added
    to the model, where \eqn{p} is the number of unique values of the
    \code{attrname} attribute. The \eqn{k}th such statistic counts the
    number of edges \eqn{(i,j)} for which \code{attrname(i) ==
    attrname(j) == value(k)}, where \code{value(k)} is the \eqn{k}th
    smallest unique value of the attrname attribute. If set to non-NULL,
    the optional \code{keep} argument should be a vector of integers
    giving the values of \code{k} that should be considered for matches;
    other values are ignored (this works for both \code{diff=FALSE} and
    \code{diff=TRUE}). For instance, to add two statistics, counting the
    matches for just the 2nd and 4th categories, use \code{nodematch}
    with \code{diff=TRUE} and \code{keep=c(2,4)}.}
   
    \item{\code{nodemix(attrname, base=NULL)}  (binary)  (dyad-independent)  (frequently-used)  (directed)  (undirected) (categorical nodal attribute)
    ,
    \code{nodemix(attrname, base=NULL, form="sum")}  (valued)  (dyad-independent)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Nodal attribute
    mixing:} The \code{attrname} argument is a character vector giving
    the names of categorical attributes in the network's vertex
    attribute list. By default, this term adds one network statistic to
    the model for each possible pairing of attribute values. The
    statistic equals the number of edges in the network in which the
    nodes have that pairing of values. (When multiple names are given, a
    statistic is added for each combination of attribute values for
    those names.) In other words, this term produces one statistic for
    every entry in the mixing matrix for the attribute(s). The ordering of
    the attribute values is alphabetical (for nominal categories) or
    numerical (for ordered categories). The optional \code{base}
    argument is a vector of integers corresponding to the pairings that
    should not be included. If \code{base} contains only negative
    integers, then these integers correspond to the only pairings that
    should be included. By default (i.e., with \code{base=NULL} or
    \code{base=0}), all pairings are included.}
    
    \item{\code{nodeocov(attrname)}  (binary)  (directed) (dyadic-independent)(quantitative nodal attribute)
    , \code{nodeocov(attrname, form="sum")}  (valued)  (directed) (dyadic-independent)
	(quantitative nodal attribute)
    }{\emph{Main effect of a covariate for
    out-edges:} 
    The \code{attrname} argument is a character string giving the name of a
    numeric (not categorical) attribute in the network's vertex attribute list.
    This term adds a single network statistic to the model equaling the total
    value of \code{attrname(i)} for all edges \eqn{(i,j)} in the network. This
    term may only be used with directed networks. For categorical attributes,
    see \code{nodeofactor}.}

    \item{\code{nodeocovar}  (valued)  (directed) (quantitative nodal attribute)
    }{\emph{Uncentered covariance of out-dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k} y_{i,j}y_{i,k}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.
    }

  
    \item{\code{nodeofactor(attrname, base=1)}  (binary)  (dyad-independent)  (directed) (categorical nodal attribute)
    ,
    \code{nodeofactor(attrname, base=1, form="sum")}  (valued)  (dyad-independent)  (categorical nodal attribute) (directed)}{\emph{Factor attribute
    effect for out-edges:} The \code{attrname} argument is a character
    string giving one or more names of categorical attributes in the
    network's vertex attribute list. This term adds multiple network
    statistics to the model, one for each of (a subset of) the unique
    values of the \code{attrname} attribute (or each combination of the
    attributes given). Each of these statistics gives the number of
    times a node with that attribute or those attributes appears as the
    node of origin of a directed tie. To include all attribute values is
    usually not a good idea -- though this may be accomplished if desired
    by setting \code{base=0} --
    because the sum of all such statistics
    equals the number of edges and hence a linear dependency would arise
    in any model also including \code{edges}. Thus, the \code{base}
    argument tells which value(s) (numbered in order according to the
    \code{sort} function) should be omitted. The default value,
    \code{base=1}, means that the smallest (i.e., first in sorted order)
    attribute value is omitted. For example, if the \dQuote{fruit}
    factor has levels \dQuote{orange}, \dQuote{apple}, \dQuote{banana},
    and \dQuote{pear}, then to add just two terms, one for
    \dQuote{apple} and one for \dQuote{pear}, then set \dQuote{banana}
    and \dQuote{orange} to the base (remember to sort the values first)
    by using \code{nodefactor("fruit", base=2:3)}. For an analogous term
    for quantitative vertex attributes, see \code{nodeocov}. }

      \item{\code{nodeosqrtcovar}  (valued)  (directed)  (non-negative)
	  (quantitative nodal attribute)
      }{\emph{Uncentered covariance of square roots of out-dyad values
	incident on each actor:} This term adds one statistic equal to \eqn{\sum_{i,j,k}
	\sqrt{y_{i,j}}\sqrt{y_{i,k}}}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{ostar(2)}}
      statistic.
    }

    \item{\code{nodesqrtcovar(center=TRUE)}  (valued)  (non-negative)  (directed)  (undirected) (quantitative nodal attribute)
    }{\emph{Covariance of square roots of dyad values
	incident on each actor:} This term adds one statistic equal to
      \eqn{\sum_{i,j,k}
	(\sqrt{y_{i,j}}\sqrt{y_{i,k}}+\sqrt{y_{k,j}}\sqrt{y_{k,j}})} if
      \code{center=FALSE}. This can be
      viewed as a valued analog of the \code{\link[=ergm-terms]{kstar(2)}}
      statistic. If \code{center=FALSE} (the default), the statistic is
	instead \eqn{\sum_{i,j,k}
	((\sqrt{y_{i,j}}-\bar{\sqrt{y}})(\sqrt{y_{i,k}}-\bar{\sqrt{y}})+(\sqrt{y_{k,j}}-\bar{\sqrt{y}})(\sqrt{y_{k,j}}-\bar{\sqrt{y}}))},
	where \eqn{\bar{\sqrt{y}}} is the mean of the square root of
	dyad values.
    }

  
    \item{\code{nsp(d)}  (binary)  (directed)  (undirected)}{\emph{Nonedgewise shared partners:} This is
    just like the \code{dsp} and \code{esp} terms, except this term adds
    one network statistic to the model for each element in \code{d}
    where the \eqn{i}th such statistic equals the number of 
    \emph{non-edges} (that is, dyads that do not have an edge) in the network
    with exactly \code{d[i]} shared partners. This term can be used with
    directed and undirected networks. For directed networks the count is
    over homogeneous shared partners only (i.e., only partners on a
    directed two-path connecting the nodes in the non-edge and in the same
    direction).}

    \item{\code{odegrange(from, to=+Inf, by=NULL, homophily=FALSE)} (binary)  (directed) (categorical nodal attribute)
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
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
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
  
    \item{\code{odegree(d, by=NULL, homophily=FALSE)}  (binary)  (directed)  (categorical nodal attribute) (frequently-used)}{\emph{Out-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of nodes in the network of out-degree \code{d[i]}, i.e. the
    number of nodes with exactly \code{d[i]} out-edges. The optional argument
    \code{by} is a character string giving the name of an attribute in the
    network's vertex attribute list. 
    If this is specified and \code{homophily} is \code{TRUE},
    then degrees are calculated using the subnetwork consisting of only
    edges whose endpoints have the same value of the \code{by} attribute.
    If \code{by} is specified and
    \code{homophily} is \code{FALSE} (the default), then separate degree 
    statistics are calculated for nodes having each separate
    value of the attribute.
    This term can only be used with directed networks; for undirected networks
    see \code{degree}. }
  
  \item{\code{odegreepopularity}  (binary)  (directed)
  }{\emph{Out-degree popularity:}
    This term adds one network statistic to the model equaling the sum
    over the actors of each actor's outdegree taken to the 3/2 power
    (or, equivalently, multiplied by its square root). This term is
    analogous to the term of Snijders et al. (2010), equation (12). This
    term can only be used with directed networks.}

    \item{\code{opentriad}  (binary)  (undirected) (triad-related)}{\emph{Open triads:} This term
    adds one statistic to the model equal to the number of 2-stars minus
    three times the number of triangles in the network. It is currently
    only implemented for undirected networks. }
  
    \item{\code{ostar(k, attrname=NULL)}  (binary)  (directed) (categorical nodal attribute)
    }{\emph{k-Outstars:} The \code{k} argument is
    a vector of distinct integers. This term adds one network statistic to the
    model for each element in \code{k}. The \eqn{i}th such statistic counts the
    number of distinct \code{k[i]}-outstars in the network, where a
    \eqn{k}-outstar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}} such that the ties
    \eqn{(N{\rightarrow}O_j)}{(N,O_j)} exist for \eqn{j=1, \dots, k}. The
    optional argument \code{attrname} is a character string giving the name of
    an attribute in the network's vertex attribute list. If this is specified
    then the count is the number of \eqn{k}-outstars where all nodes have the
    same value of the attribute. This term can only be used with directed
    networks; for undirected networks see \code{kstar}. Note that
    \code{ostar(1)} is equal to both \code{istar(1)} and \code{edges}. }
    
    \item{\code{receiver(base=1)}  (binary)  (directed)  (dyad-independent)}{\emph{Receiver effect:} 
    This term adds one network statistic for each node equal to the number of
    in-ties for that node. This measures the popularity of the node. The term
    for the first node is omitted by default because of linear dependence that
    arises if this term is used together with \code{edges}, but its coefficient
    can be computed as the negative of the sum of the coefficients of all the
    other actors. That is, the average coefficient is zero, following the
    Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
    1981). The \code{base} argument allows the user to determine which nodes'
    statistics should be omitted. The \code{base} argument can also be a vector
    of negative indices, to specify which should be added instead of deleted,
    and \code{base=0} specifies that all statistics should be included. This
    term can only be used with directed networks. For undirected networks, see
    \code{sociality}.}
    
    \item{\code{sender(base=1)}  (binary)  (directed)  (dyad-independent)}{\emph{Sender effect:} 
    This term adds one network statistic for each node equal to the number of
    out-ties for that node. This measures the activity of the node. The term for
    the first node is omitted by default because of linear dependence that
    arises if this term is used together with \code{edges}, but its coefficient
    can be computed as the negative of the sum of the coefficients of all the
    other actors. That is, the average coefficient is zero, following the
    Holland-Leinhardt parametrization of the $p_1$ model (Holland and Leinhardt,
    1981). The \code{base} argument allows the user to determine which nodes'
    statistics should be omitted. The \code{base} argument can also be a vector
    of negative indices, to specify which should be added instead of deleted,
    and \code{base=0} specifies that all statistics should be included. This
    term can only be used with directed networks. For undirected networks, see
    \code{sociality}.}

    \item{\code{simmelian}  (binary)  (directed)  (triad-related)}{\emph{Simmelian triads:} This term adds one
    statistic to the model equal to the number of Simmelian triads, as defined
    by Krackhardt and Handcock (2007). This is a complete sub-graph of size
    three. This term can only be used with directed networks.}
    
    \item{\code{simmelianties}  (binary)  (triad-related)  (directed)}{\emph{Ties in simmelian triads:} This term adds
    one statistic to the model equal to the number of ties in the network that
    are associated with Simmelian triads, as defined by Krackhardt and Handcock
    (2007). Each Simmelian has six ties in it but, because Simmelians can
    overlap in terms of nodes (and associated ties), the total number of ties in
    these Simmelians is less than six times the number of Simmelians. Hence this
    is a measure of the clustering of Simmelians (given the number of
    Simmelians). This term can only be used with directed networks.}

  \item{\code{smalldiff(attrname, cutoff)}  (binary)  (dyad-independent)  (directed)  (undirected) (quantitative nodal attribute)
  }{\emph{Number of
    ties between actors with similar (but not necessarily identical)
    attribute values:} The \code{attrname} argument is a character
    string giving the name of a quantitative attribute in the network's
    vertex attribute list. This term adds one statistic, having as its
    value the number of edges in the network for which the incident
    actors' attribute values differ less than \code{cotoff}; that is,
    number of edges between \code{i} to \code{j} such that
    \code{abs(attrname[i]-attrname[j])<cutoff}.}
     
    \item{\code{sociality(attrname=NULL, base=1)}  (binary)  (undirected) (categorical nodal attribute)
    }{\emph{Undirected degree:} 
    This term adds one network statistic for each node equal to the number of
    ties of that node. The optional \code{attrname} argument is a character
    string giving the name of an attribute in the network's vertex attribute
    list that takes categorical values. If provided, this term only counts ties
    between nodes with the same value of the attribute (an actor-specific
    version of the \code{nodematch} term). This term can only be used with
    undirected networks.  For directed networks, see \code{sender} and
    \code{receiver}. By default, \code{base=1} means that the statistic for the
    first node will be omitted, but this argument may be changed to control
    which statistics are included just as for the \code{sender} and
    \code{receiver} terms.}


    \item{\code{sum(pow=1)}  (valued)  (directed)  (undirected)}{\emph{Sum of dyad values (optionally taken
	to a power):} This term adds one statistic equal to the sum of
      dyad values taken to the power \code{pow}, which defaults to 1.
    }
      
  
    \item{\code{threepath(keep=1:4)}  (binary)  (directed)  (undirected) (triad-related)
    }{\emph{ Three-paths:}
    For an undirected network, this term adds one statistic equal to the number
    of threepaths, where a threepath is defined as a path of length three that
    traverses three distinct edges.
    Note that a threepath need not
    include four distinct nodes; in particular, a triangle counts as three
    threepaths.  For a directed network, this term adds four statistics 
    (or some subset of these four specified by the \code{keep} argument),
    one for each of the four distinct types of directed three-paths.  If the 
    nodes of the path are written from left to right such that the middle edge
    points to the right (R), then the four types are RRR, RRL, LRR, and LRL.
    That is, an RRR threepath is of the form 
    \eqn{i\rightarrow j\rightarrow k\rightarrow l}{i-->j-->k-->l}, and RRL
    threepath is of the form
    \eqn{i\rightarrow j\rightarrow k\leftarrow l}{i-->j-->k<--l}, etc.
    Like in the undirected case, there is no requirement that the nodes be
    distinct in a directed threepath.  However, the three edges must all be 
    distinct.  Thus, a mutual tie \eqn{i\leftrightarrow j}{i<-->j} does not
    count as a threepath of the form
    \eqn{i\rightarrow j\rightarrow i\leftarrow j}{i-->j-->i<--j}; however,
    in the subnetwork \eqn{i\leftrightarrow j \rightarrow k}{i<-->j-->k},
    there are two directed threepaths, one LRR
    (\eqn{k\leftarrow j\rightarrow i\leftarrow j}{k<--j-->i-->j})
    and one RRR
    (\eqn{j\rightarrow i\rightarrow j\leftarrow k}{k<--j-->i-->j}).
    }

    \item{\code{transitive}  (binary)  (directed)  (triad-related)}{\emph{Transitive triads:} 
    This term adds one statistic to the model, equal to the number of triads in
    the network that are transitive. The transitive triads are those of type
    \code{120D}, \code{030T}, \code{120U}, or \code{300} in the categorization
    of Davis and Leinhardt (1972). For details on the 16 possible triad types,
    see \code{\link[sna]{triad.classify}} in the \code{\link[sna]{sna}} package.
    Note the distinction from the \code{ttriple} term. This term can only be
    used with directed networks.}
    
    \item{\code{transitiveties(attrname=NULL)}   (binary)  (directed) (triad-related) (categorical nodal attribute)
    , \code{transitiveties(threshold=0)}   (valued)  (directed)  (undirected) (triad-related)}{\emph{Transitive ties:} 
    This term adds one statistic, equal to the number of ties 
    \eqn{i\rightarrow j}{i-->j} such that there exists a two-path from
    \eqn{i} to \eqn{j}.  (Related to the \code{ttriple} term.)
    The binary version takes a nodal attribute \code{attrname}, and, if
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
	to induce a multimodal distribution of networks.
      }

    
    \item{\code{triadcensus(d)}  (binary)  (triad-related)  (directed)  (undirected)}{\emph{Triad census:} 
    For a directed network, this term adds one network statistic for each of
    an arbitrary subset of the 16 possible types of triads categorized by
    Davis and Leinhardt (1972) as \code{003, 012, 102, 021D, 021U, 021C, 111D,
    111U, 030T, 030C, 201, 120D, 120U, 120C, 210,} and \code{300}. Note that at
    least one category should be dropped; otherwise a linear dependency will
    exist among the 16 statistics, since they must sum to the total number of
    three-node sets. By default, the category \code{003}, which is the category
    of completely empty three-node sets, is dropped. This is considered category
    zero, and the others are numbered 1 through 15 in the order given above. By
    specifying a numeric vector of integers from 0 to 15 as the \code{d}
    argument, the user may specify a set of terms to add other than the default
    value of \code{1:15}. Each statistic is the count of the corresponding triad
    type in the network. For details on the 16 types, see \code{?triad.classify}
    in the \code{{sna}} package, on which this code is based. For an undirected
    network, the triad census is over the four types defined by the number of
    ties (i.e., 0, 1, 2, and 3), and the default is to add \code{1:3}, which is
    to say that the 0 is dropped; however, this too may be controlled by
    changing the \code{d} argument to a numeric vector giving a subset of
    \eqn{\{0, 1, 2, 3\}}.}
    
    \item{\code{triangle(attrname=NULL)}  (binary)  (frequently-used)  (triad-related)  (directed)  (undirected) (categorical nodal attribute)
    }{\emph{Triangles:}
    This term adds one statistic to the model equal to the number of triangles
    in the network. For an undirected network, a triangle is defined to be any
    set \eqn{\{(i,j), (j,k), (k,i)\}} of three edges. For a directed network, a
    triangle is defined as any set of three edges \eqn{(i{\rightarrow}j)}{(i,j)}
    and \eqn{(j{\rightarrow}k)}{(j,k)} and either \eqn{(k{\rightarrow}i)}{(k,i)}
    or \eqn{(k{\leftarrow}i)}{(i,k)}. The former case is called a ``transitive
    triple'' and the latter is called a ``cyclic triple'', so in the case of a
    directed network, \code{triangle} equals \code{ttriple} plus \code{ctriple}
    --- thus at most two of these three terms can be in a model. The optional
    argument \code{attrname} restricts the count to those triples of nodes with
    equal values of the vertex attribute specified by \code{attrname}. }
    
    \item{\code{tripercent(attrname=NULL)}  (binary)  (undirected)  (triad-related) (categorical nodal attribute)
    }{\emph{Triangle percentage:} 
    This term adds one statistic to the model equal to 100 times the ratio of
    the number of triangles in the network to the sum of the number of triangles
    and the number of 2-stars not in triangles (the latter is considered a
    potential but incomplete triangle).  In case the denominator equals zero,
    the statistic is defined to be zero. For the definition of triangle, see
    \code{triangle}. The optional argument \code{attrname} restricts the counts
    (both numerator and denominator) to those triples of nodes with equal values
    of the vertex attribute specified by \code{attrname}. This is often called
    the mean correlation coefficient. This term can only be
    used with undirected networks; for directed networks, it is difficult to
    define the numerator and denominator in a consistent and meaningful way.}
    
    \item{\code{ttriple(attrname=NULL)}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    , a.k.a. \code{ttriad}  (binary)  (directed)  (triad-related) (categorical nodal attribute)
    }{\emph{Transitive triples:} 
    This term adds one statistic to the model, equal to the number of transitive
    triples in the network, defined as a set of edges \eqn{\{(i{\rightarrow}j),
    (j{\rightarrow}k), (i{\rightarrow}k)\}}{\{(i,j), (j,k), (i,k)\}}. Note that
    \code{triangle} equals \code{ttriple+ctriple} for a directed network, so at
    most two of the three terms can be in a model. The optional argument
    \code{attrname} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified then the count is over
    the number of transitive triples where all three nodes have the same value
    of the attribute. This term can only be used with directed networks.}
    
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
    1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

    \item Snijders T. A. B., G. G. van de Bunt, and
    C. E. G. Steglich. Introduction to Stochastic Actor-Based Models for
    Network Dynamics. \emph{Social Networks}, 2010, 32(1),
    44-60. \href{http://dx.doi.org/10.1016/j.socnet.2009.02.004}{\code{doi:10.1016/j.socnet.2009.02.004}}

    \item Morris M, Handcock MS, and Hunter DR. Specification of Exponential-Family Random Graph Models: Terms and Computational Aspects. \emph{Journal of Statistical Software}, 2008, 24(4), 1-24. \url{http://www.jstatsoft.org/v24/i04}
    
    
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
