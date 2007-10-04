\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\docType{methods}
\title{Terms used in Exponential Family Random Graph Models}
\description{
  The function \code{\link{ergm}} is used to fit linear exponential 
  random graph models, in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{\exp(\theta{\cdot}g(y))/c(\theta)}, where 
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\theta} is a parameter vector of the same 
  length and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  The network statistics \eqn{g(y)} are entered as terms in the 
  function call to \code{\link{ergm}}.

  This page describes the possible terms (and hence network statistics).
}
%\usage{
%ergm-terms()
%}
%\arguments{
%}
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

\section{Possible terms to represent network statistics}{
  The \code{\link{ergm}} function allows the user to explore a large number
  of potential models for their network data.  What follows
  is a list of model terms currently available by the program,
  and a brief description of each.
  In the formula for the model, the model terms are various function-like
  calls, some of which require arguments, separated by \code{+} signs.

  Additional terms can be coded up by users via 
  the \code{statnetuserterms} package.

  The terms currently available are:

    \item{\code{absdiff(attrname)}}{\emph{Absolute Difference:} The
    \code{attrname} argument is a character string giving the name of an
    attribute in the network's vertex attribute list. This term adds one network
    statistic to the model equaling the sum of
    \code{abs(attrname[i]-attrname[j])} for all edges (i,j) in the network. }

    \item{\code{altkstar(lambda, fixed=FALSE)}}{\emph{Alternating k-star:} This
    term adds one network statistic to the model equal to a weighted alternating
    sequence of k-star statistics with weight parameter \code{lambda}. This is
    the version given in Snijders, et. al (2006). We suggest using the
    \code{gwdegree} term instead. The \code{gwdegree} and \code{altkstar}
    produce mathematically equivalent models, as long as they are used together
    with the \code{edges} (or \code{kstar(1)}) term, yet the interpretation of
    the \code{gwdegree} parameters is slightly more straightforward than the
    interpretation of the \code{altkstar} parameters. For this reason, we
    recommend the use of the \code{gwdegree} instead of \code{altkstar}. See
    Section 3 and especially equation (13) of
    \url{http://www.sna.unimelb.edu.au/publications/cef4.pdf} for details. The
    optional argument \code{fixed} is a logical indicating if the scale
    parameter \code{lambda} is to be fit as a curved exponential family model
    (see Hunter and Handcock, 2006). The default is \code{TRUE}. This term can
    only be used with undirected networks.}

    \item{\code{asymmetric}}{\emph{Asymmetric Dyads:} This term adds one network
    statistic to the model, equaling the number of pairs of actors for which
    exactly one of \eqn{(i{\rightarrow}j)}{(i,j)} or
    \eqn{(j{\rightarrow}i)}{(j,i)} exists. This term can only be used with
    directed networks. }

    \item{\code{bounded.degree(d, bound)}}{\emph{Bounded Degree:} The \code{d}
    argument is a vector of distinct integers representing degrees. \code{bound}
    is a vector of upper bounds corresponding to each degree. This term adds one
    network statistic to the model for each element in \code{d}; the \eqn{i}th
    such statistic equals the minimum of \code{bound[i]} and the number of nodes
    in the network of degree exactly \code{d[i]}, i.e. with exactly \code{d[i]}
    edges. This term can only be used with undirected networks; for directed
    networks see \code{bounded.idegree} and \code{bounded.odegree}. }

    \item{\code{bounded.idegree(d, bound)}}{\emph{Bounded In-Degree:} The
    \code{d} argument is a vector of distinct integers representing in-degrees.
    \code{bound} is a vector of upper bounds corresponding to each degree. This
    term adds one network statistic to the model for each element in \code{d};
    the \eqn{i}th such statistic equals the minimum of \code{bound[i]} and the
    number of nodes in the network of in-degree exactly \code{d[i]}, i.e. with
    exactly \code{d[i]} edges. This term can only be used with directed
    networks; for undirected networks see \code{bounded.degree}. }

    \item{\code{bounded.istar(k, bound)}}{\emph{Bounded In-stars:} The \code{k}
    argument is a vector of distinct integers. \code{bound} is a vector of upper
    bounds corresponding to each degree. This term adds one network statistic to
    the model for each element in \code{k}. The \eqn{i}th such statistic counts
    the number of distinct \code{k[i]}-instars in the network, where a
    \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}}{O_1, \dots, O_k} such that the
    ties \eqn{\{(N{\lefttarrow}O_j)\}}{(O_j,N)} all exist for \eqn{j=1, \dots,
    k}. The value is the minimum of this count and \code{bound[i]}. This term
    can only be specified for directed networks; for undirected networks, see
    \code{bounded.kstar}. }

    \item{\code{bounded.kstar(k, bound)}}{\emph{Bounded k-stars:} The \code{k}
    argument is a vector of distinct integers. \code{bound} is a vector of upper
    bounds corresponding to each degree. This term adds one network statistic to
    the model for each element in \code{k}. The \eqn{i}th such statistic counts
    the number of distinct \code{k[i]}-stars in the network, where a
    \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}}{O_1, \dots, O_k} such that the tie
    \eqn{\{N,O_j\}} exists for \eqn{j=1, \dots, k}. The value is the minimum of
    this count and \code{bound[i]}. This term can only be specified for
    undirected networks; for directed networks, see \code{bounded.istar},
    \code{bounded.ostar}.}
    
    \item{\code{bounded.odegree(d, bound)}}{\emph{Bounded Out-Degree:} The
    \code{d} argument is a vector of distinct integers representing in-degrees.
    \code{bound} is a vector of upper bounds corresponding to each degree. This
    term adds one network statistic to the model for each element in \code{d};
    the \eqn{i}th such statistic equals the minimum of \code{bound[i]} and the
    number of nodes in the network of out-degree exactly \code{d[i]}, i.e. with
    exactly \code{d[i]} edges. This term can only be used with directed
    networks; for undirected networks see \code{bounded.degree}. }

    \item{\code{bounded.ostar(k, bound)}}{\emph{Bounded Out-stars:} The \code{k}
    argument is a vector of distinct integers. \code{bound} is a vector of upper
    bounds corresponding to each degree. This term adds one network statistic to
    the model for each element in \code{k}. The \eqn{i}th such statistic counts
    the number of distinct \code{k[i]}-outstars in the network, where a
    \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}}{O_1, \dots, O_k} such that the
    ties \eqn{\{(N{\rightarrow}O_j)\}}{(N,O_j)} all exist for \eqn{j=1, \dots,
    k}. The value is the minimum of this count and \code{bound[i]}. This term
    can only be specified for directed networks; for undirected networks, see
    \code{bounded.kstar}. }

    \item{\code{bounded.triangle(bound)}}{\emph{Bounded Triangles:} This term
    adds one statistic to the model equal to the minimum of the bound and the
    number of triangles in the network. For an undirected network, a triangle is
    defined to be any set of three edges \eqn{\{i,j\}, \{j,k\}, and \{k,i\}}.
    For a directed network, a triangle is defined as any set of three edges
    \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k)\}}{(i,j), (j,k)} and either
    \eqn{(k{\rightarrow}i)}{(k,i)} or \eqn{(i{\rightarrow}k)}{(i,k)}.}

    \item{\code{ctriad(attrname)}}{\emph{Cyclic Triads:} This term adds one
    statistic to the model, equal to the number of cyclic triads in the network,
    defined as a set of edges \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k),
    (k{\rightarrow}i)\}}{\{(i,j), (j,k), (k,i)\}}. Note that for all directed
    networks, \code{triangle} is equal to \code{ttriad+ctriad}, so at most two
    of these three terms can be in a model. The argument \code{attrname} is a
    character string giving the name of an attribute in the network's vertex
    attribute list. If this is specified then the count is over the number of
    cyclic triads where all three nodes have the same value of the attribute.
    This term can only be used with directed networks.}

    \item{\code{cycle(k)}}{\emph{Cycles:} The \code{k} argument is a vector of
    distinct integers. This term adds one network statistic to the model for
    each element in \code{k}; the \eqn{i}th such statistic equals the number of
    cycles in the network with length exactly \code{k[i]}. The cycle statistic
    applies to both directed and undirected graphs. For directed networks, it
    counts directed cycles of length \eqn{k}, as opposed to undirected cycles in
    the undirected case. The directed cycle terms of lengths 2 and 3 are
    equivalent to \code{mutual} and \code{ctriad} (respectively). The undirected
    cycle term of length 3 is equivalent to \code{triangle}, and there is no
    undirected cycle term of length 2. }

    \item{\code{degree(d, attrname)}}{\emph{Degree:} The \code{d} argument is a
    vector of distinct integers. This term adds one network statistic to the
    model for each element in \code{d}; the \eqn{i}th such statistic equals the
    number of nodes in the network of degree \code{d[i]}, i.e. with exactly
    \code{d[i]} edges. The term \code{attrname} is a character string giving the
    name of an attribute in the network's vertex attribute list. If this is
    specified then the degree count is the number of nodes with the same value
    of the attribute as the ego node. This term can only be used with
    undirected networks; for directed networks see \code{idegree} and
    \code{odegree}. }

    \item{\code{density}}{\emph{Density:} This term adds one network statistic
    equal to the density of the network. For undirected networks, \code{density}
    equals \code{kstar(1)} or \code{edges} divided by \eqn{n(n-1)/2}; for
    directed networks, \code{density} equals \code{edges} or \code{istar(1)} or
    \code{ostar(1)} divided by \eqn{n(n-1)}. }
 
%   \item{\code{duration}}{\emph{Duration:}  
%   This term adds network statistics to model duration times.}

    \item{\code{dsp(d)}}{\emph{Dyadwise Shared Partners:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of dyads in the network with exactly \code{d[i]} shared partners.
    This term can only be used with undirected networks. }

    \item{\code{dyadcov(x, attrname)}}{\emph{Dyadic Covariate:} If the network
    is directed \code{x} is either a (symmetric) matrix of dyadic covariates, or
    an undirected network; if the latter, optional argument \code{attrname}
    provides the name of the edge attribute to use for edge values. This term
    adds three statistics to the model, representing the (polytomous) effect of
    the given covariate on the four possible dyad states (i.e., null, out-tie,
    in-tie, mutual). The statistics are the appearance of mutual,
    upper-triangular asymmetric, and lower-triangular asymmetric dyads (with the
    null state serving as a reference category). If the network is undirected,
    \code{x} is either a matrix of edgewise covariates, or a network; if the
    latter, optional argument \code{attrname} provides the name of the edge
    attribute to use for edge values. This term adds one statistic to the model,
    representing the effect of the given covariate on the appearance of edges.
    \code{dyadcov} can be called more than once, to model the effects of
    multiple covariates. }

%   \item{\code{dyadundircov(x, attrname=NULL)}}{\emph{Edge
%   Covariate.}  \code{x} is either a matrix of edgewise
%   covariates, or a network; if the latter, optional argument
%   \code{attrname} provides the name of the edge attribute to
%   use for edge values.  This term adds one statistic to the
%   model, representing the effect of the given covariate on
%   the appearance of edges.  \code{edgecov} can be called more
%   than once, to model the effects of multiple covariates.  }

    \item{\code{edgecov(x, attrname=NULL)}}{\emph{Edge Covariate:} The \code{x}
    argument is either a matrix of edgewise covariates, or a network; if the
    latter, optional argument \code{attrname} provides the name of the edge
    attribute to use for edge values. This term adds one statistic to the model,
    representing the effect of the given covariate on the appearance of edges.
    \code{edgecov} applies to both directed and undirected networks. For
    undirected networks the covariates are also assumed to be undirected.
    \code{edgecov} can be called more than once, to model the effects of
    multiple covariates. }

    \item{\code{edges}}{\emph{Edges:} This term adds one network statistic equal
    to the number of edges in the network. For undirected networks, \code{edges}
    is equal to \code{kstar(1)}; for directed networks, \code{edges} is equal to
    both \code{ostar(1)} and \code{istar(1)}. }

    \item{\code{esp(d)}}{\emph{Edgewise Shared Partners:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic equals
    the number of edges in the network with exactly \code{d[i]} shared partners.
    This term can only be used with undirected networks. }

    \item{\code{gwdegree(decay, fixed=FALSE)}}{\emph{Geometrically Weighted
    Degree Distribution:} This term adds one network statistic to the model
    equal to the weighted degree distribution with weight parameter
    \code{decay}. This is the version given as equation (14) in
    \url{http://www.sna.unimelb.edu.au/publications/cef4.pdf}. See the ``Remark"
    in section 3 of that paper to see why it is used rather than the version
    given in Snijders, et. al (2006). The optional argument \code{fixed} is a
    logical indicating if the decay parameter is to be fit as a curved
    exponential family model. The default is \code{TRUE}. This term can only be
    used with undirected networks.}

    \item{\code{gwdsp(alpha, fixed=FALSE)}}{\emph{Geometrically Weighted
    Dyadwise Shared Partner Distribution:} This term adds one network statistic
    to the model equal to the geometrically weighted dyadwise shared partner
    distribution with weight parameter \code{alpha}. The optional argument
    \code{fixed} is a logical indicating if the scale parameter \code{alpha} is
    to be fit as a curved exponential family model. The default is \code{TRUE}.
    This term can only be used with undirected networks.}

    \item{\code{gwesp(alpha, fixed=FALSE)}}{\emph{Geometrically Weighted
    Edgewise Shared Partner Distribution:} This term adds one network statistic
    to the model equal to the geometrically weighted edgewise shared partner
    distribution with weight parameter \code{alpha}. The optional argument
    \code{fixed} is a logical indicating if the scale parameter \code{alpha} is
    to be fit as a curved exponential family model. The default is \code{TRUE}.
    This term can only be used with undirected networks.}

%%   \item{\code{hamming(x, attrname)}}{\emph{Hamming Distance:}
%%   This term adds one statistic to the model equal to the Hamming distance
%%   of the network from the network specified by \code{x}.
%%   The optional argument \code{attrname} specifies an optional network 
%%   attribute of \code{x} to use instead.}
%
%    \item{\code{hamming(x)}}{\emph{Hamming Distance:}
%    This term adds one statistic to the model equal to the Hamming distance
%    of the network from the network specified by \code{x}.}

    \item{\code{idegree(d, attrname)}}{\emph{In-degree:} The \code{d} argument
    is a vector of distinct integers. This term adds one network statistic to
    the model for each element in \code{d}; the \eqn{i}th such statistic
    equals the number of nodes in the network of in-degree \code{d[i]}, i.e. the
    number of nodes with exactly \code{d[i]} in-edges. The optional argument
    \code{attrname} is a character string giving the name of an attribute in the
    network's vertex attribute list. If this is specified then the count is over
    the number of configuration where all nodes have the same value of the
    attribute. This term can only be used with directed networks; for undirected
    networks see \code{degree}. }

    \item{\code{isolates}}{\emph{Isolates:} This term adds one statistic to the
    model equal to the number of isolates in the network. For an undirected
    network, an isolate is defined to be any node with degree zero. This term
    can only be used with undirected networks. }

    \item{\code{istar(k, attrname)}}{\emph{In-stars: } The \code{k} argument is a
    vector of distinct integers. This term adds one network statistic to the
    model for each element in \code{k}. The \eqn{i}th such statistic counts the
    number of distinct \code{k[i]}-instars in the network, where a
    \eqn{k}-instar is defined to be a node \eqn{N} and a set of \eqn{k}
    different nodes \eqn{\{O_1, \dots, O_k\}} such that the ties
    \eqn{(O_i{\rightarrow}N)}{(O_i, N)} exist for \eqn{i=1, \dots, k}. The
    optional argument \code{attrname} is a character string giving the name of
    an attribute in the network's vertex attribute list. If this is specified
    then the count is over the number of \eqn{k}-instars where all nodes have
    the same value of the attribute. This term can only be used for directed
    networks; for undirected networks see \code{kstar}. Note that
    \code{istar(1)} is equal to \code{ostar(1)} and \code{edges}. }

    \item{\code{kstar(k, attrname)}}{\emph{k-Stars: } The \code{k} argument is a
    vector of distinct integers. This term adds one network statistic to the
    model for each element in \code{k}. The \eqn{i}th such statistic counts the
    number of distinct \code{k[i]}-stars in the network, where a \eqn{k}-star is
    defined to be a node \eqn{N} and a set of \eqn{k} different nodes
    \eqn{\{O_1, \dots, O_k\}} such that the ties \eqn{\{N, O_i\}} exist for
    \eqn{i=1, \dots, k}. The optional argument \code{attrname} is a character
    string giving the name of an attribute in the network's vertex attribute
    list. If this is specified then the count is over the number of
    \eqn{k}-stars where all nodes have the same value of the attribute. This
    term can only be used for undirected networks; for directed networks, see
    \code{istar}, \code{ostar}, and \code{m2star}. Note that \code{kstar(1)} is
    isomorphic to \code{edges}. }

    \item{\code{localtriangle(attrname)}}{\emph{Triangles Within Neighborhoods:
    } This term adds one statistic to the model equal to the number of triangles
    in the network between nodes close to each other. For an undirected network,
    a local triangle is defined to be any set of three edges between nodal pairs
    \eqn{\{(i,j), (j,k), (k,i)\}} that are in the same neighborhood. For a
    directed network, a triangle is defined as any set of three edges
    \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)} and either
    \eqn{(k{\rightarrow}i)} or \eqn{(k{\leftarrow}i)} where again all nodes are
    within the same neighborhood. The argument \code{attrname} is a network or
    an adjacency matrix that specifies whether the two nodes are in the same
    neighborhood. Note that this is technically a special case of
    \code{triangle}. }

    \item{\code{m2star}}{\emph{Mixed 2-stars, a.k.a 2-paths:} This term can only
    be specified with a directed network; for undirected networks see
    \code{kstar(2)}. This term adds one statistic to the model, equal to the
    number of mixed 2-stars in the network, defined as a pair of edges
    \eqn{(i{\rightarrow}j), (j{\rightarrow}k)}{(i,j), (j,k)}. A mixed 2-star is
    sometimes called a 2-path because it is a directed path of length 2 from
    \eqn{i} to \eqn{k} via \eqn{j}.}

    \item{\code{match(attrname, diff=FALSE)}}{\emph{Uniform Homophily
    and Differential Homophily:}  This is an alias for \code{nodematch(attrname,
    diff=FALSE)}.}

    \item{\code{meandeg}}{\emph{Mean Vertex Degree:}
    This term adds one network statistic to the model equal to the
    average degree of a vertex.  Note that this term is a constant multiple of
    both \code{edges} and \code{density}.}
  
    \item{\code{mutual}}{\emph{Mutuality:} This term adds one
    network statistic to the model, equaling the number of pairs
    of actors for which \eqn{\{(i{\rightarrow}j), (j{\rightarrow}i)\}} both exist.
    This term can only be implemented in a directed network.  }

    \item{\code{nodefactor(attrname)}}{\emph{Factor Attribute Effect.}  \code{attrname}
    is a character string giving the name of an attribute in
    the network's vertex attribute list.  This term adds \eqn{p-1}
    network statistics are added to the model, where \eqn{p} is the
    number of unique values of the \code{attrname} attribute.
    The \eqn{k^{th}} such statistic is the sum of
    \eqn{I\{\code{attrname(i)==value(k+1)}\}+
         I\{\code{attrname(i)==value(k+1)}\}}
    taken over all edges \eqn{(i,j)} in the network, where
    \code{value(k+1)} is the \eqn{(k+1)^{th}} smallest unique
    value of the \code{attrname} attribute and I is the
    indicator function.  }

    \item{\code{nodemain(attrname)}}{\emph{Main Effect of a Covariate.}
    \code{attrname} is a character string giving the name of an
    attribute in the network's vertex attribute list that takes
    numeric values.  This term adds one network statistic
    to the model equaling the sum of \code{attrname(i)}
    and \code{attrname(j)} for all \eqn{(i,j)} edges in
    the network. For categorical attributes see \code{nodefactor}.}

    \item{\code{nodematch(attrname, diff=FALSE)}}{\emph{Uniform Homophily
    and Differential Homophily.}  \code{attrname} is a
    character string giving the name of an attribute in the
    network's vertex attribute list.  When \code{diff=FALSE},
    this term adds one network statistic to the model,
    which counts the number of edges \eqn{(i,j)} for which
    \code{attrname(i)==attrname(j)}.  When \code{diff=TRUE},
    \eqn{p} network statistics are added to the model, where \eqn{p} is the
    number of unique values of the \code{attrname} attribute.
    The \eqn{k^{th}} such statistic counts the number of edges \eqn{(i,j)}
    for which \code{attrname(i) == attrname(j) == value(k)}, where
    \code{value(k)} is the \eqn{k^{th}} smallest unique value of the
    attrname attribute.  }

    \item{\code{nodemix(attrname, contrast=FALSE)}}{\emph{Nodal Attribute mixing.}
    \code{attrname} is a character string giving the name of an
    attribute in the network's vertex attribute list that takes
    categorical values.  This term adds one network statistic
    to the model equaling the count of the number of ties
    between nodes with each possible pairing of the attribute values;
    that is, a term for every entry in the mixing matrix for the attribute.
    The ordered of the attribute values is alphabetic.
    If the option \code{contrast=TRUE} then a statistic for 
    the first pairing is not included, making it a defacto reference category.}

    \item{\code{odegree(d, attrname)}}{\emph{Out-degree.}  \code{d}
    is a vector of distinct integers.  This term adds one
    network statistic to the model for each element in \code{d};
    the \eqn{i^{th}} such statistic equals the number of nodes
    in the network of out-degree \code{d[i]}, i.e. the number of
    nodes with exactly \code{d[i]} out-edges.
    The optional argument \code{attrname} is a character string giving the
    name of an attribute in the network's vertex attribute list.
    If this is specified then the count is over the number of configurations
    where all nodes have the same value of the attribute.
    This term can
    only be implemented with directed networks; for undirected
    networks see \code{degree}.
    }

    \item{\code{ostar(k, attrname)}}{\emph{k-outstars.}  \code{k}
    is a vector of distinct integers.  This term adds
    one network statistic to the model for each element in
    \code{k}.  The \eqn{i^{th}} such statistic counts the
    number of distinct \code{k[i]}-outstars in the network,
    where a k-outstar is defined to be a node N and a set of
    k different nodes \eqn{\{O_1, \dots, O_k\}} such that the
    ties \eqn{\{N{\rightarrow}O_i\}} exist for i=1, \dots, k. 
    The optional argument \code{attrname} is a character string giving the
    name of an attribute in the network's vertex attribute list.
    If this is specified then the count is over the number of k-outstars
    where all nodes have the same value of the attribute.
    This term
    can only be specified for directed networks; for undirected
    networks see \code{kstar}.  
    Note that \code{ostar(1)} is
    isomorphic to \code{istar(1)} and \code{edges}.
    }

    \item{\code{receiver}}{\emph{Receiver effect.}  This term adds one network 
    statistic for each node representing the number of in-ties for that node.
    This measures the popularity of the node.
    The term for the first node is omitted because of redundancy, but the 
    coefficient can be computed as the negative of the sum of the coefficients
    of all the other actors. That is the average coefficient is zero, following
    the Holland-Leinhardt parametrization of the \eqn{p_1} model.
    This term can only be specified with a directed network. For undirected
    networks, see \code{sociality}.}

    \item{\code{receivercov(attrname)}}{\emph{Receiver covariate effect.}
    \code{attrname} is a character string giving the name of an
    attribute in the network's vertex attribute list that takes
    numeric values.  This term adds one network statistic
    to the model equaling the sum of the attribute values of the receiver 
    of all in-ties.}

    \item{\code{sender}}{\emph{Sender effect.}  This term adds one network 
    statistic for each node representing the number of out-ties for that node.
    This measures the activity of the node.
    The term for the first node is omitted because if redundancy, but the 
    coefficient can be computed as the negative of the sum of the coefficients
    of all the other actors. That is the average coefficient is zero, following
    the Holland-Leinhardt parametrization of the \eqn{p_1} model.
    This term can only be specified with a directed network. For undirected
    networks, see \code{sociality}.}

    \item{\code{sendercov(attrname)}}{\emph{Sender covariate effect.}
    \code{attrname} is a character string giving the name of an
    attribute in the network's vertex attribute list that takes
    numeric values.  This term adds one network statistic
    to the model equaling the sum of the attribute values of the sender 
    of all out-ties.}

    \item{\code{smalldiff(attrname, cutoff)}}{\emph{Small
    difference.}  \code{attrname} is a character string giving
    the name of an attribute in the network's vertex attribute
    list and \code{cutoff} is any real number. This term adds
    one network statistic to the model, equal to the number of
    edges \eqn{(i,j)} for which \code{abs(attrname(i)-attrname(j))}
    is less than or equal to \code{cutoff}.  }

    \item{\code{sociality(attrname)}}{\emph{Centralized Covariate Effect:} This
    term adds one network statistic for each node equal to the number of ties of
    that node. The optional \code{attrname} is a character string giving the
    name of an attribute in the network's vertex attribute list that takes
    categorical values. If provided, this term only counts ties between nodes
    with the same value of the attribute. This term can only be used with
    undirected networks. For directed networks, see \code{sender} and
    \code{receiver}. }

    \item{\code{triangle(attrname)}}{\emph{Triangles.}  This term
    adds one statistic to the model equals to the number
    of triangles in the network.  For an undirected network,
    a triangle is defined to be any set of three edges
    between nodal pairs \eqn{\{(i,j), (j,k), (k,i)\}}.  For a
    directed network, a triangle is defined as any set of three edges
    \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k)\}} and either 
    \eqn{(k{\rightarrow}i)} or \eqn{(k{\leftarrow}i)}. Note
    that for all directed networks, triangle is isomorphic to
    \code{ttriad+ctriad}, so at most two of the three terms
    can be in a model.
    The optional argument \code{attrname} restricts the count to those
    triples of nodes with equal values of the vertex attribute specified
    by \code{attrname}.
    }

    \item{\code{tripercent(attrname)}}{\emph{Triangle Percentage.}  
    This term
    adds one statistic to the model equals to the percentage
    of triangles in the network relative to the number of potential triangles
    (i.e., 2-stars).  
    The optional argument \code{attrname} restricts the percent to those
    triples of nodes with equal values of the vertex attribute specified
    by \code{attrname}.
    For an undirected network, a triangle is defined to be any set of three edges
    between nodal pairs \eqn{\{(i,j), (j,k), (k,i)\}}.  For a
    directed network, a triangle is defined as any set of three edges
    \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k)\}} and either 
    \eqn{(k{\rightarrow}i)} or \eqn{(k{\leftarrow}i)}. The number
    of potential triangles is the number of out 2-stars.}

    \item{\code{ttriad(attrname)}}{\emph{Transitive Triads.}  
    This term adds one statistic to the model, equal to the
    number of transitive triads in the network, defined as
    a set of edges 
    \eqn{\{(i{\rightarrow}j), (j{\rightarrow}k), (i{\rightarrow}k)\}}.
     Note that for all directed networks, \code{triangle} is isomorphic to
    \code{ttriad+ctriad}, so at most two of the three terms
    can be in a model.
    The optional \code{attrname} is a character string giving the
    name of an attribute in the network's vertex attribute list.
    If this is specified then the count is over the number of cyclic triads where
    all three nodes have the same value of the attribute.
    This term can only be specified with a directed network.}
}

\references{
\item Hunter, D. R. and M. S. Handcock (2006), Inference in curved
exponential family models for networks, Journal of Computational
and Graphical Statistics, 15: 565-583.

\item Hunter, D. R. (2007), Curved exponential family models for social
networks, Social Networks, 29: 216-230.

\item Snijders, T. A. B., P. E. Pattison, G. L. Robins, and
M. S. Handcock (2006), New specifications for exponential random graph
models, Sociological Methodology, 36(1): 99-153. 
}
\seealso{ergm, network, \%v\%, \%n\%,
         summary.ergm, print.ergm}
\examples{
\dontrun{
ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)

ergm(molecule ~ edges + kstar(2:3) + triangle
                      + nodematch("atomic type",diff=TRUE)
                      + triangle + absdiff("atomic type"))
}
}
\keyword{models}
