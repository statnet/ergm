%  File man/plot.network.ergm.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{plot.network.ergm}
\alias{plot.network.ergm}
\title{Two-Dimensional Visualization of Networks }
\description{
  \code{\link{plot.network.ergm}} produces a simple two-dimensional plot of
  the network object \code{x}. A variety of options are available
  to control vertex placement, display details, color, etc.
  The function is based on the 
  plotting capabilities of the \code{\link[network]{network}} package with additional
  pre-processing of arguments.
  Some of the capabilites require the \code{\link[latentnet]{latentnet}} package.
  See \code{\link[network]{plot.network}} in the \code{\link[network]{network}} package for details.
}
\usage{
\method{plot.network}{ergm}(x,
    attrname=NULL,
    label=network.vertex.names(x),
    coord=NULL,
    jitter=TRUE,
    thresh=0,
    usearrows=TRUE,
    mode="fruchtermanreingold",
    displayisolates=TRUE,
    interactive=FALSE,
    xlab=NULL,                                      
    ylab=NULL,
    xlim=NULL,
    ylim=NULL,
    pad=0.2,
    label.pad=0.5,
    displaylabels=FALSE,
    boxed.labels=TRUE,
    label.pos=0,
    label.bg="white",
    vertex.sides=8,
    vertex.rot=0,
    arrowhead.cex=1,
    label.cex=1,
    loop.cex=1,
    vertex.cex=1,
    edge.col=1,
    label.col=1,
    vertex.col=2,
    label.border=1,
    vertex.border=1,
    edge.lty=1,
    label.lty=NULL,
    vertex.lty=1,
    edge.lwd=0,
    label.lwd=par("lwd"),
    edge.len=0.5,
    edge.curve=0.1,
    edge.steps=50,
    loop.steps=20,
    object.scale=0.01,
    uselen=FALSE,
    usecurve=FALSE,
    suppress.axes=TRUE,
    vertices.last=TRUE,
    new=TRUE,
    layout.par=NULL,
    cex.main=par("cex.main"), 
    cex.sub=par("cex.sub"),
    seed=NULL,
    latent.control=list(maxit=500,
                        trace=0,
                        dyadsample=10000,
                        penalty.sigma=c(5,0.5), 
                        nsubsample=200),
    colornames="rainbow",
    verbose=FALSE, 
    latent=FALSE, 
    \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{\link[network]{network}}. }
  \item{attrname}{ an optional edge attribute, to be used to set edge values. }
  \item{label}{ a vector of vertex labels, if desired; defaults to the vertex labels returned by \code{\link[network]{network.vertex.names}}. }
  \item{coord}{ user-specified vertex coordinates, in an NCOL(dat)x2 matrix.  Where this is specified, it will override the \code{mode} setting. }
  \item{jitter}{ boolean; should the output be jittered? }
  \item{thresh}{ real number indicating the lower threshold for tie values.  Only ties of value >\code{thresh} are displayed.  By default, \code{thresh}=0.}
  \item{usearrows}{ boolean; should arrows (rather than line segments) be used to indicate edges? }
  \item{mode}{ the vertex placement algorithm; this must correspond to a \code{network.layout} function.
    These include 
    \code{"latent"}, \code{"latentPrior"},
     and
    \code{"fruchtermanreingold"}. }
  \item{displayisolates}{ boolean; should isolates be displayed? }
  \item{interactive}{ boolean; should interactive adjustment of vertex placement be attempted? }
  \item{xlab}{ x axis label. }
  \item{ylab}{ y axis label. }
  \item{xlim}{ the x limits (min, max) of the plot. }
  \item{ylim}{ the y limits of the plot. } 
  \item{pad}{ amount to pad the plotting range; useful if labels are being clipped. }
  \item{label.pad}{ amount to pad label boxes (if \code{boxed.labels==TRUE}), in character size units. }
  \item{displaylabels}{ boolean; should vertex labels be displayed? }
  \item{boxed.labels}{ boolean; place vertex labels within boxes? }
  \item{label.pos}{ position at which labels should be placed, relative to vertices.  \code{0} results in labels which are placed away from the center of the plotting region; \code{1}, \code{2}, \code{3}, and \code{4} result in labels being placed below, to the left of, above, and to the right of vertices (respectively); and \code{label.pos>=5} results in labels which are plotted with no offset (i.e., at the vertex positions). }
  \item{label.bg}{ background color for label boxes (if \code{boxed.labels==TRUE}); may be a vector, if boxes are to be of different colors.}
  \item{vertex.sides}{ number of polygon sides for vertices; may be given as a vector or a vertex attribute name, if vertices are to be of different types. }
  \item{vertex.rot}{ angle of rotation for vertices (in degrees); may be given as a vector or a vertex attribute name, if vertices are to be rotated differently. }
  \item{arrowhead.cex}{ expansion factor for edge arrowheads.}
  \item{label.cex}{ character expansion factor for label text. }
  \item{loop.cex}{ expansion factor for loops; may be given as a vector or a vertex attribute name, if loops are to be of different sizes. }
  \item{vertex.cex}{ expansion factor for vertices; may be given as a vector or a vertex attribute name, if vertices are to be of different sizes. }
  \item{edge.col}{ color for edges; may be given as a vector, adjacency matrix, or edge attribute name, if edges are to be of different colors. }
  \item{label.col}{ color for vertex labels; may be given as a vector or a vertex attribute name, if labels are to be of different colors. }
  \item{vertex.col}{ color for vertices; may be given as a vector or a vertex attribute name, if vertices are to be of different colors. }
  \item{label.border}{ label border colors (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different colors. }
  \item{vertex.border}{ border color for vertices; may be given as a vector or a vertex attribute name, if vertex borders are to be of different colors. }
  \item{edge.lty}{ line type for edge borders; may be given as a vector, adjacency matrix, or edge attribute name, if edge borders are to have different line types. }
  \item{label.lty}{ line type for label boxes (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different line types. }
  \item{vertex.lty}{ line type for vertex borders; may be given as a vector or a vertex attribute name, if vertex borders are to have different line types. }
  \item{edge.lwd}{ line width scale for edges; if set greater than 0, edge widths are scaled by \code{edge.lwd*dat}.  May be given as a vector, adjacency matrix, or edge attribute name, if edges are to have different line widths. }
  \item{label.lwd}{ line width for label boxes (if \code{boxed.labels==TRUE}); may be given as a vector, if label boxes are to have different line widths. }
  \item{edge.len}{ if \code{uselen==TRUE}, curved edge lengths are scaled by \code{edge.len}. }
  \item{edge.curve}{ if \code{usecurve==TRUE}, the extent of edge curvature is controlled by \code{edge.curv}.  May be given as a fixed value, vector, adjacency matrix, or edge attribute name, if edges are to have different levels of curvature.  }
  \item{edge.steps}{ for curved edges (excluding loops), the number of line segments to use for the curve approximation. }
  \item{loop.steps}{ for loops, the number of line segments to use for the curve approximation. }
  \item{object.scale}{ base length for plotting objects, as a fraction of the linear scale of the plotting region. Defaults to 0.01.}
  \item{uselen}{ boolean; should we use \code{edge.len} to rescale edge lengths? }
  \item{usecurve}{ boolean; should we use \code{edge.curve}? }
  \item{suppress.axes}{ boolean; suppress plotting of axes? }
  \item{vertices.last}{ boolean; plot vertices after plotting edges? }
  \item{new}{ boolean; create a new plot?  If \code{new==FALSE}, vertices and edges will be added to the existing plot. }
  \item{layout.par}{ parameters to the \code{network.layout} function specified in \code{mode}. }
  \item{cex.main}{ Character expansion for the plot title.}
  \item{cex.sub}{ Character expansion for the plot sub-title.}
  \item{seed}{ Integer for seeding random number generator.  See \code{\link[base]{set.seed}}.}
  \item{latent.control}{A list of parameters to control the \code{latent} and \code{latentPrior} models, \code{dyadsample} determines the size above which to sample the latent dyads; see \code{\link{ergm}} and \code{\link[stats]{optim}} for details.}
  \item{colornames}{A vector of color names that can be selected by index for
the plot. By default it is \code{colors()}.}
  \item{verbose}{logical; if this is \code{TRUE}, we will print out more information as we run the function.}
  \item{latent}{logical; use a two-dimensional latent space model based on the MLE fit.  See documentation for \code{ergmm()} in \code{\link[latentnet]{latentnet}}.}
  \item{\dots}{ additional arguments to \code{\link{plot}}. }
}
\details{
\code{\link[network]{plot.network}} is a version of the standard network visualization tool
within the \code{sna} package.  By means of clever selection
of display parameters, a fair amount of display flexibility
can be obtained.  Network layout -- if not specified directly
using \code{coord} -- is determined via one of the various
available algorithms.  These are (briefly) as follows:
\enumerate{

\item \code{latentPrior}: Use a two-dimensional latent space model based on a
Bayesian minimum Kullback-Leibler fit. 
See documentation for \code{latent()} in \code{\link{ergm}}.

\item \code{random}: Vertices are placed (uniformly) randomly within a square region about the origin.

\item \code{circle}: Vertices are placed evenly about the unit circle.

\item \code{circrand}: Vertices are placed in a ``Gaussian donut,'' with distance from the origin following a normal distribution and angle relative to the X axis chosen (uniformly) randomly.

\item \code{eigen}, \code{princoord}: Vertices are placed via (the real components of) the first two eigenvectors of: 
  \enumerate{
  \item \code{eigen}: the matrix of correlations among (concatenated) rows/columns of the adjacency matrix

  \item \code{princoord}: the raw adjacency matrix.
  }

\item \code{mds}, \code{rmds}, \code{geodist}, \code{adj}, \code{seham}: Vertices are placed by a metric MDS.  The distance matrix used is given by:
  \enumerate{
  \item \code{mds}: absolute row/column differences within the adjacency matrix

  \item \code{rmds}: Euclidean distances between rows of the adjacency matrix
  
  \item \code{geodist}: geodesic distances between vertices within the network
  
  \item \code{adj}: \eqn{(\max A)-A}{(max A)-A}, where \eqn{A}{A} is the raw adjacency matrix
  
  \item \code{seham}: structural (dis)equivalence distances 
     (i.e., as per \code{sedist} in the package \code{sna})
     based on the Hamming metric
  }
  
\item \code{spring}, \code{springrepulse}: Vertices are placed
      using a simple spring embedder.  Parameters for the embedding
      model are given by \code{embedder.params}, in the following
      order: vertex mass; equilibrium extension; spring coefficient;
      repulsion equilibrium distance; and base coefficient of
      friction.  Initial vertex positions are in random order around
      a circle, and simulation proceeds -- increasing the coefficient
      of friction by the specified base value per unit time -- until
      ``motion'' within the system ceases.  If \code{springrepulse}
      is specified, then an inverse-cube repulsion force between
      vertices is also simulated; this force is calibrated so as to
      be exactly equal to the force of a unit spring extension at
      a distance specified by the repulsion equilibrium distance.  }
}
\value{
None.
}
\references{ Wasserman, S., and Faust, K.  (1994).  ``Social Network Analysis: Methods and Applications.''  Cambridge: Cambridge University Press. }
\author{Carter T. Butts \email{buttsc@uci.edu} }
%\note{  }

 %~Make other sections like WARNING with \section{WARNING }{....} ~
\section{Requires}{\code{mva}}

\seealso{\code{\link{plot}} }

\examples{
data(florentine)
plot(flomarriage)  #Plot the Florentine Marriage data
plot(network(10))  #Plot a random network
\dontrun{plot(flomarriage,interactive="points")}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
