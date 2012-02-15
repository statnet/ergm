\name{simulate.stergm}
\alias{simulate.stergm}
\title{
Draw from the distribution of an Separable Temporal Exponential Family Random Graph Model
}
\description{
\code{\link[stats]{simulate}} is used to draw from separable temporal exponential family 
  random network models in their natural parameterizations.
  See \code{\link{stergm}} for more information on these models.
}
\usage{
  \method{simulate}{stergm}(object, nsim = 1, seed = NULL, 
      theta.form = object$coef.form, theta.diss = object$coef.diss, 
      time.burnin = 0, time.interval = 1, MH.burnin = object$MH.burnin, 
      constraints = object$constraints, stergm.order = object$stergm.order, 
      control = control.simulate.stergm(), toggles = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{object}{an \R object. Either
   a \code{\link{formula}} or an \code{\link{ergm}} object.
   The \code{\link{formula}} should be of the form \code{y ~ <model terms>},
   where \code{y} is a network object or a matrix that can be
   coerced to a \code{\link[network]{network}} object.  For the details 
   on the possible
   \code{<model terms>}, see \code{\link{ergm-terms}}.  To create a
   \code{\link[network]{network}} object in \R, use the 
   \code{network()} function,
   then add nodal attributes to it using the \code{\%v\%}
   operator if necessary.
}
  \item{nsim}{Number networks to return from the dynamic process, by
    advancing it \code{time.burnin+time.interval*(nsim-1)} time steps forward.
}
  \item{seed}{Random number integer seed.
   The default is \code{sample(10000000, size=1)}.
}
  \item{theta.form}{When given either
    a \code{\link{formula}} or an object of class \code{ergm}, 
    \code{theta.form} are the parameters
    from which the post-formation network is drawn.
}
\item{theta.diss}{
  When given either
    a \code{\link{formula}} or an object of class \code{ergm}, 
    \code{theta.diss} are the parameters
    from which the post-dissolution network is drawn.
}
\item{time.burnin}{
  Number of time steps to discard before starting to collect statistics.
}
  \item{time.interval}{
Number of time steps between successive recordings of network statistics.
}
  \item{MH.burnin}{
Number of Metropolis-Hastings steps per formation/dissolution ERGM draw.
}
  \item{constraints}{
A one-sided formula specifying one or more constraints
    on the support of the distribution of the networks being
    simulated. See the documentation for a similar argument for
    \code{\link{ergm}} for more information. For
    \code{simulate.formula.stergm}, defaults to no constraints. For
    \code{simulate.stergm}, defaults to using the same constraints as
    those with which \code{object} was fitted.
}
\item{stergm.order}{The "order" in which formation and dissolution phases occur:
    "FormThenDiss" and "DissThenForm" for sequential models, "FormAndDiss"
    for the simultaneous model (the default), and "FormOnly" and
    "DissOnly" for debugging.
}
  \item{control}{
A list of control parameters for algorithm
   tuning. Constructed using \code{\link{control.simulate.stergm}}.
}
\item{toggles}{
  Whether or not to return a list of time-stamped toggles. Defaults to
  \code{TRUE}, but only makes sense if \code{time.burnin==0} and \code{time.interval==1}.
}
  \item{verbose}{If this is \code{TRUE}, we will print out more information as 
  we run the program.
}
  \item{\dots}{Further arguments passed to or used by methods.
}
}
\details{
  The dynamic process is run forward and the results are returned.
}
\value{
  If control parameter \code{final=TRUE}, returns the final
  network. Otherwise, returns an object of class \code{network.series}
  containing the initial network and the list of toggles.
}
