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
  \method{simulate}{stergm}(object, nsim=1, seed=NULL,
                            coef.form=object$formation.fit$coef,
                            coef.diss=object$dissolution.fit$coef,
                            monitor = NULL,
                            time.burnin=0, 
                            time.interval=1,
                            control=control.simulate.stergm(),
                            statsonly=FALSE,
                            toggles=!statsonly,
                            verbose=FALSE, 
                            \dots)
}
\arguments{
  \item{object}{an \R object of type \code{\link{stergm}}. 
  }
  \item{nsim}{Number networks to return from the dynamic process, by
    advancing it \code{time.burnin+time.interval*(nsim-1)} time steps forward.
  }
  \item{seed}{Random number integer seed.  
    See \code{\link[base]{set.seed}}.
  }
  \item{coef.form}{Parameters for the model from which the 
    post-formation network is drawn.
  }
  \item{coef.diss}{As \code{coef.form}, but for the post-dissolution network.
  }
  \item{monitor}{A one-sided formula specifying one or more terms whose
    value is to be monitored. These terms are appeneded to the model, along
    with a coefficient of 0, so their statistics are returned.
  }
  \item{time.burnin}{
    Number of time steps to discard before starting to collect statistics.
  }
  \item{time.interval}{
    Number of time steps between successive recordings of network statistics.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.simulate.stergm}}.
  }
  \item{statsonly}{Logical:  If TRUE, return only the network statistics, 
    not the network(s) themselves.
  }
  \item{toggles}{Logical:  If TRUE, return a list of time-stamped toggles.
    Only makes sense if \code{time.burnin==0} and \code{time.interval==1}.
  }
  \item{verbose}{Logical:  If TRUE, extra information is printed
    as the Markov chain progresses. 
  }
  \item{\dots}{Further arguments passed to or used by methods.
  }
}
\details{
  The dynamic process is run forward and the results are returned.
}
\value{
  If control parameter \code{final=TRUE}, returns the final
  network. Otherwise, returns an object of class \code{network.list}
  containing the initial network and the list of toggles.
}
