%  File ergm/man/simulate.stergm.Rd
%  Part of the statnet package, http://statnet.org
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) in
%    http://statnet.org/attribution
%
%  Copyright 2012 the statnet development team
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{simulate.stergm}
\alias{simulate.stergm}
\alias{simulate.network}
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
                             coef.form=object$formation.fit$coef,coef.diss=object$dissolution.fit$coef,
                             monitor = object$targets,
                             time.slices, time.burnin=0, time.interval=1,
                             control=control.simulate.stergm(),
                             statsonly=time.burnin>0||time.interval>1,
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             verbose=FALSE, 
                             \dots)
  \method{simulate}{network}(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form,coef.diss,
                             monitor = NULL,
                             time.slices, time.burnin=0, time.interval=1,
                             control=control.simulate.stergm(),
                             statsonly=time.burnin>0||time.interval>1,
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             verbose=FALSE,
                             \dots)

}
\arguments{
  \item{object}{an \R object of type \code{\link{stergm}} giving a model
    fit or of type \code{\link{network}} giving the initial network.
  }
  \item{formation, dissolution}{
    One-sided \code{\link{ergm}}-style formulas for the formation and
    dissolution models, respectively.
  }
  
  \item{nsim}{Number of replications (separate chains of networks) of
    the process to run and return.
  }
  \item{seed}{Random number integer seed.  
    See \code{\link[base]{set.seed}}.
  }
  
  \item{coef.form}{Parameters for the model from which the 
    post-formation network is drawn.
  }
  \item{coef.diss}{As \code{coef.form}, but for the post-dissolution network.
  }
  \item{monitor}{Either a one-sided formula specifying one or more terms whose
    value is to be monitored, or a string containing \code{"formation"} or
    \code{"dissolution"}, to monitor their respective terms, or
    \code{"all"} to monitor distinct terms from both.
  }
  \item{time.slices}{
    Number of networks or network statistics to return from each replication of the dynamic process. See below for
    return types.
    }
  \item{time.burnin}{
    Number of time steps to discard before starting to collect network
    statistics. Actual network will only be returned if \code{time.burnin==0}.
  }
  \item{time.interval}{
    Number of time steps between successive recordings of network
    statistics. Actual network will only be returned if \code{time.interval==1}.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.simulate.stergm}} or \code{\link{control.simulate.network}}.
  }
  \item{statsonly}{Logical:  If TRUE, return only the network statistics, 
    not the network(s) themselves. Defaults to FALSE if both
    \code{time.burnin==0} and \code{time.interval==1},  and TRUE otherwise.
  }
  \item{stats.form, stats.diss}{Logical: Whether to return
    formation/dissolution model statistics. This is not the recommended
    method: use \code{monitor} argument instead.}

  \item{verbose}{Logical:  If TRUE, extra information is printed
    as the Markov chain progresses. 
  }
  \item{\dots}{Further arguments passed to or used by methods.
  }
}
\details{
  The dynamic process is run forward and the results are returned.
}
\value{If \code{nsim>1}, 
  
  If \code{statsonly==TRUE}, then if \code{stats.form==FALSE} and
  \code{stats.diss==FALSE}, returns an \code{\link{mcmc}} matrix with
  monitored statistics, and if either of them is \code{TRUE}, returns a
  list containing elements \code{stats} for statistics specified in the
  \code{monitor} argument, and \code{stats.form} and \code{stats.diss}
  for the respective formation and dissolution statistics. Finally, if
  \code{nsim>1}, an \code{\link{mcmc.list}} (or list of them) of the
  statistics is returned instead.

  If \code{statsonly==FALSE}, a \code{\link[networkDynamic]{networkDynamic}} object
  representing the process, with ties present in the initial network
  having onset \code{-Inf} and ties present at the end of the simulation
  having terminus \code{+Inf}. Additionally, attributes
  (\code{\link{attr}}, not network attributes) are attached as follows:

  \item{formation, dissolution}{Formation and dissolution formulas
    used in the simulation,
    respectively.}
  \item{stats, stats.form, stats.diss}{Network statistics as above.}
  \item{coef.form, coef.diss}{Coefficients used in the simulation.}
  \item{start, end}{Start and end times of the simulation.}

  
  
}
