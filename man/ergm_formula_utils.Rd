%  File man/ergm_formula_utils.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm.formula.utils}
\alias{ergm.getmodel}
\alias{ergm.getnetwork}
\alias{ergm.getterms}
%\alias{ergm.rhs.formula}
\alias{ergm.update.formula}
%\alias{model.transform.formula}
%\alias{set.offset.formula}
%\alias{unset.offset.formula}
\alias{offset.info.formula}
\alias{remove.offset.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Functions for Querying, Validating and Extracting from ERGM Formulas
}
\description{
These are all functions that are generally not called directly by users, but may be employed by other depending packages.  
}
\usage{
ergm.getmodel(formula, nw, response = NULL, silent = FALSE, role = "static", ...)

ergm.getnetwork(form, loopswarning = TRUE)

ergm.getterms(formula)

offset.info.formula(object, response = NULL) 

remove.offset.formula(object, response = NULL)

# Deprecated. Use nonsimp.update.formula() from statnet.common package.
ergm.update.formula(object, new, ..., from.new = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula of the form \code{network ~ model.term(s)}
}
  \item{nw}{
the network of interest
}
  \item{response}{
charcter, name of edge attribute containing edge weights
}
  \item{silent}{
logical, whether to print the warning messages from the initialization of each model term; default=FALSE
}
  \item{role}{ A hint about how the model will be used. Used primarily for dynamic
    network models. }
  \item{\dots}{
additional parameters for model formulation
}
  \item{form}{same as formula, a formula of the form \code{'network ~ model.term(s)'}}
  \item{loopswarning}{whether warnings about loops should be printed (T or F);default=TRUE}
  
  \item{object}{formula object to be updated}
  \item{new}{new formula to be used in updating}
  \item{from.new}{logical or character vector of variable names. controls how environment of formula gets updated.}
}
\details{
\itemize{
  \item{The \code{ergm.getmodel} function parses the given formula, and initiliazes each ergm
 term via the \code{InitErgmTerm} functions to create a \code{model.ergm} object for the
 given network}
 
 \item{The \code{ergm.getnetwork} function ensures that the network in a given formula
 is valid; if so, the network is returned; if not, execution is halted with warnings}
 
 \item{The \code{ergm.getterms} function returns the terms of a given formula and ensures that the formula is  indeed a formula with the necessary \code{~} operator}
 
 \item{\code{ergm.update.formula} (DEPRECATED: use
 \code{\link[statnet.common]{nonsimp.update.formula}} instead) is a reimplementation of \code{\link{update.formula}} that does not simplify.  Note that the resulting formula's environment is set as follows. If
 \code{from.new==FALSE}, it is set to that of object. Otherwise, a new
 sub-environment of object, containing, in addition, variables in new
 listed in from.new (if a character vector) or all of new (if TRUE).}
  
 \item{\code{offset.info.formula} returns the offset vectors associated with a formula.}
 
 \item{\code{remove.offset.formula} deletes all \code{\link{offset}} terms in an ERGM formula.}
}
}
\value{
\code{ergm.getmodel} returns a 'model.ergm' object as a list containing:
  \itemize{
    \item{ formula       :  the formula inputted to \code{\link{ergm.getmodel}}}
    \item{ coef.names    :  a vector of coefficient names}
    \item{ offset        :  a logical vector of whether each term was "offset", i.e. fixed}
    \item{ terms         :  a list of terms and 'term components' initialized by the 
                          appropriate \code{InitErgmTerm.X} function. } 
    \item{ network.stats0:  NULL always??}
    \item{ etamap        :  the theta -> eta mapping as a list returned from <ergm.etamap> }
    \item{ class         :  the character string "model.ergm"} 
  }
\code{ergm.getnetwork} returns:
  \itemize{
    \item{the network from the formula IF (i) the formula was correctly structured
      and (ii) the network is found within the formula's enviornment}  
  }
  \code{ergm.getterms} returns:
  \itemize{
    \item{the terms object associated with the formula and returned by the 
         native R function \code{\link{terms}}>. see \code{\link[stats]{terms.object}} for details about
         the components}
     }
     \code{terms.list.formula} returns a list of formula terms, each of
     witch having an additional attribute \code{"sign"}.
     
  \code{ergm.update.formula}, \code{remove.offset.formula} and 
}

