\name{ergm.formula.utils}
\alias{ergm.getmodel}
\alias{ergm.getnetwork}
\alias{ergm.getterms}
%\alias{ergm.rhs.formula}
\alias{ergm.update.formula}
%\alias{model.transform.formula}
%\alias{set.offset.formula}
\alias{term.list.formula}
%\alias{unset.offset.formula}
\alias{offset.info.formula}
\alias{remove.offset.formula}
\alias{append.rhs.formula}
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

ergm.update.formula(object, new, ..., from.new = FALSE)

term.list.formula(rhs, sign = +1)

offset.info.formula(object, response = NULL) 

remove.offset.formula(object, response = NULL)

append.rhs.formula(object, newterms, keep.onesided = FALSE)
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
  \item{rhs}{a formula-style call containing the right hand side of
  formula, obtained by \code{fmla[[3]]} for a two-sided formula and
  \code{fmla[[2]]} for a one-sided formula.}
  \item{sign}{an internal parameter used by \code{term.list.formula}
    when calling itself recursively}
  \item{newterms}{list of terms (names) to append to the formula }
  \item{keep.onesided}{if the initial formula is one-sided, keep it
    whether to keep it one-sided or whether to make the initial formula
    the new LHS}
}
\details{
\itemize{
  \item{The \code{ergm.getmodel} function parses the given formula, and initiliazes each ergm
 term via the \code{InitErgmTerm} functions to create a \code{model.ergm} object for the
 given network}
 
 \item{The \code{ergm.getnetwork} function ensures that the network in a given formula
 is valid; if so, the network is returned; if not, execution is halted with warnings}
 
 \item{The \code{ergm.getterms} function returns the terms of a given formula and ensures that the formula is  indeed a formula with the necessary \code{~} operator}
 
 \item{\code{ergm.update.formula} is a reimplementation of \code{\link{update.formula}} that does not simplify.  Note that the resulting formula's environment is set as follows. If
 \code{from.new==FALSE}, it is set to that of object. Otherwise, a new
 sub-environment of object, containing, in addition, variables in new
 listed in from.new (if a character vector) or all of new (if TRUE).}
 
 \item{\code{term.list.formula} returns a list containing terms in a given
   formula, handling \code{+} and \code{-} operators and parentheses,
   and keeping track of whether a term has a plus or a minus sign. }
 
 \item{\code{offset.info.formula} returns the offset vectors associated with a formula.}
 
 \item{\code{remove.offset.formula} deletes all \code{\link{offset}} terms in an ERGM formula.}
 \item{\code{append.rhs.formula}  appends a list of terms to the RHS of a formula. If the formula is one-sided, the RHS becomes the LHS, if \code{keep.onesided==FALSE} (the default).}
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
\code{append.rhs.formula} each return an updated formula object
}

\examples{
## append.rhs.formula

append.rhs.formula(y~x,list(as.name("z1"),as.name("z2"))) # y~x+z1+z2
append.rhs.formula(~y,list(as.name("z"))) # y~z
append.rhs.formula(~y+x,list(as.name("z"))) # y+x~z
append.rhs.formula(~y,list(as.name("z")),TRUE) # ~y+z
}