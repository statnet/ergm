\name{ergm.formula.utils}
\alias{ergm.getmodel}
\alias{ergm.getnetwork}
\alias{ergm.getterms}
\alias{ergm.rhs.formula}
\alias{ergm.update.formula}
\alias{model.transform.formula}
\alias{set.offset.formula}
\alias{term.list.formula}
\alias{unset.offset.formula}
\alias{offset.info.formula}
\alias{remove.offset.formula}
\alias{append.rhs.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Functions for Manipulating and Extracting from ergm Formulas
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ergm.getmodel(formula, nw, response = NULL, silent = FALSE, role = "static", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{nw}{
%%     ~~Describe \code{nw} here~~
}
  \item{response}{
%%     ~~Describe \code{response} here~~
}
  \item{silent}{
%%     ~~Describe \code{silent} here~~
}
  \item{role}{
%%     ~~Describe \code{role} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, nw, response = NULL, silent = FALSE, role = "static", 
    ...) 
{
    if ((dc <- data.class(formula)) != "formula") 
        stop(paste("Invalid formula of class ", dc), call. = FALSE)
    if (formula[[1]] != "~") 
        stop("Formula must be of the form 'network ~ model'.", 
            call. = FALSE)
    if (length(formula) < 3) 
        stop(paste("No model specified for network ", formula[[2]]), 
            call. = FALSE)
    v <- term.list.formula(formula[[3]])
    formula.env <- environment(formula)
    model <- structure(list(formula = formula, coef.names = NULL, 
        offset = NULL, terms = NULL, networkstats.0 = NULL, etamap = NULL), 
        class = "model.ergm")
    termroot <- if (is.null(response)) 
        "InitErgm"
    else "InitWtErgm"
    for (i in 1:length(v)) {
        if (is.call(v[[i]])) {
            if (v[[i]][[1]] == "offset") {
                if (length(v[[i]][[2]]) <= 1) {
                  v[[i]] <- as.call(v[[i]][2])
                }
                else {
                  v[[i]] <- as.call(v[[i]][[2]])
                }
                model$offset <- c(model$offset, TRUE)
            }
            else {
                model$offset <- c(model$offset, FALSE)
            }
            args = v[[i]]
            args[[1]] = as.name("list")
            fname <- paste(termroot, "Term.", v[[i]][[1]], sep = "")
            newInitErgm <- exists(fname, envir = formula.env, 
                mode = "function")
            v[[i]] <- call(ifelse(newInitErgm, fname, paste(termroot, 
                ".", v[[i]][[1]], sep = "")))
        }
        else {
            fname <- paste(termroot, "Term.", v[[i]], sep = "")
            newInitErgm <- exists(fname, envir = formula.env, 
                mode = "function")
            v[[i]] <- call(ifelse(newInitErgm, fname, paste(termroot, 
                ".", v[[i]], sep = "")))
            model$offset <- c(model$offset, FALSE)
            args = list()
        }
        if (!newInitErgm) {
            v[[i]][[2]] <- nw
            names(v[[i]])[2] <- ""
            v[[i]][[3]] <- model
            names(v[[i]])[3] <- ""
            v[[i]][[4]] <- args
            dotdotdot <- c(if (!is.null(response)) list(response = response), 
                list(role = role), list(...))
            for (j in seq_along(dotdotdot)) {
                if (is.null(dotdotdot[[j]])) 
                  next
                v[[i]][[4 + j]] <- dotdotdot[[j]]
                names(v[[i]])[4 + j] <- names(dotdotdot)[j]
            }
            if (!exists(as.character(v[[i]][[1]]), envir = formula.env, 
                mode = "function")) {
                stop("The term ", substring(as.character(v[[i]][[1]]), 
                  first = nchar(termroot) + 2), " does not exist for this type of ERGM. Are you sure you have the right name?\n", 
                  call. = FALSE)
            }
            if (silent) {
                silentwarnings <- capture.output(model <- eval(v[[i]], 
                  formula.env))
            }
            else {
                model <- eval(v[[i]], formula.env)
            }
            if (is.null(model$terms[[length(model$terms)]]$pkgname)) 
                model$terms[[length(model$terms)]]$pkgname <- which.package.InitFunction(v[[i]][[1]], 
                  formula.env)
        }
        else {
            v[[i]][[2]] <- nw
            names(v[[i]])[2] <- ""
            v[[i]][[3]] <- args
            names(v[[i]])[3] <- ""
            dotdotdot <- c(if (!is.null(response)) list(response = response), 
                list(role = role), list(...))
            for (j in seq_along(dotdotdot)) {
                if (is.null(dotdotdot[[j]])) 
                  next
                v[[i]][[3 + j]] <- dotdotdot[[j]]
                names(v[[i]])[3 + j] <- names(dotdotdot)[j]
            }
            outlist <- eval(v[[i]], formula.env)
            if (is.null(outlist$pkgname)) 
                outlist$pkgname <- which.package.InitFunction(v[[i]][[1]], 
                  formula.env)
            model <- updatemodel.ErgmTerm(model, outlist)
        }
    }
    model$etamap <- ergm.etamap(model)
    ergm.MCMC.packagenames(unlist(sapply(model$terms, "[[", "pkgname")))
    class(model) <- "ergm.model"
    model
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
