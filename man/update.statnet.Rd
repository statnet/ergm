\name{update.statnet}
\alias{update.statnet}
\title{Updates or Installs the Component Packages of statnet}
\usage{
update.statnet(object, \dots,
            contriburl = "http://csde.washington.edu/statnet",
            repos = getOption("repos"), type = getOption("pkgType"))
}
\description{
 Updates or, if necessary, installs the component packages of 
\code{\link{statnet}}. Many of the packages are on \code{CRAN}. The rest are on
 the \code{\link{statnet}} webpage (\url{http://csde.washington.edu/statnet}).
 This command should be used immediately after installing \code{\link{statnet}}
 from \code{CRAN} to install the component packages.
 Consult the webpage for more information.
}
\arguments{
  \item{object}{an optional vector of the names of the packages to install. 
By default it updates or installs all the component packages of \code{\link{statnet}}.}
  \item{\dots}{Additional arguments (currently none).}
  \item{contriburl}{URL of the contrib section of the
    repositories. By default this is the URL of the \code{\link{statnet}} webpage
    (\url{http://csde.washington.edu/statnet}). Overrides argument \code{repos}.
#ifdef unix
    Can be \code{NULL} to install from local \file{.tar.gz} files.
#endif
#ifdef windows
    Can also be \code{NULL} to install from local \file{.zip} files.
#endif
  }
% \item{contriburl}{URL of the \code{\link{statnet}} webpage
%       (\url{http://csde.washington.edu/statnet})}
  \item{repos}{character vector, the base URL(s) of the repositories
    to use, i.e., the URL of the CRAN master such as
    \code{"http://cran.r-project.org"} or its Statlib mirror,
    \code{"http://lib.stat.cmu.edu/R/CRAN"}.
#ifdef unix
    Can be \code{NULL} to install from local \file{.tar.gz} files.
#endif
#ifdef windows
    Can be \code{NULL} to install from local \file{.zip} files.
#endif
  }
% \item{ask}{logical indicating whether to ask user before packages
%   are actually downloaded and installed, or the character string
%   \code{"graphics"}, which brings up a widget to allow the user to
%   (de-)select from the list of packages which could be updated.  The
%   latter only works on systems with a GUI version of
%   \code{\link{select.list}}, and is otherwise equivalent to \code{ask
%     = TRUE}.
% }
  \item{type}{character, indicating the type of \code{\link{statnet}}
    package to download and install.
    Possible values are \code{"base"} (for packages essential to \code{\link{statnet}}),
    \code{"recommended"} (for packages needed for a significant part of
     \code{\link{statnet}} functionality) and \code{"optional"}
    (for packages with specialized \code{\link{statnet}} functionality).
  }
}
\details{
  While this function has a number of options, few users will change the defaults
  and they are not intended to be altered except in rare circumstances.

  The main function is to install a functional set of
  \code{\link{statnet}} packages.  First a list
  of all packages found in the default library
  (the first directory in \code{.libPaths()})
% \code{lib.loc}
  is created and compared with those available at the \code{\link{statnet}}
  repositories (\code{CRAN} and \url{http://csde.washington.edu/statnet}).
  If \code{ask = TRUE} (the default for "recommended" and "optional", 
  but not for "base") packages with a newer version are
  reported and for each one the user can specify if it should be
  updated.
}
\value{
  \code{\link{update.statnet}} has no return value.
}
\section{Warning}{
  This function may fail in rare circumstances. If so, send email to
  \code{handcock@u.washington.edu} who will fix it.
}
\seealso{
  \code{\link{install.packages}}.
}
\examples{
\dontrun{
update.statnet()
}
}
\keyword{utilities}
}
