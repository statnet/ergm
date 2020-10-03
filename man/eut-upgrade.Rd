%  File man/eut-upgrade.Rd in package ergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2003-2020 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{eut-upgrade}
\alias{eut-upgrade}

\title{Updating \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} prior to 3.1}

\description{Explanation and instructions for updating custom ERGM terms
  developed prior to the release of \code{\link[=ergm-package]{ergm}}
  version 3.1 (including 3.0--999 preview release) to be used with versions 3.1 or later.}

\section{Explanation}{

  \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} --- Statnet's
  mechanism enabling users to write their own ERGM terms --- comes in a
  form of an R package containing files for the user to put their own
  statistics into (i.e., \code{changestats.user.h},
  \code{changestats.user.c}, and \code{InitErgmTerm.user.R}), as well as
  some boilerplate to support them (e.g., \code{edgetree.h},
  \code{edgetree.c}, \code{changestat.h}, \code{changestat.c}, etc.).

  Although the \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} API
  is stable, recent developments in ergm have necessitated the
  boilerplate files in ergm.userterms to be updated.  To reiterate, the
  user-written statistic source code (\code{changestats.user.h},
  \code{changestats.user.c}, and \code{InitErgmTerm.user.R}) can be used
  without modification, but other files that came with the package need
  to be changed.

  To make things easier in the future, we have implemented a mechanism
  (using R's LinkingTo API, in case you are wondering) that will keep
  things in sync in releases after the upcoming one. However, for the
  upcoming release, we need to transition to this new mechanism.

}

\section{Instructions}{

  The transition entails the following steps. They only need to be done
  once for a package. Future releases will keep up to date
  automatically.

  \enumerate{
    
    \item Download the up-to-date \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} source from CRAN
    using, e.g., \code{\link{download.packages}} and unpack it.

    \item Copy the R and C files defining the user-written terms
    (usually \code{changestats.user.h}, \code{changestats.user.c}, and
    \code{InitErgmTerm.user.R}) and \emph{only} those files from the old
    \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}} source code to
    the new. Do \emph{not} copy the boilerplate files that you did not
    modify.
	
    \item If you have customized the package \code{DESCRIPTION} file
    (e.g., to change the package name) or \code{zzz.R} (e.g., to change
    the startup message), modify them as needed in the updated
    \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}}, but do
    \emph{not} simply overwrite them with their old versions.

    \item Make sure that your \code{\link[=ergm-package]{ergm}}
    installation is up to date, and rebuild
    \code{\link[ergm.userterms:ergm.userterms-package]{ergm.userterms}}.
    
  }
}
