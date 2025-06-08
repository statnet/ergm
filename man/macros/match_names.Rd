%  File man/macros/match_names.Rd in package ergm, part of the Statnet suite of
%  packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free, open
%  source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2025 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\matchnames}{If the vector is named, its names will be matched to the corresponding #1 names, and if the named vector has a single coefficient without a name, it will be used for the unmatched #1s. In particular, \code{\link{setNames}(x,"")} will be treated as a vector of \code{x}s.}
