%  File man/macros/uid-algo.Rd in package ergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2003-2023 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\UIDalgo}{concatenating the UNIX time (\code{\link{Sys.time}()}) to maximum available precision, process ID (\code{\link{Sys.getpid}()}), and a counter that starts at \code{-\link{.Machine}$integer.max} and increments by 1 with every call}