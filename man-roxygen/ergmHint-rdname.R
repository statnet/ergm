#  File man-roxygen/ergmHint-rdname.R in package ergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
#' <% name <- if(startsWith(name, "'")) substr(name, 2, 1000) else name %>
#' @rdname <%= ergm:::.term.rdname("ergmHint", name) %>
