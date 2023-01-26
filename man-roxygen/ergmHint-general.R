#  File man-roxygen/ergmHint-general.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2023 Statnet Commons
################################################################################
#' <% name <- if(startsWith(name, "'")) substr(name, 2, 1000) else name %>
#' @name <%= name %>-ergmHint
#' @rdname <%= ergm:::.term.rdname("ergmHint", name) %>
#' @aliases <%= name %>-ergmConstraint
#' @seealso [`ergmHint`] for index of constraints and hints currently visible to the package.
#'
#' \Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmHint", "<%= name %>", "subsection")}
