#  File man-roxygen/ergmProposal-general.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' <% name <- if(startsWith(name, "'")) substr(name, 2, 1000) else name %>
#' @name <%= name %>-ergmProposal
#' @rdname <%= ergm:::.term.rdname("ergmProposal", name) %>
#' @details
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsHtml(ergm:::.buildProposalsList(proposal="<%= name %>"))}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsText(ergm:::.buildProposalsList(proposal="<%= name %>"))}}
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatProposalsLatex(ergm:::.buildProposalsList(proposal="<%= name %>"))}}
#' @keywords internal
#' @seealso [`ergmProposal`] for index of proposals currently visible to the package.
#'
#' \Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmProposal", "<%= name %>", "subsection")}
