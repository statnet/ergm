#  File man-roxygen/ergmTerm-sp-to-dsp.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
#' @note For directed networks, only outgoing two-path ("OTP") shared partners
#'   are counted.  In other words, for a <%= kind %> in a directed graph, the number of
#'   shared partners counted by `<%= fn %>` is the number of nodes `k` that have edges
#'   `i -> k -> j`.  (These may also be called homogeneous shared partners.)  To
#'   count other types of shared partners instead, see `<%= see %>`.
