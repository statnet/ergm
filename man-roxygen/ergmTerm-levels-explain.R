#  File man-roxygen/ergmTerm-levels-explain.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################
#' @param levels this optional argument controls which levels of the attribute
#'   should be included and which should be excluded. (See Specifying Vertex
#'   attributes and Levels (`?nodal_attributes`) for details.) For example, if
#'   the "fruit" attribute has levels "orange", "apple", "banana", and "pear",
#'   then to include just two levels, one for "apple" and one for "pear", use
#'   any of `<%= fn %>("fruit", levels=-(2:3))`, `<%= fn %>("fruit", levels=c(1,4))`, and
#'   `<%= fn %>("fruit", levels=c("apple", "pear"))`. Note: if you are using numeric
#'   values to specify the levels of a character variable, the levels will
#'   correspond to the alphabetically sorted character levels.
