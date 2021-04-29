#' @param levels The optional `levels` argument controls which levels of the attribute
#'   should be included and which should be excluded. (See Specifying Vertex
#'   attributes and Levels (`?nodal_attributes`) for details.) For example, if
#'   the "fruit" attribute has levels "orange", "apple", "banana", and "pear",
#'   then to include just two levels, one for "apple" and one for "pear", use
#'   any of `#1("fruit", levels=-(2:3))`, `#1("fruit", levels=c(1,4))`, and
#'   `#1("fruit", levels=c("apple", "pear"))`. Note: if you are using numeric
#'   values to specify the levels of a character variable, the levels will
#'   correspond to the alphabetically sorted character levels.
