
# Arguments and outputs are identical to the binary version, except for the C routine name.
InitWtErgmTerm..subnets <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm..subnets
  #' @importFrom utils modifyList
  modifyList(f(...), list(name="_wtsubnets"))
}

# Arguments and outputs are identical to the binary version, except for the C routine names.
InitWtErgmTerm.N <- function(...){
  # Rename the function to avoid the extra nesting level in the
  # diagnostic messages.
  f <- InitErgmTerm.N
  term <- f(...)
  term$name <- switch(term$name,
                      MultiNet = "wtMultiNet",
                      MultiNets = "wtMultiNets")
  term
}
