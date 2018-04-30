
# Arguments and outputs are identical to the binary version, except for the C routine name.
InitWtErgmTerm..subnets <- function(...){
  modifyList(InitErgmTerm..subnets(...), list(name="_wtsubnets"))
}

  # Arguments and outputs are identical to the binary version, except for the C routine names.
InitWtErgmTerm.N <- function(...){
  term <- InitErgmTerm.N(...)
  term$name <- switch(term$name,
                      MultiNet = "wtMultiNet",
                      MultiNets = "wtMultiNets")
  term
}
