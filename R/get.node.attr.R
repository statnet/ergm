#  File R/get.node.attr.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2025 Statnet Commons
################################################################################
###############################################################################
# The <get.node.attr> function returns the vector of covariates for the given
# network and specified attribute if the attribute exists - execution will
# halt if the attribute is not correctly given as a single string or is not 
# found in the vertex attribute list; optionally <get.node.attr> will also 
# check that return vector is numeric, halting execution if not
#
# --PARAMETERS--
#   nw          : a network object
#   attrname    : the name of a nodal attribute, as a character string
#   functionname: the name of the calling function; this is only used for
#                 the warning messages that accompany a halt
#   numeric     : whether to halt execution if the return vector is not
#                 numeric; default=FALSE
#   
# --RETURNED--
#   out:  the vector of 'attrname' covariates
#
###############################################################################



#' Retrieve and check assumptions about vertex attributes (nodal covariates) in
#' a network
#' 
#' The \code{get.node.attr} function returns the vector of nodal covariates for
#' the given network and specified attribute if the attribute exists -
#' execution will halt if the attribute is not correctly given as a single
#' string or is not found in the vertex attribute list; optionally
#' \code{get.node.attr} will also check that return vector is numeric, halting
#' execution if not. The purpose is to validate assumptions before passing
#' attribute data into an ergm term.
#' 
#' 
#' @param nw a [`network`] object
#' @param attrname the name of a nodal attribute, as a character string
#' @param functionname the name of the calling function a character string;
#' this is only used for the warning messages that accompany a halt
#' @param numeric logical, whether to halt execution if the return vector is
#' not numeric; default=FALSE
#' @return returns the vector of 'attrname' covariates for the vertices in the
#' network
#' @seealso [get.vertex.attribute()] for a version without
#' the checking functionality
#' @examples
#' 
#' data(faux.mesa.high)
#' get.node.attr(faux.mesa.high,'Grade')
#' @keywords internal
#' @export get.node.attr
get.node.attr <- function(nw, attrname, functionname=NULL, numeric=FALSE) {  
  ergm_get_vattr(attrname, nw, accept=if(numeric)"numeric"else"character")
}

#' @name nodal_attributes
#' @title Specifying nodal attributes and their levels
#'
#' @description This document describes the ways to specify nodal
#'   attributes or functions of nodal attributes and which levels for
#'   categorical factors to include. For the helper functions to
#'   facilitate this, see [`nodal_attributes-API`].
#'
#' @param object,l,a,n,into `COLLAPSE_SMALLEST`, `LARGEST`, and
#'   `SMALLEST` are technically functions but they are generally not
#'   called in a standard fashion but rather as a part of an vertex
#'   attribute specification or a level specification as described
#'   below. The above usage examples are needed to pass \R's package
#'   checking without warnings; please disregard them, and refer to the
#'   sections and examples below instead.
#'
#' @section Specifying nodal attributes:
#'
#' Term nodal attribute arguments, typically called `attr`, `attrs`, `by`, or
#' `on` are interpreted as follows: \describe{
#' 
#' \item{a character string}{Extract the vertex attribute with
#' this name.}
#' 
#' \item{a character vector of length > 1}{Extract the vertex
#' attributes and paste them together, separated by dots if the term
#' expects categorical attributes and (typically) combine into a
#' covariate matrix if it expects quantitative attributes.}
#' 
#' \item{a function}{The function is called on the LHS network and
#' additional arguments to [ergm_get_vattr()], expected to return a
#' vector or matrix of appropriate dimension. (Shorter vectors and
#' matrix columns will be recycled as needed.)}
#' 
#' \item{a formula}{The expression on the RHS of the formula is
#' evaluated in an environment of the vertex attributes of the
#' network, expected to return a vector or matrix of appropriate
#' dimension. (Shorter vectors and matrix columns will be recycled as
#' needed.) Within this expression, the network itself accessible as
#' either `.` or `.nw`. For example,
#' `nodecov(~abs(Grade-mean(Grade))/network.size(.))` would return the
#' absolute difference of each actor's "Grade" attribute from its
#' network-wide mean, divided by the network size.}
#'
#' \item{an `AsIs` object created by `I()`}{Use as is, checking only
#' for correct length and type.}
#' 
#' }
#'
#' Any of these arguments may also be wrapped in or piped through
#' `COLLAPSE_SMALLEST(attr, n, into)` or, `attr %>%
#' COLLAPSE_SMALLEST(n, into)`, a convenience function that will
#' transform the attribute by collapsing the smallest `n` categories
#' into one, naming it `into`. Note that `into` must be of the same
#' type (numeric, character, etc.) as the vertex attribute in
#' question. If there are ties for `n`th smallest category, they will
#' be broken in lexicographic order, and a warning will be issued.
#'
#' The name the nodal attribute receives in the statistic can be
#' overridden by setting a an [attr()]-style attribute `"name"`.
#'
#' @section Specifying categorical attribute levels and their ordering:
#'
#' For categorical attributes, to select which levels are of interest
#' and their ordering, use the argument `levels`.  Selection of nodes (from
#' the appropriate vector of nodal indices) is likewise handled as the
#' selection of levels, using the argument `nodes`.  These arguments are interpreted
#' as follows: \describe{
#'
#' \item{an expression wrapped in [I()]}{Use the given list of levels
#' as is.}
#' 
#' \item{a numeric or logical vector}{Used for indexing of a list of
#' all possible levels (typically, unique values of the attribute) in
#' default older (typically lexicographic), i.e.,
#' `sort(unique(attr))[levels]`. In particular, `levels=TRUE` will
#' retain all levels. Negative values exclude. Another special value
#' is `LARGEST`, which will refer to the most frequent category, so,
#' say, to set such a category as the baseline, pass
#' `levels=-LARGEST`. In addition, `LARGEST(n)` will refer to the `n`
#' largest categories. `SMALLEST` works analogously. If there are ties
#' in frequencies, they will be broken in lexicographic order, and a
#' warning will be issued. To specify numeric or logical levels
#' literally, wrap in [I()].}
#'
#'\item{[`NULL`]}{Retain all possible levels; usually equivalent to
#' passing `TRUE`.}
#'
#' \item{a character vector}{Use as is.}
#' 
#' \item{a function}{The function is called on the list of unique
#' values of the attribute, the values of the attribute themselves,
#' and the network itself, depending on its arity. Its return value is
#' interpreted as above.}
#'
#' \item{a formula}{The expression on the RHS of the formula is
#' evaluated in an environment in which the network itself is
#' accessible as `.nw`, the list of unique values of the attribute as
#' `.` or as `.levels`, and the attribute vector itself as
#' `.attr`. Its return value is interpreted as above.}
#'
#' \item{a matrix}{For mixing effects (i.e., `level2=` arguments), a
#' matrix can be used to select elements of the mixing matrix, either
#' by specifying a logical (`TRUE` and `FALSE`) matrix of the same
#' dimension as the mixing matrix to select the corresponding cells or
#' a two-column numeric matrix indicating giving the coordinates of
#' cells to be used.}
#' 
#' }
#' 
#' Note that `levels`, `nodes`, and others often have a default that is sensible for the
#' term in question.
#' 
#' @aliases attr attrname on by attrs node.attr nodal.attr vertex.attr node.attribute nodal.attribute vertex.attribute
#' @examples
#' library(magrittr) # for %>%
#'
#' data(faux.mesa.high)
#' 
#' # Activity by grade with a baseline grade excluded:
#' summary(faux.mesa.high~nodefactor(~Grade))
#' # Name overrides:
#' summary(faux.mesa.high~nodefactor("Form"~Grade)) # Only for terms that don't use the LHS.
#' summary(faux.mesa.high~nodefactor(~structure(Grade,name="Form")))
#' # Retain all levels:
#' summary(faux.mesa.high~nodefactor(~Grade, levels=TRUE)) # or levels=NULL
#' # Use the largest grade as baseline (also Grade 7):
#' summary(faux.mesa.high~nodefactor(~Grade, levels=-LARGEST))
#' # Activity by grade with no baseline smallest two grades (11 and
#' # 12) collapsed into a new category, labelled 0:
#' table(faux.mesa.high %v% "Grade")
#' summary(faux.mesa.high~nodefactor((~Grade) %>% COLLAPSE_SMALLEST(2, 0),
#'                                   levels=TRUE))
#'
#' # Handling of tied frequencies
#' faux.mesa.high %v% "Plans" <-
#'     sample(rep(c("College", "Trade School", "Apprenticeship", "Undecided"), c(80,80,20,25)))
#' summary(faux.mesa.high ~ nodefactor("Plans", levels = -LARGEST))
#' 
#' # Mixing between lower and upper grades:
#' summary(faux.mesa.high~mm(~Grade>=10))
#' # Mixing between grades 7 and 8 only:
#' summary(faux.mesa.high~mm("Grade", levels=I(c(7,8))))
#' # or
#' summary(faux.mesa.high~mm("Grade", levels=1:2))
#' # or using levels2 (see ? mm) to filter the combinations of levels,
#' summary(faux.mesa.high~mm("Grade",
#'         levels2=~sapply(.levels,
#'                         function(l)
#'                           l[[1]]%in%c(7,8) && l[[2]]%in%c(7,8))))
#'
#' # Here are some less complex ways to specify levels2. This is the
#' # full list of combinations of sexes in an undirected network:
#' summary(faux.mesa.high~mm("Sex", levels2=TRUE))
#' # Select only the second combination:
#' summary(faux.mesa.high~mm("Sex", levels2=2))
#' # Equivalently,
#' summary(faux.mesa.high~mm("Sex", levels2=-c(1,3)))
#' # or
#' summary(faux.mesa.high~mm("Sex", levels2=c(FALSE,TRUE,FALSE)))
#' # Select all *but* the second one:
#' summary(faux.mesa.high~mm("Sex", levels2=-2))
#' # Select via a mixing matrix: (Network is undirected and
#' # attributes are the same on both sides, so we can use either M or
#' # its transpose.)
#' (M <- matrix(c(FALSE,TRUE,FALSE,FALSE),2,2))
#' summary(faux.mesa.high~mm("Sex", levels2=M)+mm("Sex", levels2=t(M)))
#' # Select via an index of a cell:
#' idx <- cbind(1,2)
#' summary(faux.mesa.high~mm("Sex", levels2=idx))
#' # Or, select by specific attribute value combinations, though note
#' # the names 'row' and 'col' and the order for undirected networks:
#' summary(faux.mesa.high~mm("Sex",
#'                           levels2 = I(list(list(row="M",col="M"),
#'                                            list(row="M",col="F"),
#'                                            list(row="F",col="M")))))
#' # Note the warning: in an undirected network with identical row and
#' # column attributes, the mixing matrix is symmetric and only the
#' # upper triangle (where row < column) is valid, so the [M,F] cell
#' # will get a statistic of 0 with a warning.
#'
#' # mm() term allows two-sided attribute formulas with different attributes:
#' summary(faux.mesa.high~mm(Grade~Race, levels2=TRUE))
#' # It is possible to have collapsing functions in the formula; note
#' # the parentheses around "~Race": this is because a formula
#' # operator (~) has lower precedence than pipe (|>):
#' summary(faux.mesa.high~mm(Grade~(~Race) %>% COLLAPSE_SMALLEST(3,"BWO"), levels2=TRUE))
#'
#' # Some terms, such as nodecov(), accept matrices of nodal
#' # covariates. An certain R quirk means that columns whose
#' # expressions are not typical variable names have their names
#' # dropped and need to be adjusted. Consider, for example, the
#' # linear and quadratic effects of grade:
#' Grade <- faux.mesa.high %v% "Grade"
#' colnames(cbind(Grade, Grade^2)) # Second column name missing.
#' colnames(cbind(Grade, Grade2=Grade^2)) # Can be set manually,
#' colnames(cbind(Grade, `Grade^2`=Grade^2)) # even to non-variable-names.
#' colnames(cbind(Grade, Grade^2, deparse.level=2)) # Alternatively, deparse.level=2 forces naming.
#' rm(Grade)
#' \dontshow{
#' options(warn=1) # Print warnings immediately.
#' }
#' # Therefore, the nodal attribute names are set as follows:
#' summary(faux.mesa.high~nodecov(~cbind(Grade, Grade^2))) # column names dropped with a warning
#' summary(faux.mesa.high~nodecov(~cbind(Grade, Grade2=Grade^2))) # column names set manually
#' summary(faux.mesa.high~nodecov(~cbind(Grade, Grade^2, deparse.level=2))) # using deparse.level=2
#'
#' # Activity by grade with a random covariate. Note that setting an attribute "name" gives it a name:
#' randomcov <- structure(I(rbinom(network.size(faux.mesa.high),1,0.5)), name="random")
#' summary(faux.mesa.high~nodefactor(I(randomcov)))
NULL

#' @name nodal_attributes-API
#' @title Helper functions for specifying nodal attribute levels
#'
#' @description These functions are meant to be used in `InitErgmTerm` and other
#' implementations to provide the user with a way to extract nodal
#' attributes and select their levels in standardized and flexible
#' ways described under [`nodal_attributes`].
#'
#' @param object An argument specifying the nodal attribute to select
#'   or which levels to include.
#' @param nw Network on the LHS of the formula.
#' @param attr A vector of length equal to the number of nodes,
#'   specifying the attribute vector.
#' @param levels Starting set of levels to use; defaults to the sorted
#'   list of unique attributes.
#' @param bip Bipartedness mode: affects either length of attribute
#'   vector returned or the length permited: `"n"` for full network,
#'   `"b1"` for first mode of a bipartite network, `"b2"` for the
#'   second, and `"a"` for not adjusting.
#' @param multiple Handling of multiple attributes or matrix or data
#'   frame output. See the Details section for the specification.
#' @param accept A character vector listing permitted data types for
#'   the output. See the Details section for the specification.
#' @param a arguments to `LARGEST`, which is actually a function
#'   that gets processed as a function level spec does.
#' @param ... Additional argument to the functions of network or to
#'   the formula's environment.
#'
#' @details The `accept` argument is meant to allow the user to
#'   quickly check whether the output is of an *acceptable* class or
#'   mode. Typically, if a term accepts a character (i.e.,
#'   categorical) attribute, it will also accept a numeric one,
#'   treating each number as a category label. For this reason, the
#'   following outputs are defined:
#' \describe{
#'
#' \item{`"character"`}{Accept any mode or class (since it can
#' be converted to character).}
#' 
#' \item{`"numeric"`}{Accept real, integer, or logical.}
#' 
#' \item{`"logical"`}{Accept logical.}
#' 
#' \item{`"integer"`}{Accept integer or logical.}
#' 
#' \item{`"natural"`}{Accept a strictly positive integer.}
#' 
#' \item{`"0natural"`}{Accept a nonnegative integer or logical.}
#' 
#' \item{`"nonnegative"`}{Accept a nonnegative number or logical.}
#'
#' \item{`"positive"`}{Accept a strictly positive number or logical.}
#'
#' \item{`"index"`}{Accept input appropriate for selecting from an unnamed vector: an integer or a logical; positive integers are returned as they are (`bip` ignored), logicals are right-sized, and negative integers reverse the selection (as with vector indexing).}
#'
#' }
#'
#' Given that, the `multiple` argument controls how passing multiple
#' attributes or functions that result in vectors of appropriate
#' dimension are handled: \describe{
#' 
#' \item{`"paste"`}{Paste together with dot as the separator.}
#' 
#' \item{`"stop"`}{Fail with an error message.}
#'
#' \item{`"matrix"`}{Construct and/or return a matrix whose rows correspond to vertices.}
#'
#' }
#'
#'
NULL

#' @rdname nodal_attributes-API
#' @export
ERGM_GET_VATTR_MULTIPLE_TYPES <- c("paste", "matrix", "stop")

#' @rdname nodal_attributes-API
#'
#' @description `ergm_get_vattr` extracts and processes the specified
#'   nodal attribute vector. It is strongly recommended that
#'   [check.ErgmTerm()]'s corresponding
#'   `vartype="function,formula,character"` (using the
#'   `ERGM_VATTR_SPEC` constant).
#' 
#' @return `ergm_get_vattr` returns a vector of length equal to the
#'   number of nodes giving the selected attribute function or, if
#'   `multiple="matrix"`, a matrix whose number of row equals the
#'   number of nodes. Either may also have an attribute `"name"`, which
#'   controls the suggested name of the attribute combination.
#'
#' @examples
#' data(florentine)
#' ergm_get_vattr("priorates", flomarriage)
#' ergm_get_vattr(~priorates, flomarriage)
#' ergm_get_vattr(~cbind(priorates, priorates^2), flomarriage, multiple="matrix")
#' ergm_get_vattr(c("wealth","priorates"), flomarriage)
#' ergm_get_vattr(c("wealth","priorates"), flomarriage, multiple="matrix")
#' ergm_get_vattr(~priorates>30, flomarriage)
#' ergm_get_vattr(~TRUE, flomarriage, accept="index")
#' ergm_get_vattr(~-(2:12), flomarriage, accept="index")
#' (a <- ergm_get_vattr(~cut(priorates,c(-Inf,0,20,40,60,Inf),label=FALSE)-1, flomarriage))
#' @keywords internal
#' @export
ergm_get_vattr <- function(object, nw, accept="character", bip=c("n","b1","b2","a"), multiple=if(accept=="character") "paste" else "stop", ...){
  bip <- match.arg(bip)
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)
  UseMethod("ergm_get_vattr")
}

.handle_multiple <- function(a, multiple){
  name <- attr(a, "name")
  if(!is.list(a)) a <- list(a)
  a <- do.call(cbind, a)
  structure(
    if(ncol(a)>1)
      switch(multiple,
             paste =  apply(a, 1, paste, collapse="."),
             matrix = a,
             stop = ergm_Init_stop("This term does not accept multiple vertex attributes or matrix vertex attribute functions."))
    else c(a),
    name = name)
}

.rightsize_vattr <- function(a, nw, bip, accept){
  if(bip=="a") return(a)

  if(!is.bipartite(nw)) bip <- "n"

  name <- attr(a, "name")

  if(accept=="index" && is.numeric(a)){
    a <- a[a!=0]
    if(all(a>0)) return(a)

    if(!is.null(nrow(a))) ergm_Init_stop("Subtractive (negative) index matrices are not supported at this time.")
    # Now it's negative indices.
    l <- switch(bip,
                n=seq_len(network.size(nw)),
                b1=seq_len(nw%n%"bipartite"),
                b2=seq_len(network.size(nw)-nw%n%"bipartite")+nw%n%"bipartite")
    if(bip=="b2" && any(-a > nw%n%"bipartite")) a <- a + nw%n%"bipartite"
    if(!all(-a%in%seq_along(l))) ergm_Init_warning("Vertex index is out of bound.")

    structure(l[a], name=name)
  }else{
    rep_len_warn <- function(x, length.out){
      if(length.out%%NVL(nrow(x), length(x))) ergm_Init_warning("Network size or bipartite group size is not a multiple of the length of vertex attributes.")
      if(is.null(nrow(x))) rep_len(x, length.out) else apply(x, 2, rep_len, length.out)
    }

    structure(
      if(bip=="n") rep_len_warn(a, network.size(nw))
      else if(NVL(nrow(a), length(a))==network.size(nw)) # Input vector is n-long, need to trim.
        switch(bip,
               b1 = if(is.null(nrow(a))) a[seq_len(nw%n%"bipartite")] else a[seq_len(nw%n%"bipartite"),,drop=FALSE],
               b2 = if(is.null(nrow(a))) a[-seq_len(nw%n%"bipartite")] else a[-seq_len(nw%n%"bipartite"),,drop=FALSE])
      else # Othewise, recycle until the right length.
        rep_len_warn(a, switch(bip,
                               b1=nw%n%"bipartite",
                               b2=network.size(nw)-nw%n%"bipartite")),
      name = name)
  }
}

.check_acceptable <- function(x, accept=c("character", "numeric", "logical", "integer", "natural", "0natural", "nonnegative", "index"), xspec=NULL){
  accept <- match.arg(accept)

  ACCNAME <- list(character = "a character",
                  logical = "a logical",
                  numeric = "a numeric or logical",
                  integer = "an integer or logical",
                  natural = "a natural (positive integer) numeric",
                  `0natural` = "a nonnegative integer or logical",
                  nonnegative = "a nonnegative numeric or logical",
                  positive = "a positive numeric or logical",
                  index = "an integer (of the same sign) or logical")
  OK <-
    if(accept == "character") TRUE
    else if(!is.numeric(x) && !is.logical(x)) FALSE
    else switch(accept,
                numeric = TRUE,
                logical = all(x %in% c(FALSE, TRUE)),
                integer = all(round(x)==x),
                natural = all(round(x)==x) && x>0,
                `0natural` = all(round(x)==x) && x>=0,
                nonnegative = x>=0,
                positive = x>0,
                index = is.logical(x) || (all(round(x)==x) && (all(x>0) || all(x<0))))

  if(!OK) ergm_Init_stop("Attribute ", NVL3(xspec, paste0(sQuote(paste(deparse(.),collapse="\n")), " ")), "is not ", ACCNAME[[accept]], " vector as required.")
  if(any(is.na(x))) ergm_Init_stop("Attribute ", NVL3(xspec, paste0(sQuote(paste(deparse(.),collapse="\n")), " ")), "has missing data, which is not currently supported by ergm.")
  if(is.matrix(x) && !is.null(cn <- colnames(x))){
    if(any(cn=="")){
      ergm_Init_warning("Attribute specification ", NVL3(xspec, paste0(sQuote(paste(deparse(.),collapse="\n")), " ")), "is a matrix with some column names set and others not; you may need to set them manually. See example(nodal_attributes) for more information.")
      colnames(x) <- NULL
    }
  }
  x
}

#' @rdname nodal_attributes-API
#' @importFrom purrr map pmap_chr map_chr discard
#' @importFrom rlang set_names
#' @export
ergm_get_vattr.AsIs <- function(object, nw, accept="character", bip=c("n","b1","b2","a"), multiple=if(accept=="character") "paste" else "stop", ...){
  bip <- match.arg(bip)
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  object %>% .handle_multiple(multiple=multiple) %>%
    .rightsize_vattr(nw, bip, accept) %>% structure(name=attr(object,"name")) %>%
    structure(class = class(object) %>% discard(~(.=="AsIs")) %>%
    map_chr(identity)) %>% .check_acceptable(accept=accept, xspec=object)
}

#' @rdname nodal_attributes-API
#' @importFrom purrr map pmap_chr
#' @importFrom rlang set_names
#' @export
ergm_get_vattr.character <- function(object, nw, accept="character", bip=c("n","b1","b2","a"), multiple=if(accept=="character") "paste" else "stop", ...){
  bip <- match.arg(bip)
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  missing_attr <- setdiff(object, list.vertex.attributes(nw))
  if(length(missing_attr)){
    ergm_Init_stop(paste.and(sQuote(missing_attr)), " is/are not valid nodal attribute(s).")
  }

  object %>% map(~nw%v%.) %>% set_names(object) %>% .handle_multiple(multiple=multiple) %>%
    .rightsize_vattr(nw, bip, accept) %>% structure(name=paste(object, collapse=".")) %>%
    .check_acceptable(accept=accept, xspec=object)
}


#' @rdname nodal_attributes-API
#' @export
ergm_get_vattr.function <- function(object, nw, accept="character", bip=c("n","b1","b2","a"), multiple=if(accept=="character") "paste" else "stop", ...){
  bip <- match.arg(bip)
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  args <- list()
  for(aname in c("accept", "bip", "multiple"))
    if('...' %in% names(formals(object)) || aname %in% names(formals(object)))
      args[[aname]] <- get(aname)
  args <- c(list(nw), list(...), args)

  ergm_Init_try({
    a <- do.call(object, args)
    while(is(a,'formula')||is(a,'function')) a <- ergm_get_vattr(a, nw, accept=accept, bip=bip, multiple=multiple, ...)
    a %>% .rightsize_vattr(nw, bip, accept) %>% .handle_multiple(multiple=multiple) %>%
      structure(., name=NVL(attr(.,"name"), strtrim(despace(paste(deparse(body(object)),collapse="\n")),80)))
  }) %>% .check_acceptable(accept=accept)
}


#' @rdname nodal_attributes-API
#' @importFrom purrr map when
#' @importFrom tibble lst
#' @export
ergm_get_vattr.formula <- function(object, nw, accept="character", bip=c("n","b1","b2","a"), multiple=if(accept=="character") "paste" else "stop", ...){
  bip <- match.arg(bip)
  multiple <- match.arg(multiple, ERGM_GET_VATTR_MULTIPLE_TYPES)

  a <- list.vertex.attributes(nw)
  vlist <- c(a %>% map(~nw%v%.) %>% set_names(a),
             lst(`.`=nw, .nw=nw, ...))

  e <- ult(object)
  ergm_Init_try({
    a <- eval(e, envir=vlist, enclos=environment(object))
    while(is(a,'formula')||is(a,'function')) a <- ergm_get_vattr(a, nw, accept=accept, bip=bip, multiple=multiple, ...)
    a %>% .rightsize_vattr(nw, bip, accept) %>% .handle_multiple(multiple=multiple) %>%
      structure(., name=NVL(attr(.,"name"), if(length(object)>2) eval_lhs.formula(object) else despace(paste(deparse(e),collapse="\n"))))
  }) %>% .check_acceptable(accept=accept, xspec=object)
}

#' @rdname nodal_attributes-API
#'
#' @description `ergm_attr_levels` filters the levels of the
#'   attribute.  It is strongly recommended that [check.ErgmTerm()]'s
#'   corresponding
#'   `vartype="function,formula,character,numeric,logical,AsIs,NULL"` (using the
#'   `ERGM_LEVELS_SPEC` constant).
#' 
#' @return `ergm_attr_levels` returns a vector of levels to use and their order.
#' @examples
#' ergm_attr_levels(NULL, a, flomarriage)
#' ergm_attr_levels(-1, a, flomarriage)
#' ergm_attr_levels(1:2, a, flomarriage)
#' ergm_attr_levels(I(1:2), a, flomarriage)
#' @export
ergm_attr_levels <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  UseMethod("ergm_attr_levels")
}

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.numeric <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  levels[object]
}

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.logical <- ergm_attr_levels.numeric

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.AsIs <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  object
}

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.character <- ergm_attr_levels.AsIs

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.NULL <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  levels
}

#' @rdname nodal_attributes-API
#'
#' @note `ergm_attr_levels.matrix()` expects `levels=` to be a
#'   [`list`] with each element having length 2 and containing the
#'   values of the two categorical attributes being crossed. It also
#'   assumes that they are in the same order as the user would like
#'   them in the matrix.
#' @export
ergm_attr_levels.matrix <- function(object, attr, nw, levels=sort(unique(attr)), ...){

  # This should get the levels in the right order.
  ol <- levels %>% map(1L) %>% unique
  nol <- length(ol)
  il <- levels %>% map(2L) %>% unique
  nil <- length(il)

  # Construct a matrix indicating where on the levels list does each
  # element go. Then, indexing elements of m with either a logical
  # matrix or a two-column matrix of cell indices will produce a list
  # of level indices selected along with 0s, which can then be
  # dropped.
  ol2c <- match(levels%>%map(1L), ol)
  il2c <- match(levels%>%map(2L), il)
  m <- matrix(0L, nol, nil)
  m[cbind(ol2c,il2c)] <- seq_along(levels)

  sel <- switch(mode(object),
                logical = { # Binary matrix
                  if(any(dim(object)!=c(nol,nil))) ergm_Init_stop("Level combination selection binary matrix should have dimension ", nol, " by ", nil, " but has dimension ", nrow(object), " by ", ncol(object), ".") # Check dimension.
                  if(!is.directed(nw) && !is.bipartite(nw) && identical(ol,il)) object <- object | t(object) # Symmetrize, if appropriate.
                  object
                },
                numeric = { # Two-column index matrix
                  if(ncol(object)!=2) ergm_Init_stop("Level combination selection two-column index matrix should have two columns but has ", ncol(object), ".")
                  if(!is.directed(nw) && !is.bipartite(nw) && identical(ol,il)) object <- rbind(object, object[,2:1,drop=FALSE]) # Symmetrize, if appropriate.
                  object
                },
                ergm_Init_stop("Level combination selection matrix must be either numeric or logical.")
                )

  sel <- m[sel] %>% keep(`!=`,0L) %>% sort %>% unique
  levels[sel]
}

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.function <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  object <- if('...' %in% names(formals(object))) object(levels, attr, nw, ...)
            else switch(length(formals(object)),
                        object(levels),
                        object(levels, attr),
                        object(levels, attr, nw))
  ergm_attr_levels(object, attr, nw, levels, ...)
}

#' @rdname nodal_attributes-API
#' @export
ergm_attr_levels.formula <- function(object, attr, nw, levels=sort(unique(attr)), ...){
  vlist <- lst(`.`=levels, .levels=levels, .attr=attr, .nw=nw, ...)
  e <- ult(object)
  object <- eval(e, envir=vlist, enclos=environment(object))  
  ergm_attr_levels(object, attr, nw, levels, ...)
}

#' @rdname nodal_attributes-API
#' @export
ERGM_VATTR_SPEC <- "function,formula,character,AsIs"

#' @rdname nodal_attributes-API
#' @export
ERGM_VATTR_SPEC_NULL <- "function,formula,character,AsIs,NULL"

#' @rdname nodal_attributes-API
#' @export
ERGM_LEVELS_SPEC <- "function,formula,character,numeric,logical,AsIs,NULL,matrix"

rank_cut <- function(x, n, tie_action = c("warning", "error"), top = FALSE){    
  ordrank <- if(top) function(r) length(x) + 1 - r else identity
  s1 <- ordrank(rank(x, ties.method="min")) <= n
  s2 <- ordrank(rank(x, ties.method="max")) <= n

  if(identical(s1, s2)) which(s1)
  else{
    tie_action <- match.arg(tie_action)
    msg <- paste0("Levels ", paste.and(sQuote(names(x)[s1!=s2])), " are tied.")
    switch(tie_action,
           error = ergm_Init_stop(msg, " Specify explicitly."),
           warning = {
             ergm_Init_warning(msg, " Using the order given.")
             which(ordrank(rank(x, ties.method="first")) <= n)
           })
  }
}

levels_cut <- function(x, n, lvls = sort(unique(x)), top = FALSE, ...){
  f <- setNames(tabulate(match(x, lvls)), lvls)
  sel <- rank_cut(f, n, top=top, ...)
  if(missing(lvls)) lvls[sel] else sel
}

#' @rdname nodal_attributes
#' @export
LARGEST <- structure(function(l, a){
  if(!missing(a)) levels_cut(a, 1, l, top=TRUE) # passed as levels=LARGEST
  else{ # passed as levels=LARGEST(n): return a function
    n <- l
    structure(function(l, a){
      levels_cut(a, n, l, top=TRUE)
    }, class = c("ergm_levels_spec_function", "function"))
  }
}, class = c("ergm_levels_spec_function", "function"))

#' @rdname nodal_attributes
#' @export
SMALLEST <- structure(function(l, a){
  if(!missing(a)) levels_cut(a, 1, l) # passed as levels=SMALLEST
  else{ # passed as levels=SMALLEST(n): return a function
    n <- l
    structure(function(l, a){
      levels_cut(a, n, l)
    }, class = c("ergm_levels_spec_function", "function"))
  }
}, class = c("ergm_levels_spec_function", "function"))

#' @noRd
#' @export
`-.ergm_levels_spec_function` <- function(e1, e2){
  if(!missing(e2)) stop("Addition and subtraction of ERGM level specifications is not supported at this time.")

  structure(
    function(levels, attr, nw, ...)
      - if('...' %in% names(formals(e1))) e1(levels, attr, nw, ...)
        else switch(length(formals(e1)),
                    e1(levels),
                    e1(levels, attr),
                    e1(levels, attr, nw)),
    class = c("ergm_levels_spec_function", "function"))
}


#' @rdname nodal_attributes
#' @export
COLLAPSE_SMALLEST <- function(object, n, into){
  attr <- object
  function(...){
    vattr <- ergm_get_vattr(attr, ...)
    smallest <- levels_cut(vattr, n)
    vattr[vattr %in% smallest] <- into
    vattr
  }
}
