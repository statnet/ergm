#  File R/ergm.errors.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Sensible error and warning messages by `ergm` initializers
#'
#' These functions use traceback and pattern matching to find which
#' `ergm` initializer caused the problem, and prepend this information
#' to the specified message. They are not meant to be used by
#' end-users, but may be useful to developers.
#' 
#' @param ... Objects that can be coerced (via [paste0()]) into a
#'   character vector, concatenated into the message.
#' @param default.loc Optional name for the source of the error, to be
#'   used if an initializer cannot be autodetected.
#'
#' @importFrom rlang abort
#' @seealso [stop()], [abort()]
#' @name ergm-errors
#' @keywords internal
#' @export
ergm_Init_abort <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format.traceback()
  abort(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...))
}

#' @rdname ergm-errors
#' @seealso [warning()], [warn()]
#' @importFrom rlang warn
#' @export
ergm_Init_warn <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format.traceback()
  warn(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...))
}

#' @rdname ergm-errors
#' @seealso [message()], [inform()]
#' @importFrom rlang inform
#' @export
ergm_Init_inform <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format.traceback()
  inform(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...))
}

format.traceback <- function(x){
  if(nrow(x)==0) return(NULL)
  x <- x[nrow(x):1,]
  x <- paste0(ifelse(x$valued,"valued ", ""),
              x$type, " ",
              sQuote(x$name),
              ifelse(nchar(x$pkg),
                     paste0(' in package ', sQuote(x$pkg)), ""))

  if(length(x)==1) x
  else paste0(x[1], ' (', paste('called from', x[-1], collapse=', '), ')')
}

traceback.Initializers <- function(){
  pat <- '^((?<pkg>[^:]+):::?)?Init(?<valued>Wt)?Ergm(?<type>Term|Proposal|Reference|Constraint)\\.(?<name>.*)$'
  traceback.search(pat, perl=TRUE) %>% map(regexpr_list, pat) %>% do.call(rbind,.)
}

# Search back in time through sys.calls() to find the name of the last
# function whose name matches a regex.
#' @importFrom purrr map_chr
traceback.search <- function(pattern, ...) {
  sys.calls() %>%
    as.list() %>%
    map(~.[[1]]) %>%
    map(deparse) %>%
    map_chr(paste, collapse="\n") %>%
    grep(pattern, ., value=TRUE, ...)
}

regexpr_list <- function(x, pat){
  m <- attributes(regexpr(pat, x, perl=TRUE))
  mapply(function(s, l, n){substr(x, s, s+l-1)}, c(m$capture.start)[-1], c(m$capture.length)[-1], SIMPLIFY=FALSE) %>% set_names(m$capture.names[-1]) %>% within({valued<-(valued=="Wt"); type<-tolower(type)})
}
