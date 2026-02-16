#  File R/ergm.errors.R in package ergm, part of the Statnet suite of packages
#  for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2026 Statnet Commons
################################################################################
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
#' @param call.,call See [stop()] and [abort()] respectively; note the
#'   different defaults.
#'
#' @note At this time, the \CRANpkg{rlang} analogues
#'   `ergm_Init_stop()`, `ergm_Init_warning()`, and `ergm_Init_message()`
#'   all concatenate their arguments like their base \R
#'   counterparts. This may change in the future, and if you wish to
#'   retain their old behavior, please switch to their base \R
#'   analogues `ergm_Init_stop()`, `ergm_Init_warning()`, and
#'   `ergm_Init_message()`.
#' @importFrom rlang abort
#' @seealso [stop()], [abort()]
#' @name ergm-errors
#' @keywords internal
#' @export
ergm_Init_abort <- function(..., default.loc=NULL, call=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  abort(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...), call=call)
}

#' @rdname ergm-errors
#' @export
ergm_Init_stop <- function(..., call. = FALSE, default.loc=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  stop(paste0('In ', NVL(loc, default.loc, "unknown function"), ': '), ..., call.=call.)
}

#' @rdname ergm-errors
#' @seealso [warning()], [warn()]
#' @importFrom rlang warn
#' @export
ergm_Init_warn <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  warn(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...))
}

#' @rdname ergm-errors
#' @export
ergm_Init_warning <- function(..., call. = FALSE, default.loc=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  warning(paste0('In ', NVL(loc, default.loc, "unknown function"), ': '), ..., call.=call.)
}

#' @rdname ergm-errors
#' @seealso [message()], [inform()]
#' @importFrom rlang inform
#' @export
ergm_Init_inform <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  inform(paste0('In ', NVL(loc, default.loc, "unknown function"), ': ', ...))
}

#' @rdname ergm-errors
#' @export
ergm_Init_message <- function(..., default.loc=NULL){
  loc <- traceback.Initializers() %>% format_traceback()
  message(paste0('In ', NVL(loc, default.loc, "unknown function"), ': '), ...)
}

#' @describeIn ergm-errors A helper function that evaluates the
#'   specified expression in the caller's environment, passing any
#'   errors to [ergm_Init_stop()].
#' @param expr Expression to be evaluated (in the caller's
#'   environment).
#' @seealso [try()], [tryCatch()]
#' @export
ergm_Init_try <- function(expr){
  expr <- substitute(expr)
  tryCatch(eval(expr, parent.frame(1)),
           error = function(e) ergm_Init_stop(e$message))
}

format_traceback <- function(x){
  if(EVL(nrow(x)==0,TRUE)) return(NULL)
  x <- as.data.frame(x)[rev(seq_len(nrow(x))),,drop=FALSE]
  x <- paste0(ifelse(x$valued,"valued ", ""),
              x$type, " ",
              sQuote(x$name),
              ifelse(nchar(x$pkg),
                     paste0(' in package ', sQuote(x$pkg)), ""))

  if(length(x)==1) x
  else paste0(x[1], ' (', paste('called from', x[-1], collapse=', '), ')')
}

traceback.Initializers <- function(){
  pat <- '^((?<pkg>[^:]+):::?)?`?Init(?<valued>Wt)?Ergm(?<type>Term|Proposal|Reference|Constraint)\\.(?<name>[^`]*)`?$'
  traceback.search(pat, perl=TRUE) %>% map(regexpr_list, pat) %>% do.call(rbind,.)
}

# Search back in time through sys.calls() to find the name of the last
# function whose name matches a regex.
traceback.search <- function(pattern, ...) {
  sys.calls() %>%
    as.list() %>%
    map(~.[[1]]) %>%
    map_chr(deparse1) %>%
    grep(pattern, ., value=TRUE, ...)
}

regexpr_list <- function(x, pat){
  m <- attributes(regexpr(pat, x, perl=TRUE))
  Map(function(s, l, n){substr(x, s, s+l-1)}, c(m$capture.start)[-1], c(m$capture.length)[-1]) %>% setNames(m$capture.names[-1]) %>% within({valued<-(valued=="Wt"); type<-tolower(type)})
}
