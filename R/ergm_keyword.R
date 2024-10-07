#  File R/ergm_keyword.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2024 Statnet Commons
################################################################################

DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS <- list('name'=20, 'short'=8, 'description'=30, 'popular'=8, 'package'=8)
DISPLAY_LATEX_KW_INDEX_PCT_WIDTHS <- c(0.15, 0.10, 0.6, 0.05, 0.10)

#' Dynamic ERGM keyword registry
#'
#' A function to manage dynamic ERGM keywords. To register a keyword, call the function with all parameters
#' provided. To fetch all registered keywords, call the function with no parameters specified.
#'
#' @param name full name of the keyword
#' @param short abbreviation of the keyword name
#' @param description description of the keyword
#' @param popular logical to indicate if a keyword is popular
#' @param package package the keyword is first defined in
#' @return Returns a dataframe with the following columns:
#'   - name
#'   - short
#'   - description
#'   - popular
#'   - package
#' @keywords models internal
#' @export
ergm_keyword <- local({
  cache <- data.frame(name=c(), short=c(), description=c(), popular=c(), package=c())

  function(name=NULL, short=NULL, description=NULL, popular=NULL, package=NULL) {
    if (all(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      return(cache)
    } else if (any(is.null(name), is.null(short), is.null(description), is.null(popular), is.null(package))) {
      stop("All arguments are needed to register ", sQuote("ergm"), " keyword")
    } else if (!is.logical(popular)) {
      stop("Logical value expected for argument 'popular'")
    } else {
      cache <<- rbind.data.frame(cache, data.frame(name=name, short=short, description=description, popular=popular, package=package, stringsAsFactors=FALSE))
    }
  }
})

.formatKeywordsText <- function(df) {
  line_wrap <- function(lines, max_width) {
    lines <- unlist(strsplit(sapply(strsplit(lines, '\n'), stringr::str_wrap, max_width), '\n'))

    out <- c()
    for (line in lines) {
      while (nchar(line) > max_width) {
        out <- c(out, substr(line, 1, max_width))
        line <- substr(line, max_width + 1, nchar(line))
      }
      out <- c(out, line)
    }
    out
  }

  pad_lines <- function(lines, max_lines) {
    c(lines, rep('', max_lines - length(lines)))
  }

  df[, 'popular'] = ifelse(df[, 'popular'], 'o', '')
  out <- sprintf('|%s|\n', paste(stringr::str_pad(names(DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS), DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS, side='right', pad='-'), collapse='|'))
  empty_row <- sprintf('|%s|\n', paste(stringr::str_pad(rep('', length(DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS)), DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS), collapse='|'))

  r <- list()
  for (i in seq_len(dim(df)[1])) {
    print(df[i, ])
    for (c in names(DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS)) {
      r[[c]] <- if (df[i, c] != '') line_wrap(df[i, c], DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS[[c]]) else character(0)
    }

    max_lines <- max(sapply(r, length))
    for (c in names(DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS)) {
      r[[c]] <- pad_lines(r[[c]], max_lines)
    }

    for (j in seq_len(max_lines)) {
      out <- sprintf('%s|%s|\n', out, paste(stringr::str_pad(sapply(r, "[[", j), DISPLAY_TEXT_KW_INDEX_MAX_WIDTHS, side='right'), collapse='|'))
    }
    out <- paste(out, empty_row, sep='')
  }

  sprintf('\\preformatted{%s}', out)
}

.formatKeywordsLatex <- function(df) {
  sprintf('\\out{%s}',
    knitr::kable(df, 'latex', escape=FALSE, longtable=TRUE, align=sprintf('p{%.1f\\textwidth}', DISPLAY_LATEX_KW_INDEX_PCT_WIDTHS), vline="") %>%
      gsub(' *\n *', ' ', .) %>%
      gsub('\\\\ ', '\\\\\\\\ ', .))
}

.formatKeywordsHtml <- function(df) {
  sprintf('\\out{%s}', knitr::kable(df, 'html', escape=FALSE, table.attr='class="termtable"'))
}

#' Keywords defined for Exponential-Family Random Graph Models
#'
#' @name ergmKeyword
#' @aliases ergm-keywords keywords-ergm ergm.keywords keywords.ergm
#' @description This collects all defined keywords defined for the ERGM and derived packages
#'
#' @section Possible keywords defined by the ERGM and derived packages:
#'
#' \ergmCSS
#'
#' \if{latex}{\Sexpr[results=rd,stage=render]{ergm:::.formatKeywordsLatex(ergm::ergm_keyword())}}
#' \if{text}{\Sexpr[results=rd,stage=render]{ergm:::.formatKeywordsText(ergm::ergm_keyword())}}
#' \if{html}{\Sexpr[results=rd,stage=render]{ergm:::.formatKeywordsHtml(ergm::ergm_keyword())}}
#'
#' @keywords models
NULL
