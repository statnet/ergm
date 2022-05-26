#  File R/ergm-terms-index.R in package ergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2003-2022 Statnet Commons
################################################################################

SUPPORTED_TERM_TYPES <- c('Term', 'Constraint', 'Reference', 'Hint', 'Proposal')
SUPPORTED_TERM_TYPE_REGEX <- sprintf('-ergm(%s)(.Rd)?', paste(SUPPORTED_TERM_TYPES, collapse='|'))

DISPLAY_TEXT_INDEX_MAX_WIDTHS <- list('Term'=25, 'Pkg'=5, 'Description'=33, 'Concepts'=12)
DISPLAY_TEXT_MAX_WIDTH <- sum(unlist(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) + length(DISPLAY_TEXT_INDEX_MAX_WIDTHS) - 1
DISPLAY_LATEX_INDEX_PCT_WIDTHS <- c(0.3, 0.05, 0.5, 0.1)
DISPLAY_LATEX_TOC_PCT_WIDTHS <- function(n_concepts) c(2.4, rep(.7, n_concepts))

.fsub <- function(x, pattern, replacement, fixed = TRUE, ...) gsub(pattern, replacement, x, fixed = fixed, ...)

# Return the index entry for a single term in the new format
.parseTerm <- function(name, pkg, pkg_name) {
  doc <- pkg[[name]]
  name <- substr(name, 1, nchar(name) - 3)
  tags <- sapply(doc, attr, 'Rd_tag')
  comps <- strsplit(name, '-')[[1]]
  type <- comps[length(comps)]

  if (type != "ergmProposal") {
    concepts <- doc[tags == '\\concept'] %>% unlist
    raw_usage <- doc[tags == '\\usage'] %>%
      unlist %>%
      paste(collapse="") %>%
      gsub("\n *# *"," ", .) %>%
      strsplit("\n") %>%
      .[[1]] %>%
      trimws()

    if (type %in% c('ergmTerm')) {
      usages <- list()
      for (usage_line in regmatches(raw_usage, regexec("^(binary|valued): *(.+)$", raw_usage))) {
        usages[[length(usages) + 1]] <- list('type'=usage_line[2], 'usage'=usage_line[3])
        concepts <- c(concepts, usage_line[2])
      }
    } else {
      usages <- lapply(raw_usage, function(u) list(type=NULL, usage=u))
    }

    if (length(usages) == 0) {
      return(NULL)
    }

    return(list(
      link=name,
      name=comps[1],
      type=type,
      alias=doc[tags == '\\alias'] %>% unlist,
      package=pkg_name,
      usages=usages,
      title=doc[tags == '\\title'] %>% unlist %>% paste(collapse='') %>% trimws(),
      description=doc[tags == '\\description'] %>% unlist %>% paste(collapse='') %>% trimws(),
      details=doc %>% unlist %>% paste(collapse='') %>% trimws(),
      concepts=if (!is.null(concepts)) unique(concepts) else c(),
      keywords=c()))
  } else {
    ps <- ergm_proposal_table()
    ps <- .filterProposals(ps, proposal=comps[1])

    proposals = list()
    for (i in seq_len(nrow(ps))) {
      if (!stringr::str_detect(ps$Constraints[i], '[|&]') || stringr::str_detect(ps$Constraints[i], '\\+')) {
        constraints <- strsplit(ps$Constraints[i], '\\+')[[1]]
        constraints <- paste0(ifelse(constraints == '.dyads', '|', '&'), constraints)
      } else {
        constraints <- strsplit(ps$Constraints[i], '(?<=.)(?=[|&])', perl=TRUE)[[1]]
      }
      constraints <- lapply(constraints, function(c) list(
        name=substr(c, 2, stringr::str_length(c)),
        enforce=substr(c, 1, 1) == '&'))
      proposals[[length(proposals) + 1]] <- list(
        Proposal=ps$Proposal[i],
        Reference=ps$Reference[i],
        Enforces=constraints %>% keep("enforce") %>% map("name") %>% unlist,
        May_Enforce=constraints %>% discard("enforce") %>% map("name") %>% unlist,
        Priority=ps$Priority[i],
        Weight=ps$Weights[i],
        Class=ifelse(ps$Class[i] == 'c', 'cross-sectional', 'last-toggle')
      )
    }
    return(list(
      link=name,
      name=comps[1],
      type=type,
      alias=doc[tags == '\\alias'] %>% unlist,
      package=pkg_name,
      title=doc[tags == '\\title'] %>% unlist %>% paste(collapse='') %>% trimws(),
      description=doc[tags == '\\description'] %>% unlist %>% paste(collapse='') %>% trimws(),
      details=doc %>% unlist %>% paste(collapse='') %>% trimws(),
      rules=proposals))
  }
}

#' A simple dictionary to cache loaded terms
#'
#' @param the name of the term
#' @param env the environment name for the function; if `NULL`, look
#'   up in cache, otherwise insert or overwrite.
#'
#' @return A character string giving the name of the environment
#'   containing the function, or `NULL` if not in cache.
#' @noRd
ergmTermCache <- local({
  cache <- lapply(SUPPORTED_TERM_TYPES, function(x) list())
  names(cache) <- paste('ergm', SUPPORTED_TERM_TYPES, sep='')
  pkglist <- character(0) # Current list of packages.

  # Reset the cache and update the list of watched packages.
  unload <- function(pkg_name, ...) {
    pkglist <<-.packages(TRUE)

    for (term_type in SUPPORTED_TERM_TYPES) {
      for (term in names(cache[[term_type]])) {
        if (cache[[term_type]][[term]]$package == pkg_name) {
          cache[[term_type]][[term]] <<- NULL
        }
      }
    }
  }

  # Crawl all loaded packages for terms, parse them once and store in the singleton store
  load <- function(pkg_name) {
    pkg <- tools::Rd_db(pkg_name)
    all_doco <- attributes(pkg)$names
    converted <- all_doco[grep(SUPPORTED_TERM_TYPE_REGEX, all_doco)]

    for (term in lapply(converted, .parseTerm, pkg, pkg_name)) {
      if (!is.null(term)) {
        cache[[term$type]][[term$link]] <<- term
      }
    }
  }

  # Check if new namespaces have been added.
  checknew <- function() {
    loaded_packages <- loadedNamespaces()
    revdeps <- c("ergm", tools::dependsOnPkgs("ergm"))

    term_packages <- tryCatch({
      db <- utils::hsearch_db(package=revdeps)$Base
      unique(db$Package[grep(SUPPORTED_TERM_TYPE_REGEX, db$Topic)])
    }, error = function(e) {
      message("Error querying the list of term packages when indexing: ", sQuote(toString(e)), ";  skipping.")
      character(0)
    })

    for (pkg_name in intersect(loaded_packages, term_packages)) {
      if (!pkg_name %in% pkglist) {
        load(pkg_name)

        setHook(packageEvent(pkg_name, "detach"), unload)
        setHook(packageEvent(pkg_name, "onUnload"), unload)
      }
    }
    cache <<- lapply(cache, function(terms) terms[sort(names(terms))])

    pkglist <<- loaded_packages
  }

  function (term_type) {
    checknew()

    cache[[term_type]]
  }
})

#' Filter a list of terms, currently by keywords/concepts
#'
#' @param terms a term list returned by `ergmTermCache()`
#' @param keywords a function with one argument or a formula understood by [purrr::as_mapper()] that takes a character vector containing a term's keywords/concepts and returns `TRUE` or `FALSE`
#' @param packages a character vector containing the packages to search through
#' @param ... additional arguments, currently unused
#'
#' @return a filtered term list
#' @noRd
.filterTerms <- function(terms, keywords = NULL, packages=NULL, ...) {
  if (!is.null(keywords)) {
    keywords <- as_mapper(keywords)
    keep <- terms %>% map("concepts") %>% map_lgl(keywords)
    terms <- terms[keep]
  }

  if (!is.null(packages)) {
    terms <- terms[(terms %>% map("package")) %in% packages]
  }

  terms
}

# Parse a usage string and remove the default arguments from the string
.removeDefaultArguments <- function(usage) {
  expr <- str2lang(usage)

  if (!inherits(expr, 'call')) {
    return(usage)
  } else {
    expr <- NVL3(names(expr), ifelse(. != "", ., as.character(expr)), as.character(expr))
    return(sprintf("%s(%s)", expr[1], paste(expr[2:length(expr)], collapse=", ")))
  }
}

#' Constructs a data frame containing term information, suitable for typesetting in help files and vignettes
#'
#' @param term_type character string giving the type of term, currently `"ergmTerm"`, `"ergmConstraint"`, or `"ergmReference"`
#' @param ... further arguments passed to `.filterTerms()`
#' @param display.keywords keywords to not put into the table column (usually because they are redundant)
#'
#' @return a data frame with columns for usage, package, title, concepts, and link
#' @noRd
.buildTermsDataframe <- function(term_type, ..., display.keywords = setdiff(ergm::ergm_keyword()$name, c("binary", "valued"))) {
  terms <- ergmTermCache(term_type)

  terms <- .filterTerms(terms, ...)

  if (length(terms) == 0) return(NULL)

  sapply(terms, function(term) {
    if (!is.null(term$usages[[1]]$type)) {
      usage <- paste(sprintf('`%s` (%s)',
        sapply(term$usages, "[[", 'usage') %>% sapply(., .removeDefaultArguments) %>% trimws,
        sapply(term$usages, "[[", 'type')), collapse='\n')
    } else {
      usage <- paste(sprintf('`%s`', sapply(term$usages, '[[', 'usage')), collapse='\n')
    }

    term$concepts <- intersect(term$concepts, display.keywords)

    c(Term = usage, Package = term$package, Description = term$title, Concepts = if (length(term$concepts) > 0) paste(term$concepts, collapse='\n') else '', Link = term$link)
  }) %>% t() %>% as.data.frame(stringsAsFactors=FALSE)
}

#' Constructs a data frame containing term proposals, suitable for typesetting in help files and vignettes
#'
#' @param proposal only include rules relevant to the given proposal
#'
#' @return a data frame with columns for usage, package, title, concepts, and link
#' @noRd
.buildProposalsList <- function(proposal) {
  proposals <- ergmTermCache("ergmProposal")
  if (!missing(proposal)) {
    proposals <- proposals[[paste0(proposal, '-ergmProposal')]]$rules
  } else {
    proposals <- proposals %>% map("rules") %>% flatten()
  }

  if (length(proposals) == 0) return(NULL)
  names(proposals) <- seq_along(proposals)
  proposals
}


# terms : a list structure of the documentation data
# ... : further arguments passed to `.filterTerms()`
# display.keywords : an optional vector of column names to print and include
.termMatrix <- function(term_type, ..., display.keywords=NULL) {
  terms <- ergmTermCache(term_type)

  terms <- .filterTerms(terms, ...)

  # if list of display.keywords not supplied, generate it
  # otherwise, use the display.keywords (and column order) provided
  if (is.null(display.keywords)) {
    display.keywords <- unique(unlist(sapply(terms,'[[','concepts')))
  }

  if(length(display.keywords)==0) return(NULL)

  if (length(terms) == 0) return(NULL)

  # figure out which terms are members of each cat
  membership <- lapply(display.keywords, function(cat) {
    # return true for terms that match cat
    sapply(terms, function(term) cat %in% term$concepts)
  })

  df <- data.frame(membership)

  concepts <- ergm_keyword()
  display.keywords <- concepts[match(display.keywords, concepts$name), 'short']
  colnames(df) <- display.keywords
  rownames(df) <- NULL
  df$Link <- sapply(terms,'[[','link')
  df$Term <- sapply(terms,'[[','name')

  df[c('Term', display.keywords, 'Link')]
}

# output listings of terms, grouped by keywords
.termToc <- function(term_type, ..., display.keywords=NULL) {
  terms <- ergmTermCache(term_type)

  terms <- .filterTerms(terms, ...)

  # if list of display.keywords not supplied, generate it
  # otherwise, use the display.keywords (and column order) provided
  if (is.null(display.keywords)) {
    display.keywords <- unique(unlist(sapply(terms,'[[','concepts')))
  }

  if(length(display.keywords)==0) return(NULL)

  if (length(terms) == 0) return(NULL)

  ret <- list()
  for (cat in display.keywords) {
    #find ones with matching cats
    matchingTerms <- terms[sapply(terms, function(term) any(term$concepts==cat))]
    ret[[cat]] <- list(
      link=sapply(matchingTerms, '[[', 'link'),
      name=sapply(matchingTerms, '[[', 'name')
    )
  }

  ret
}

.filterProposals <- function(proposals, proposal=NULL, ...) {
  if (!is.null(proposal)) {
    proposals <- proposals[proposals$Proposal == proposal,]
  }

  proposals
}

PROPOSAL_NOT_IN_TABLE <- "This proposal is not referenced in the lookup table."

.formatProposalsHtml <- function(df, keepProposal=FALSE) {
  if (NROW(df) == 0) return(paste0("\\out{<p>", PROPOSAL_NOT_IN_TABLE, "</p>}"))

  for (i in seq_along(df)) {
    df[[i]]$Proposal <- sprintf('<a href="../help/%1$s-ergmProposal">%1$s</a>', df[[i]]$Proposal)
    df[[i]]$Reference <- sprintf('<a href="../help/%1$s-ergmReference">%1$s</a>', df[[i]]$Reference)
    df[[i]]$Enforces <- if (length(df[[i]]$Enforces) > 0) paste(sprintf('<a href="../help/%1$s-ergmConstraint">%1$s</a>', df[[i]]$Enforces), collapse=' ') else ""
    df[[i]]$May_Enforce <- if (length(df[[i]]$May_Enforce) > 0) paste(sprintf('<a href="../help/%1$s-ergmConstraint">%1$s</a>', df[[i]]$May_Enforce), collapse=' ') else ""
  }

  df <- do.call(rbind.data.frame, df)

  if (!keepProposal) {
    df <- df[,colnames(df)!="Proposal"]
  }

  css <- '<style>th,td {padding:3px 10px}</style>'
  sprintf("\\out{%s%s}", css, knitr::kable(df, 'html', escape=FALSE, row.names=FALSE))
}

.formatProposalsLatex <- function(df, keepProposal=FALSE) {
  if (NROW(df) == 0) return(paste0("\\out{",PROPOSAL_NOT_IN_TABLE,"}"))

  for (i in seq_along(df)) {
    df[[i]]$Enforces <- if (length(df[[i]]$Enforces) > 0) paste(df[[i]]$Enforces, collapse=' ') else ""
    df[[i]]$May_Enforce <- if (length(df[[i]]$May_Enforce) > 0) paste(df[[i]]$May_Enforce, collapse=' ') else ""
  }

  df <- as.data.frame(do.call(rbind.data.frame, df))
  names(df) <- c('Proposal', 'Reference', 'Enforces', 'May Enforce', 'Priority', 'Weight', 'Class')
  df[, 'Class'] = ifelse(df[, 'Class'] == 'cross-sectional', 'c', 't')

  if (!keepProposal) {
    df <- df[,colnames(df)!="Proposal"]
  }

  sprintf("\\out{%s}", knitr::kable(df, 'latex', escape=FALSE, row.names=FALSE, longtable=TRUE, vline="") %>%
    gsub(' *\n *', ' ', .) %>%
    gsub('\\\\ ', '\\\\\\\\ ', .))

}

.formatProposalsText <- function(df, keepProposal=FALSE) {
  if (NROW(df) == 0) return(PROPOSAL_NOT_IN_TABLE)

  for (i in seq_along(df)) {
    df[[i]]$Enforces <- if (length(df[[i]]$Enforces) > 0) paste(df[[i]]$Enforces, collapse=' ') else ""
    df[[i]]$May_Enforce <- if (length(df[[i]]$May_Enforce) > 0) paste(df[[i]]$May_Enforce, collapse=' ') else ""
  }

  df <- do.call(rbind.data.frame, df)

  if (!keepProposal) {
    df <- df[,colnames(df)!="Proposal"]
  }

  sprintf('\\preformatted{%s}', paste(knitr::kable(df, 'pipe'), collapse='\n'))
}

# output listings of terms, grouped by keywords
#' @importFrom magrittr "%>%" "%<>%"
.formatIndexText <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term %<>% .fsub('binary', 'bin') %>% .fsub('valued', 'val') %>% .fsub('`', '')

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

  colnames(df)[[2]] <- 'Pkg'
  out <- sprintf('|%s|\n', paste(stringr::str_pad(names(DISPLAY_TEXT_INDEX_MAX_WIDTHS), DISPLAY_TEXT_INDEX_MAX_WIDTHS, side='right', pad='-'), collapse='|'))
  empty_row <- sprintf('|%s|\n', paste(stringr::str_pad(rep('', length(DISPLAY_TEXT_INDEX_MAX_WIDTHS)), DISPLAY_TEXT_INDEX_MAX_WIDTHS), collapse='|'))

  r <- list()
  for (i in seq_len(dim(df)[1])) {
    for (c in names(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) {
      r[[c]] <- if (df[i, c] != '') line_wrap(df[i, c], DISPLAY_TEXT_INDEX_MAX_WIDTHS[[c]]) else c()
    }

    max_lines <- max(sapply(r, length))
    for (c in names(DISPLAY_TEXT_INDEX_MAX_WIDTHS)) {
      r[[c]] <- pad_lines(r[[c]], max_lines)
    }

    for (j in seq_len(max_lines)) {
      out <- sprintf('%s|%s|\n', out, paste(stringr::str_pad(sapply(r, "[[", j), DISPLAY_TEXT_INDEX_MAX_WIDTHS, side='right'), collapse='|'))
    }
    out <- paste(out, empty_row, sep='')
  }

  sprintf('\\preformatted{%s}', out)
}

.formatMatrixText <- function(df) {
  if(is.null(df)) return(NULL)

  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  sprintf('\\preformatted{%s}', paste(knitr::kable(df, 'pipe'), collapse='\n'))
}

.formatTocText <- function(toc) {
  if(is.null(toc)) return(NULL)

  out <- ''
  for (cat in names(toc)) {
    out <- sprintf('%s\n\n%s:\n%s', out, cat,
      paste(paste('   ', stringr::str_wrap(paste(toc[[cat]]$name, collapse=', '), DISPLAY_TEXT_MAX_WIDTH)), collapse='\n'))
  }
  sprintf('\\preformatted{%s}', out)
}

# Format the table for text output
.formatIndexLatex <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term %<>% .fsub('binary', 'bin') %>% .fsub('valued', 'val')
  df$Term %<>% .fsub('\n', ' \\newline ') %>% gsub('`([^`]*)`', '\\1', .) %>%
    strsplit(' ') %>% sapply(., function(x) x %>% .fsub('_', '\\_') %>% sprintf('\\code{%s}', .) %>% paste(collapse=' ')) %>%
    gsub('\\\\code\\{\\(([^(]*)\\)\\}', '(\\1)', .) %>% .fsub('\\code{\\newline}', '\\newline') %>%
    paste('\\\\raggedright \\\\allowbreak', .)
  df$Link <- NULL
  sprintf('\\out{%s}',
    knitr::kable(df, 'latex', escape=FALSE, row.names=FALSE, longtable=TRUE, align=sprintf('p{%.1f\\textwidth}', DISPLAY_LATEX_INDEX_PCT_WIDTHS), vline="") %>%
      gsub(' *\n *', ' ', .) %>%
      gsub('\\\\ ', '\\\\\\\\ ', .))
}

.formatMatrixLatex <- function(df) {
  if(is.null(df)) return(NULL)

  df$Link <- NULL

  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], 'o', '')
  }

  sprintf('\\out{%s}',
    knitr::kable(df, 'latex', escape=FALSE, row.names=FALSE, longtable=TRUE, align=sprintf('p{%.1fcm}', DISPLAY_LATEX_TOC_PCT_WIDTHS(dim(df)[2] - 1)), vline="") %>%
      gsub(' *\n *', ' ', .) %>%
      gsub('\\\\ ', '\\\\\\\\ ', .))
}

.formatTocLatex <- function(toc) {
  if(is.null(toc)) return(NULL)

  out <- '\\out{\\begin{description}'
  for (cat in names(toc)) {
    out <- sprintf('%s\n\\item[%s] %s', out, cat, paste(toc[[cat]]$name, collapse=', '))
  }
  paste(out, "\n\\end{description}}")
}

# Format the table for text output
.formatIndexHtml <- function(df) {
  if(is.null(df)) return(NULL)

  df$Term %<>% .fsub('binary', 'bin') %>% .fsub('valued', 'val')

  # Hack! because HTML code has to be wrapped \out{} which prevents \link{} from being parsed, the link has
  # to be constructed with the actual address. To find the address, generate the correct link with 
  # \link[=absdiff-ergmTerm]{test} and check that it works with \out{<a href="../help/absdiff-ergmTerm">test</a>}.
  # This address may change from an upstream R-studio change

  df$Term <- sprintf(gsub('`([^`(]*)([^`]*)`', '<span class="code"><a href="../help/%1$s" id="%1$s">\\1\\2</a></span>', .fsub(df$Term, '\n', '<br />\n')), df$Link)
  df$Link <- NULL

  css <- '<style>.striped th,.striped td {padding:3px 10px} .striped tbody tr:nth-child(odd) {background: #eee} .striped .code {font-family: monospace} .matrix td {align: center} .matrix th,.matrix td {padding-right:5px; width: 75px}</style>'
  sprintf('\\out{%s%s}', css, knitr::kable(df, 'html', escape=FALSE, row.names=FALSE, table.attr='class="striped"'))
}

.formatMatrixHtml <- function(df, wrapRdTags=TRUE) {
  if(is.null(df)) return(NULL)

  df$Term <- sprintf('<a href="#%s">%s</a>', df$Link, df$Term)
  df$Link <- NULL
  for (c in colnames(df)[-1]) {
    df[[c]] <- ifelse(df[[c]], '&#10004;', '')
  }

  out <- knitr::kable(df, 'html', escape=FALSE, row.names=FALSE, table.attr='class="matrix"')
  if(wrapRdTags) {
    out <- sprintf('\\out{%s}', out)
  }
  out
}

.formatTocHtml <- function(toc, wrapRdTags=TRUE) {
  if(is.null(toc)) return(NULL)

  out <- paste('Jump to keyword:', paste(sprintf('<a href="#cat_%s">%s</a>', gsub(' ', '_', names(toc)), names(toc)), collapse=' '))
  for (cat in names(toc)) {
    out <- sprintf('%s<h3><a id="cat_%s">%s</a></h3>%s', out, gsub(' ', '_', cat), cat,
      paste(sprintf('<a href="#%s">%s</a>', toc[[cat]]$link, toc[[cat]]$name), collapse=' '))
  }

  if(wrapRdTags) {
    out <- sprintf('\\out{%s}', out)
  }

  out
}

search.ergmTermType <-function(term.type, search, net, keywords, name, packages) {
  if (!missing(net)){
    if(!is.network(net)){
      stop("the 'net' argument must be the network argument that applicable terms are to be searched for")
    }
  }

  if(missing(keywords)){
    keywords<-character(0)
  }
  if (!missing(net)){
    if(is.directed(net)){
      keywords<-c(keywords,'directed')
    } else {
      keywords<-c(keywords,'undirected')
    }
    if(is.bipartite(net)){
      keywords<-c(keywords,'bipartite')
    }
  }

  terms <- ergmTermCache(term.type)

  found<-rep(TRUE,length(terms))

  # if name is specified, restrict to terms with that name
  if(!missing(name)){
    for (t in seq_along(terms)){
      term<-terms[[t]]
      found[t]<-any(term$name==name || term$link==name)
    }
  }

  # restrict by keywords
  for (t in which(found)){
    term<-terms[[t]]

    if(!all(keywords%in%term$concepts)){
      found[t]<-FALSE
    }
  }

  # next (optionally) restrict by search matches
  if (!missing(search)){
    for (t in which(found)){
      term<-terms[[t]]
      # if we don't find the search in the text grep, mark it as false
      if(length(grep(search,term$details,ignore.case=TRUE))==0){
        found[t]<-FALSE
      }
    }
  }

  # optionally restrict by package
  if (!missing(packages)) {
    for (t in which(found)) {
      term <- terms[[t]]
      if (!term$package %in% packages) {
        found[t] <- FALSE
      }
    }
  }

  # if term name was specified, print out all the matching terms
  # otherwise,  loop over the remaining terms to format output as condensed
  output<-list()
  if(!missing(name)){
    if(sum(found)==0){
      cat("No terms named '",name,"' were found. Try searching with search='",name,"'instead.\n",sep='')
    } else {
      cat("Definitions for term(s) ",name,":\n")
      for (t in which(found)){
        term<-terms[[t]]
        output<-c(output, term)
        cat(sprintf('%s\n    %s: %s\n    Keywords: %s\n\n',
          term$usages[[1]]$usage,
          term$title,
          term$description,
          paste(term$concepts, collapse=', ')))
      }
    }
  }else{
    for (t in which(found)){
      term<-terms[[t]]
      unique_usages <- term$usages %>% map("usage") %>% unlist %>% unique
      if (term$type %in% c('ergmTerm')) {
        if (!missing(keywords) && any(c("binary", "valued") %in% keywords)) {
          term_type_to_match <- intersect(keywords, c('binary', 'valued'))
          term_types <- sapply(unique_usages, function(usage) paste(term$usages[(term$usages %>% map("usage") %>% unlist) == usage] %>% map("type") %>% unlist %>% unique %>% intersect(term_type_to_match), collapse=", "))

          unique_usages <- unique_usages[term_types != '']
          term_types <- term_types[term_types != '']
        } else {
          term_types <- sapply(unique_usages, function(usage) paste(term$usages[(term$usages %>% map("usage") %>% unlist) == usage] %>% map("type") %>% unlist %>% unique, collapse=", "))
        }

        unique_usages <- paste0(unique_usages, " (", term_types, ")")
      }
      outText <- paste0(paste(unique_usages, collapse="\n"), "\n    ", term$title, "\n")
      output<-c(output,outText)
    }
    cat("Found ",length(output)," matching ergm terms:\n")
    cat(paste(output,collapse='\n'))
  }
  invisible(output)
}

# function to look up the set of terms applicable for a specific network

#' Search the ERGM terms
#' 
#' Searches through the \code{\link{ergm.terms}} help page and prints out a
#' list of terms appropriate for the specified network's structural
#' constraints, optionally restricting by additional keywords and search term
#' matches.
#' 
#' Uses \code{\link{grep}} internally to match the search terms against the term
#' description, so \code{search} is currently matched as a single phrase.
#' Keyword tags will only return a match if all of the specified tags are
#' included in the term.
#' 
#' @param search optional character search term to search for in the text of the
#' term descriptions. Only matching terms will be returned. Matching is case
#' insensitive.
#' @param net a network object that the term would be applied to, used as
#' template to determine directedness, bipartite, etc
#' @param keywords optional character vector of keyword tags to use to
#' restrict the results (i.e. 'curved', 'triad-related')
#' @param name optional character name of a specific term to return
#' @param packages optional character vector indicating the subset of packages in which to search
#' @return prints out the name and short description of matching terms, and
#' invisibly returns them as a list.  If \code{name} is specified, prints out
#' the full definition for the named term.
#' @author skyebend@uw.edu
#' @seealso See also \code{\link{ergm.terms}} for the complete documentation
#' @examples
#' \donttest{
#' # find all of the terms that mention triangles
#' search.ergmTerms('triangle')
#' 
#' # two ways to search for bipartite terms:
#' 
#' # search using a bipartite net as a template
#' myNet<-network.initialize(5,bipartite=3)
#' search.ergmTerms(net=myNet)
#' 
#' # or request the bipartite keyword
#' search.ergmTerms(keywords='bipartite')
#' 
#' # search on multiple keywords
#' search.ergmTerms(keywords=c('bipartite','dyad-independent'))
#' 
#' # print out the content for a specific term
#' search.ergmTerms(name='b2factor')
#'
#' # request the bipartite keyword in the ergm package
#' search.ergmTerms(keywords='bipartite', packages='ergm')
#' }
#' @importFrom utils capture.output
#' @export search.ergmTerms
search.ergmTerms <- function(search, net, keywords, name, packages) {
  search.ergmTermType("ergmTerm", search, net, keywords, name, packages)
}

#' Search the ERGM References
#'
#' Searches through the \code{\link{ergm.references}} help page and prints out a
#' list of terms appropriate for the specified network's structural
#' constraints, optionally restricting by additional keywords and search term
#' matches.
#'
#' Uses \code{\link{grep}} internally to match the search terms against the term
#' description, so \code{search} is currently matched as a single phrase.
#' Keyword tags will only return a match if all of the specified tags are
#' included in the term.
#'
#' @param search optional character search term to search for in the text of the
#' term descriptions. Only matching terms will be returned. Matching is case
#' insensitive.
#' @param keywords optional character vector of keyword tags to use to
#' restrict the results (i.e. 'curved', 'triad-related')
#' @param name optional character name of a specific term to return
#' @param packages optional character vector indicating the subset of packages in which to search
#' @return prints out the name and short description of matching terms, and
#' invisibly returns them as a list.  If \code{name} is specified, prints out
#' the full definition for the named term.
#' @author skyebend@uw.edu
#' @seealso See also \code{\link{ergm.references}} for the complete documentation
#' @examples
#' \donttest{
#' # find all of the references that mention dyad
#' search.ergmReferences('dyad')
#'
#' # search on multiple keywords
#' search.ergmReferences(keywords=c('binary','discrete'))
#'
#' # print out the content for a specific reference
#' search.ergmReferences(name='Bernoulli')
#'
#' # request the bipartite keyword in the ergm package
#' search.ergmReferences(keywords='binary', packages='ergm')
#' }
#' @importFrom utils capture.output
#' @export search.ergmReferences
search.ergmReferences <- function(search, keywords, name, packages)
  search.ergmTermType("ergmReference", search=search, keywords=keywords, name=name, packages=packages)

#' Search the ERGM Constraint
#'
#' Searches through the \code{\link{ergm.constraints}} help page and prints out a
#' list of constraints appropriate for the specified network's structural
#' constraints, optionally restricting by additional keywords and search term
#' matches.
#'
#' Uses \code{\link{grep}} internally to match the search terms against the constraint
#' description, so \code{search} is currently matched as a single phrase.
#' Keyword tags will only return a match if all of the specified tags are
#' included in the constraint
#'
#' @param search optional character search term to search for in the text of the
#' constraint descriptions. Only matching constraints will be returned. Matching is case
#' insensitive.
#' @param keywords optional character vector of keyword tags to use to
#' restrict the results (i.e. 'curved', 'triad-related')
#' @param name optional character name of a specific constraint to return
#' @param packages optional character vector indicating the subset of packages in which to search
#' @return prints out the name and short description of matching constraints, and
#' invisibly returns them as a list.  If \code{name} is specified, prints out
#' the full definition for the named constraint.
#' @author skyebend@uw.edu
#' @seealso See also \code{\link{ergm.constraints}} for the complete documentation
#' @examples
#' \donttest{
#' # find all of the constraint that mention degrees
#' search.ergmConstraints('degree')
#'
#' # search on multiple keywords
#' search.ergmConstraints(keywords=c('directed','dyad-independent'))
#'
#' # print out the content for a specific constraint
#' search.ergmConstraints(name='b1degrees')
#'
#' # request the bipartite keyword in the ergm package
#' search.ergmConstraints(keywords='directed', packages='ergm')
#' }
#' @importFrom utils capture.output
#' @export search.ergmConstraints
search.ergmConstraints <- function(search, keywords, name, packages)
  search.ergmTermType("ergmConstraint", search=search, keywords=keywords, name=name, packages=packages)

#' Search the ERGM Proposals
#'
#' Searches through the \code{\link{ergm.proposals}} help page and prints out a
#' list of proposals appropriate for the specified network's structural
#' constraints, optionally restricting by additional keywords and search term
#' matches.
#'
#' Uses \code{\link{grep}} internally to match the search terms against the proposal
#' description, so \code{search} is currently matched as a single phrase.
#' Keyword tags will only return a match if all of the specified tags are
#' included in the proposal
#'
#' @param search optional character search term to search for in the text of the
#' term descriptions. Only matching terms will be returned. Matching is case
#' insensitive.
#' @param name optional character name of a specific proposal to return
#' @param reference optional character to limit proposals to only those with a matching reference
#' @param constraints optional character vector to limit proposals to only those with matching hard or soft constraints
#' @param packages optional character vector indicating the subset of packages in which to search
#' @return prints out the name and short description of matching terms, and
#' invisibly returns them as a list.  If \code{name} is specified, prints out
#' the full definition for the named term.
#' @seealso See also \code{\link{ergm.proposals}} for the complete documentation
#' @examples
#' \donttest{
#' # find all of the proposals that mention triangles
#' search.ergmProposals('MH algorithm')
#'
#' # print out the content for a specific proposals
#' search.ergmProposals(name='randomtoggle')
#'
#' # find all proposals with required or optional constraints
#' search.ergmProposals(constraints='.dyads')
#'
#' # find all proposals with references
#' search.ergmProposals(reference='Bernoulli')
#'
#' # request proposals that mention triangle in the ergm package
#' search.ergmProposals('MH algorithm', packages='ergm')
#' }
#' @importFrom utils capture.output
#' @export search.ergmProposals
search.ergmProposals <- function(search, name, reference, constraints, packages) {
  proposals <- .buildProposalsList()
  
  terms <- ergmTermCache("ergmProposal")
  
  found <- rep(TRUE,length(terms))

  # if name is specified, restrict to terms with that name
  if (!missing(name)) {
    for (t in seq_along(terms)) {
      term<-terms[[t]]
      found[t]<-any(term$name==name || term$link==name)
    }
  }

   # next (optionally) restrict by search matches
  if (!missing(search)) {
    for (t in which(found)) {
      term <-terms[[t]]
      # if we don't find the search in the text grep, mark it as false
      if (length(grep(search,term$details, ignore.case=TRUE)) == 0) {
        found[t]<-FALSE
      } 
    }
  }

  # optionally restrict by references
  if (!missing(reference)) {
    for (t in which(found)) {
      term <-terms[[t]]
      if (! reference %in% (term$rules %>% map("Reference"))) {
        found[t]<-FALSE
      }
    }
  }

  # optionally restrict by constraints
  if (!missing(constraints)) {
    for (constraint in constraints) {
      for (t in which(found)) {
        term <-terms[[t]]
        if (!constraint %in% (term$rules %>% map('Enforces') %>% unlist) && !constraint %in% (term$rules %>% map("May_Enforce") %>% unlist)) {
          found[t]<-FALSE
        }
      }
    }
  }

  # optionally restrict by package
  if (!missing(packages)) {
    for (t in which(found)) {
      term <- terms[[t]]
      if (!term$package %in% packages) {
        found[t] <- FALSE
      }
    }
  }
  
  # if term name was specified, print out all the matching terms
  # otherwise,  loop over the remaining terms to format output as condensed
  output<-list()
  if(!missing(name)){
    if(sum(found)==0){
      cat("No proposals named '",name,"' were found. Try searching with search='",name,"'instead.\n",sep='')
    } else {
      cat("Definitions for proposal(s) ",name,":\n")
      for (t in which(found)){
        term<-terms[[t]]
        output<-c(output, list(term))
        cat(sprintf('%s\n    %s: %s\n',
          term$name,
          term$title,
          term$description))
        for (rule in term$rules) {
          cat(sprintf('    Reference: %s Class: %s\n%s%s\n',
            rule$Reference,
            rule$Class,
            if (length(rule$Enforces) > 0) paste("    Enforces:", paste(rule$Enforces, collapse=" "), "\n") else "",
            if (length(rule$May_Enforce) > 0) paste("    May Enforce:", paste(rule$May_Enforce, collapse=" "), "\n") else ""))
        }
      }
    }
  }else{
    for (t in which(found)){
      term<-terms[[t]]
      outText <- sprintf('%s\n    %s\n', term$name, term$title)
      output<-c(output,outText)
    }
    cat("Found ",length(output)," matching ergm proposals:\n")
    cat(paste(output,collapse='\n'))
  }
  invisible(output)
}
