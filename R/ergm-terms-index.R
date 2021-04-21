#  File R/ergm-terms-index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons
#######################################################################

# Cleans up the usage clause in the Roxygen documentation and return a list of length 3 vectors (full string, binary|valued, usage)
.parseUsage <- function(raw_usage) {
    usage <- raw_usage %>% 
        unlist %>% 
        paste(collapse="") %>%
        gsub("\n *# *"," ", .) %>%
        strsplit("\n") %>%
        .[[1]] %>%
        trimws()
    return(regmatches(usage, regexec("^(binary|valued): *(.+)$", usage)))
}

# Return the index entry for a single term in the new format
.indexTerm <- function(name, pkg, pkg_name) {
    doc <- pkg[[name]]
    tags <- sapply(doc, attr, 'Rd_tag')
    name <- .parseUsage(doc[tags == '\\usage'])[[1]][[3]]
    ret <- sprintf('%s (%s)\n    %s\n\n    Keywords: %s\n    %s\n\n',
                   name,
                   pkg_name,
                   doc[tags == (if ('\\description' %in% tags) '\\description' else '\\title')] %>% unlist %>% paste(collapse='') %>% trimws(),
                   doc[tags == '\\concept'] %>% unlist %>% paste(collapse=' '),
                   'link')
    return(ret)
}

# Generate a list of the unconverted terms. This will be deprecated once ergm-terms.Rd is converted to the ergm format
.indexUnconvertedTerms <- function(name, pkg, pkg_name) {
    doc <- pkg[[name]]
    tags <- sapply(doc, attr, 'Rd_tag')
    aliases <- doc[tags == '\\alias']
    ret <- aliases[5:length(aliases)] %>% unlist %>% paste(collapse='\n')
    return(ret)
}

# Generate the dynamic index text
.generateDynamicIndex <- function() {
    entries <- c()
    for (pkg_name in .packages(TRUE)) {
        pkg <- tools::Rd_db(pkg_name)
        all_doco <- attributes(pkg)$names
        converted <- all_doco[which(endsWith(all_doco, '-ergmTerm.Rd'))]
        entries <- c(entries, sapply(converted, .indexTerm, pkg, pkg_name))
    }
                       
    entries <- c(entries,
             rep('=', 80), '\nOther Terms:\n\n',
             .indexUnconvertedTerms('ergm-terms.Rd', tools::Rd_db('ergm'), 'ergm'))
                       
    return(paste(entries, collapse=''))
}

#' An index of Ergm terms
#'
#' @name ergmTerm
#' @docType package
#' @section Term index:
#'
#' \Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex()}
NULL
