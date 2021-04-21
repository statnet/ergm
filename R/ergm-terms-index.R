#  File R/ergm-terms-index.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2021 Statnet Commons
#######################################################################

# Return the index entry for a single term in the new format
.parseTerm <- function(name, pkg, pkg_name) {
    doc <- pkg[[name]]
    tags <- sapply(doc, attr, 'Rd_tag')

	raw_usage <- doc[tags == '\\usage'] %>% 
        unlist %>% 
        paste(collapse="") %>%
        gsub("\n *# *"," ", .) %>%
        strsplit("\n") %>%
        .[[1]] %>%
        trimws()
	usages <- list()
	for (usage_line in regmatches(raw_usage, regexec("^(binary|valued): *(.+)$", raw_usage))) {
		usages[[length(usages) + 1]] <- list(type=usage_line[2], value=usage_line[3])
	}

    ret <- list(
		name=substr(name, 1, nchar(name) - 3),
		alias=doc[tags == '\\alias'] %>% unlist,
		package=pkg_name,
		usages=usages,
		title=doc[tags == '\\title'] %>% unlist,
		description=doc[tags == '\\description'] %>% unlist %>% paste(collapse='') %>% trimws(),
		concepts=doc[tags == '\\concept'] %>% unlist)
    return(ret)
}

.ergmTermIndexStore <- list()
.ergmTermIndexParsedPackages <- c()
# Crawl all loaded packages for terms, parse them once and store in the singleton store
.crawlLoadedLibraries <- function() {
	loaded_libraries <- .packages(TRUE)
	for (pkg_name in loaded_libraries[!loaded_libraries %in% .ergmTermIndexParsedPackages]) {
        pkg <- tools::Rd_db(pkg_name)
        all_doco <- attributes(pkg)$names
        converted <- all_doco[which(endsWith(all_doco, '-ergmTerm.Rd'))]
		if (length(converted) > 0) {
			assign('.ergmTermIndexStore', c(.ergmTermIndexStore, lapply(converted, .parseTerm, pkg, pkg_name)), envir=.GlobalEnv)
		}
	}

	assign('.ergmTermIndexParsedPackages', loaded_libraries, envir=.GlobalEnv) 
}

# Generate the index entry for a single term
.genTermEntry <- function(term) {
    ret <- sprintf('%s (%s)\n    %s\n\n    Keywords: %s\n    %s\n\n',
                   term$name,
                   term$package,
				   if (!is.null(term$description)) term$description else term$title,
				   paste(term$concept, collapse=' '),
                   'link')
    return(ret)
}

# Generate the dynamic index text
.generateDynamicIndex <- function() {
	# TODO: Better to add this to the hooks so that this runs once when the package is loaded, then as new packages are loaded
	.crawlLoadedLibraries()

    return(paste(sapply(.ergmTermIndexStore, .genTermEntry), collapse=''))
}

#' An index of Ergm terms
#'
#' @name ergmTerm
#' @docType package
#' @section Term index:
#'
#' \Sexpr[results=rd,stage=render]{ergm:::.generateDynamicIndex()}
NULL
