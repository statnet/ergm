#  File ergm/R/statnetShowDoc.R
#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
# Copyright 2003 Mark S. Handcock, University of Washington
#                David R. Hunter, Penn State University
#                Carter T. Butts, University of California - Irvine
#                Steven M. Goodreau, University of Washington
#                Martina Morris, University of Washington
# Copyright 2007 The statnet Development Team
######################################################################
statnetShowDoc <- function (what="LICENSE", type = c("html", "txt"),
                            package="statnet", title="License Information") 
{
    type <- match.arg(type)
    html_viewer <- function(path) {
        if (.Platform$OS.type == "windows"){
          shell.exec(chartr("/", "\\", path))
        }else{
          htmlb <- statnetbrowseURL(paste("file://", URLencode(path), sep = ""))
          if(htmlb!=0){
           htmlb <- try(system(paste("lynx file://", URLencode(path), sep = "")))
          }
        }
    }
    paste. <- function(x, ext) paste(x, ext, sep = ".")
    if (length(what) != 1 || !is.character(what)) {
        message("   statnetShowDoc() should be used with a character string argument specifying\n   a documentation file")
        return(invisible())
    }
    pkgpath <- .find.package(package)
    path <- file.path(pkgpath, what)
    if (file.exists(path)) {
      if (type == "txt") {
        file.show(path, title=title)
        return(invisible(path))
      }
      if (file.exists(path)) {
        html_viewer(path)
        return(invisible(path))
      }
    }
    stop(gettextf("no documentation for '%s' found in package '%s'", 
         what, package), domain = NA)
}
statnetbrowseURL <- function (url, browser = getOption("browser")) 
{
    shQuote <- function(string) paste("\"", gsub("\\$", "\\\\$", 
        string), "\"", sep = "")
    if (!is.character(url) || !(length(url) == 1) || !nzchar(url)) 
        stop("'url' must be a non-empty character string")
    if (!is.character(browser) || !(length(browser) == 1) || 
        !nzchar(browser)) 
        stop("'browser' must be a non-empty character string")
    if (.Platform$GUI == "AQUA" || length(grep("^(localhost|):", 
        Sys.getenv("DISPLAY"))) > 0) 
        isLocal <- TRUE
    else isLocal <- FALSE
    quotedUrl <- shQuote(url)
    remoteCmd <- if (isLocal) 
        switch(basename(browser), "gnome-moz-remote" = , open = quotedUrl, 
            galeon = paste("-x", quotedUrl), kfmclient = paste("openURL", 
                quotedUrl), netscape = , mozilla = , opera = , 
            firefox = {
                paste("-remote \"openURL(", gsub("([,)$])", "%\\1", 
                  url), ")\"", sep = "")
            }, quotedUrl)
    else quotedUrl
    try(system(paste(browser, remoteCmd, "> /dev/null 2>&1 ||", 
                     browser, quotedUrl)))
}
