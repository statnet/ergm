#  File ergm/R/zzz.R
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
#
# .First.lib is run when the package is loaded.
#

.First.lib <- function(lib, pkg){
    library.dynam("ergm", pkg, lib)
    DESCpath <- file.path(system.file(package="ergm"), "DESCRIPTION")
    info <- read.dcf(DESCpath)
    cat('\nergm:', info[,"Title"], 
        '\nVersion', info[,"Version"], 'created on', info[,"Date"], '\n') 
    cat(paste("copyright (c) 2003, Mark S. Handcock, University of Washington\n",
"                    David R. Hunter, Penn State University\n",
"                    Carter T. Butts, University of California-Irvine\n",
"                    Steven M. Goodreau, University of Washington\n",
"                    Martina Morris, University of Washington\n",sep=""))
    cat('Type help(package="ergm") to get started.\n\n')
#    cat(paste('If utilization of "ergm" results in outcomes which will be published,\n',
#    'please specify the version of "ergm" you used and cite it;\n',
#    'see citation("ergm")\n'))
    cat('Based on "statnet" project software (http://statnetproject.org).\n',
        'For license and citation information see http://statnetproject.org/attribution\n',
        'or type citation("ergm").\n')
#   cat('Please cite it when you use it!\n')
#   cat('To cite, see citation("ergm")\n')
#   require(network, quietly=TRUE)
}

.Last.lib <- function(libpath){
  library.dynam.unload("ergm",libpath)
}
