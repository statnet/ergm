# Must remove the .dep from the end of this filename to compile the CRAN version!

######################################################################
# File name: zzz.R
######################################################################
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# 
# For license and attribution information see
#    http://statnetproject.org/attribution
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). 
######################################################################
#
# .First.lib is run when the package is loaded.
#

.First.lib <- function(lib, pkg){
    library.dynam("ergm", pkg, lib)
    ehelp <- library(help="ergm",lib.loc=NULL,character.only=TRUE)$info[[1]]
    cat(paste(substring(ehelp[4],first=16),"\n",
              "Version ",substring(ehelp[2],first=16),
              " created on ",
               substring(ehelp[3],first=16),".\n", sep=""))
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
