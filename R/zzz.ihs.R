######################################################################
#
# zzz.R
#
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# 
# For license information see http://statnetproject.org/license
#
# We have invested a lot of time and effort in creating 'statnet',
# for use by other researchers. We require that the attributions
# in the software are retained (even if only pieces of it are used),
# and that there is attribution when the package is loaded (e.g., via
# "library" or "require"). This is to stop
# "rebadging" of the software. 
#
# Cite us!
#
# To cite see http://statnetproject.org/cite
#
# .First.lib is run when the package is loaded.
#
######################################################################

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
    cat('Based on "statnet" project software (http://statnetproject.org),\n',
    '  for license and citation information see citation("ergm")\n')
    cat('Please cite it when you use it!\n')
#   cat('To cite, see citation("ergm")\n')
#   require(network, quietly=TRUE)
}

if(exists("ergm.ihs")) { # Look here for a list of functions that exist in two
  # distinct versions now, one for CRAN and one not.  This makes debugging 
  # twice as much work for these functions.  :(
  ergm.mple <- ergm.mple.ihs
  ergm.initialfit <- ergm.initialfit.ihs
  ergm.mainfitloop <- ergm.mainfitloop.ihs
  ergm<-ergm.ihs
}
