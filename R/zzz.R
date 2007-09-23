######################################################################
#
# zzz.R
#
# copyright (c) 2003, Mark S. Handcock, University of Washington
#                     David R. Hunter, Penn State University
#                     Carter T. Butts, University of California - Irvine
#                     Martina Morris, University of Washington
# written December 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/ergm package
#
# .First.lib is run when the package is loaded with library(ergm)
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
    cat('Type help(package="ergm") to get started.\n')
    cat('To cite, see citation("ergm")\n')
    require(network, quietly=TRUE)
}
