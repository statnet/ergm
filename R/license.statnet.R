###############################################################################
# The <license.statnet> function prints the licensing information for any of
# the statnet packages
#
# --PARAMETERS--
#   package: the name of one of the statnet packages; default="statnet"
#
# --RETURNED--
#   NULL
#
###############################################################################

license.statnet <- function (package="statnet") {

statnet.packages <- c("statnet","ergm","degreenet","latentnet","network","networksis","sna")
package <- pmatch(package, statnet.packages)
if(is.na(package)){
   stop(gettextf("The package name should be one of %s", paste(dQuote(statnet.packages), 
            collapse = ", ")), domain = NA)
}
package <- statnet.packages[package]

#if(package %in% c("ergm","degreenet","networksis")
if(package %in% c("network","sna")){
cat(
" --------------------------------------------------\n",
paste("License for the 'statnet' component package '",package,"'\n",sep=""),
"--------------------------------------------------\n",
"\n",
" This software is distributed under the GPL-2 license.  It is free\n",
"and open source.\n")
}else{
if(package=="statnet"){
cat(
" --------------------------------------------------\n",
"License for the 'statnet' package\n",
"--------------------------------------------------\n",
"\n")}else{
cat(
" --------------------------------------------------\n",
paste("License for the 'statnet' component package '",package,"'\n",sep=""),
"--------------------------------------------------\n",
"\n")
}
cat(
" This software is distributed under the GPL-3 license.  It is free,\n",
"open source, and has the following attribution requirements\n",
"(GPL Section 7):\n",
"\n",
paste("(a) you agree to retain in '",package,"' and any modifications to\n",sep=""),
paste("    '",package,"' the copyright, author attribution and URL\n",sep=""),
"    information as provided at http://statnetproject.org/attribution\n",
paste("(b) you agree that '",package,"' and any modifications to '",package,"' will,\n",sep=""),
"    when used, display the attribution:\n",
"\n",
"      Based on 'statnet' project software (http://statnetproject.org).\n",
"      For license and citation information see\n",
"      http://statnetproject.org/attribution\n\n",
"--------------------------------------------------\n",
"\n",
"What does this mean?\n",
"====================\n",
"\n",
paste("If you are modifying '",package,"' or adopting any source code from\n",sep=""),
paste("'",package,"' for use in another application, you must ensure that the\n",sep=""),
"copyright and attributions mentioned in the license above appear\n",
"in the code of your modified version or application.  These\n",
"attributions must also appear when the package is loaded\n",
"(e.g., via 'library' or 'require').",
"\n\n\n",
"Enjoy!\n",
"\n",
"Copyright 2003 Mark S. Handcock, University of Washington\n",
"               David R. Hunter, Penn State University\n",
"               Carter T. Butts, University of California-Irvine\n",
"               Steven M. Goodreau, University of Washington\n",
"               Martina Morris, University of Washington\n",
"\n",
"Copyright 2010 The statnet Development Team\n")
}
}
