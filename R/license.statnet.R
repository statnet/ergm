license.statnet <- function (package="ergm") {
#if(.Platform$OS.type="unix" & .Platform$GUI!="AQUA"){
# cat("This is the license.  For citation information see the R console window.
#}

#if (getRversion() >= 2.5) {
#  RShowDoc("LICENSE",package="ergm")
#} else {
#  cat("Full license avaialable in LICENSE file in 'ergm' package subdirectory.\n")
#}

#if (.Platform$OS.type == "windows"){
# shell.exec(chartr("/", "\\", path))
#}else{
#htmlb <- try(statnetbrowseURL("http://statnetproject.org/attribution/ergm.shtml"))
# if(htmlb!=0){
# htmlb <- try(statnetbrowseURL("http://statnetproject.org/attribution/ergm.shtml", browser="lynx"))
# }
#}

#statnetShowDoc("LICENSE",package=package)

cat(
" --------------------------------------------------\n",
"License for the 'statnet' component package 'ergm'\n",
"--------------------------------------------------\n",
"\n",
"This software is distributed under the GPL-3 license.  It is free,\n",
"open source, and has the following attribution requirements\n",
"(GPL Section 7):\n",
"\n",
"(a) you agree to retain in 'ergm' and any modifications to\n",
"    'ergm' the copyright, author attribution and URL\n",
"    information as provided at http://statnetproject.org/attribution\n",
"(b) you agree that 'ergm' and any modifications to 'ergm' will,\n",
"    when used, display the attribution:\n",
"\n",
"      Based on 'statnet' project software (http://statnetproject.org).\n",
"      For license and citation information see\n",
"      http://statnetproject.org/attribution\n\n",
"(c) you agree that 'ergm' and any modifications to 'ergm' will display\n",
"    the citation information, as provided in the original CITATION\n",
"    file, when the 'citation' function in executed with the name\n",
"    of the package given.\n",
"(d) you agree to retain in the documentation contained within\n",
"    this software, and any modifications of it, the author\n",
"    attribution (GPL Section 7, clause c).\n",
"--------------------------------------------------\n",
"\n",
"What does this mean?\n",
"====================\n",
"\n",
"If you are modifying 'ergm' or adopting any source code from\n",
"'ergm' for use in another application, you must ensure that the\n",
"copyright and attributions mentioned in the license above appear\n",
"in the code of your modified version or application.  These\n",
"attributions must also appear when the package is loaded\n",
"(e.g., via 'library' or 'require').",
"\n\n\n",
"Enjoy!\n",
"\n",
"Mark S. Handcock, University of Washington\n",
"David R. Hunter, Penn State University\n",
"Carter T. Butts, University of California-Irvine\n",
"Steven M. Goodreau, University of Washington\n",
"Martina Morris, University of Washington\n",
"\n",
"The 'statnet' development team\n",
"\n",
"Copyright 2007\n")
#   cat("\nThis software is distributed under the terms of the GNU General\n")
#   cat("Public License Version 2, June 1991.  The terms of this license\n")
#   cat("are in a file called COPYING which you should have received with\n")
#   cat("this software and which can be displayed by RShowDoc(\"COPYING\").\n")
#   cat("\n")
#   cat("If you have not received a copy of this file, you can obtain one\n")
#   cat("at http://www.R-project.org/licenses/.\n")
#   cat("\n")
#   cat("A small number of files (the API header files listed in\n")
#   cat("R_DOC_DIR/COPYRIGHTS) are distributed under the\n")
#   cat("Lesser GNU General Public LIcense version 2.1.\n")
#   cat("This can be displayed by RShowDoc(\"COPYING.LIB\"),\n")
#   cat("or obtained at the URI given.\n")
#   cat("\n")
#   cat("'Share and Enjoy.'\n\n")
}
license.statnet1 <- function (package="ergm", 
         title="License information for the 'statnet' component package ergm") {
  pkgpath <- try(.find.package(package))
  if(!inherits(pkgpath,"try-error")){
   path <- file.path(pkgpath, "LICENSE")
   if (file.exists(path)) {
    file.show(path, title=title)
    return(invisible(path))
   }
  }
  stop(gettextf("\nno license found in package '%s'.\n Based on 'statnet' project software (http://statnetproject.org).\n For license and citation information see http://statnetproject.org/attribution",
        package), domain = NA)
}
