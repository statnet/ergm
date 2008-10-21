#  File ergm/R/control.simulate.R
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
control.simulate<-control.simulate.formula<-function(prop.weights="default",
                                                     prop.args=NULL,
                                                     drop=FALSE,
                                                     summarizestats=FALSE,
                                                     maxchanges=1000000,
                                                     packagenames="ergm",
                                                     parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}

control.simulate.ergm<-function(prop.weights=NULL,
                                prop.args=NULL,
                                drop=FALSE,
                                summarizestats=FALSE,
                                maxchanges=1000000,
                                packagenames="ergm",
                                parallel=0){
  control<-list()
  for(arg in names(formals(sys.function())))
    control[[arg]]<-get(arg)
  control
}
