InitErgmConstraint.TNT<-function(lhs.nw, ...){
   if(length(list(...)))
     ergm_Init_abort(paste("TNT hint does not take arguments at this time."))
   list(dependence = FALSE, priority=10, impliedby=c("edges", "degrees", "edges", "idegrees", "odegrees", "b1degrees", "b2degrees", "idegreedist", "odegreedist", "degreedist", "b1degreedist", "b2degreedist"))
}

InitErgmConstraint.Strat <- function(lhs.nw, attr=NULL, pmat=NULL, empirical=NULL, ...) {
   if(...length()) ergm_Init_abort(paste0("Unrecognised argument(s) ", paste.and(names(list(...)), oq="'", cq="'"),".")) 

   list(dependence = FALSE, priority=10, attr=attr, pmat=pmat, empirical=empirical)
}
