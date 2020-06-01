InitErgmConstraint.TNT<-function(lhs.nw, ...){
   if(length(list(...)))
     ergm_Init_abort(paste("TNT hint does not take arguments at this time."))
   list(dependence = FALSE, priority=10, impliedby=c("edges", "degrees", "edges", "idegrees", "odegrees", "b1degrees", "b2degrees", "idegreedist", "odegreedist", "degreedist", "b1degreedist", "b2degreedist"))
}
