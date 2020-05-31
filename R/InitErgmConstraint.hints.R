InitErgmConstraint.TNT<-function(lhs.nw, ...){
   if(length(list(...)))
     ergm_Init_abort(paste("TNT hint does not take arguments at this time."))
   list(dependence = FALSE, priority=1)
}
