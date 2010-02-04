## Uses package Rglpk

# Add library(Rglpk) call to this function.
is.inCH <- function(q, p) {
  library(Rglpk) # Needed in order to solve the linear program
	R = NROW(p)
	C=length(q)+1
	ans <- Rglpk_solve_LP(obj=c(as.vector(q),-1), mat=cbind(rbind(q,p),-1), 
						  dir=as.vector(rep("<=",R+1)), rhs=as.vector(c(1,rep(0,R))), 
						  max=TRUE, bounds=list(lower=list(ind=1:C,val=rep(-Inf,C))))
	if(ans$optimum==0){x<-TRUE}  #if the max is zero, the point q is in the CH of the points p
	else{x<-FALSE}              #if the max is strictly positive, q is not in the CH of p
	x
}




## Uses package glpk.
#
#is.inCH <- function(q, p, method="simplex") {
#  lp = lpx_create_prob()
#  lpx_set_prob_name(lp, "sample")
# lpx_set_obj_dir(lp, LPX_MAX)
# nr = NROW(p)+1
#  lpx_add_rows(lp, nr)
#  for(i in 1:(nr-1)) {
#    lpx_set_row_name(lp, i, paste("u",i,sep=""))
#    lpx_set_row_bnds(lp, i, LPX_UP, 0.0, 0.0)
#  }
#  lpx_set_row_name(lp, nr, "s")
#  lpx_set_row_bnds(lp, nr, LPX_UP, 0.0, 1.0)
#  nc = length(q)+1
#  lpx_add_cols(lp, nc)
#  for(i in 1:(nc-1)) {
#    lpx_set_col_name(lp, i, paste("x",i,sep=""))
#    lpx_set_col_bnds(lp, i, LPX_FR, 0.0, 0.0)
#    lpx_set_obj_coef(lp, i, q[i])
# }
#  lpx_set_col_name(lp, nc, "x0")
#  lpx_set_col_bnds(lp, nc, LPX_FR, 0.0, 0.0)
#  lpx_set_obj_coef(lp, nc, -1.0)
#  ia = rep(1:nr,nc)
#  ja = rep(1:nc,each=nr)
#  ar = c(as.vector(rbind(p,q)), rep(-1,nr))
#  lpx_load_matrix(lp, sum(ar!=0), ia[ar!=0], ja[ar!=0], ar[ar!=0])
#  if (method=="interior") {
#    lpx_interior(lp)
#    obj <- lpx_ipt_obj_val(lp)
#    z<-1:(nc-1)
#    for(i in 1:(nc-1))
#      z[i] <- lpx_ipt_col_prim(lp, i)
#    z0 <- lpx_ipt_col_prim(lp,nc)
#  } else {
#    lpx_simplex(lp)
#    obj <- lpx_get_obj_val(lp)
#    z<-1:(nc-1)
#    for(i in 1:(nc-1))
#      z[i] <- lpx_get_col_prim(lp, i)
#    z0 <- lpx_get_col_prim(lp,nc)
#  }
##  browser()
#  lpx_delete_prob(lp)
#  obj<=0
#}

