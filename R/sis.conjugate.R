sis.conjugate <- function(col.sums, nrow)
{
   col.sum.mat<-matrix(rep(col.sums,nrow), nrow=length(col.sums), ncol=nrow)
   j.mat<-t(matrix(rep(1:nrow,length(col.sums)), nrow=nrow, ncol=length(col.sums)))

   conjugate<-apply((col.sum.mat-j.mat)>=0,2,sum)

   conjugate.mat<-t(matrix(rep(conjugate,length(col.sums)), nrow=nrow, ncol=length(col.sums)))

   for(i in 2:nrow(conjugate.mat))
   {
      conjugate.mat[i,]<-conjugate.mat[i-1,]
      conjugate.mat[i,(1:col.sums[i-1])]<-conjugate.mat[(i-1),(1:col.sums[i-1])]-1
   }

   return(conjugate.mat)
}  
