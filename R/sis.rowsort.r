sis.rowsort <- function(row.sums, nrow)
{
   ordered.row.sums<-row.sums[order(row.sums, decreasing=T)]
   ord<-order(row.sums, decreasing=T)
   output<-list(sums=ordered.row.sums, order=ord)

   output
}  
