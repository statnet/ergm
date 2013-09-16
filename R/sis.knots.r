sis.knots <- function(row.sums, conjugate)
{

   v<-rep(0,length(row.sums)-1)
   k<-rep(0,length(row.sums)-1)

   for(j in 1:length(row.sums))
   {
      required.ones<-sum(row.sums[1:j]-conjugate[1:j])

      if(required.ones>0)
      {
         ##### The value of v[j] corresponds to the number    #####
         ##### number of 1's required between k[j] and k[j-1] ##### 
         v[j]<-required.ones

         ##### k=1 indicates a knot.  #####
         ##### k=0 indicates no knot. #####
         k[j]<-j
      }
   }

   ##### Eliminate values of k that are not knots, #####
   ##### and eliminate the corresponding v         #####
   v<-v[k>0]
   k<-k[k>0]
print(cbind(k,v))

   if(length(k)>1)
   {
      for(j in 2:length(k))
      {
         ##### If v[j]-v[i]>=k[j]-k[i], #####
         ##### remove knot k[i]         #####
         remove.knots<-(1:(j-1))*((v[j]-v[1:(j-1)])>=(k[j]-k[1:(j-1)]))

         remove.knots<-remove.knots[remove.knots>0]

         k[remove.knots]<-0             
      }
   }

   ##### Remove specified knots ##### 
   ##### and corresponding v    #####
   v<-v[k>0]
   k<-k[k>0]

   output<-list(knots=k, v=v)

   output
}  
