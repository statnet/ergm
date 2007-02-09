#include "sisfuns.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

/********************************************
 * sissort() sorts the rows into descending *
 * order by row sum before sampling each    *
 * column.                                  *
 ********************************************/
void sissort(int *rowsums, int *order, int nrow)
{	
   int current_first,
       index_of_max;   

   double max_copy;

   for(current_first=0; current_first<(nrow-1); current_first++) 
   {
      index_of_max=get_max_range(rowsums, current_first, nrow-1);

      max_copy=rowsums[index_of_max];
      rowsums[index_of_max]=rowsums[current_first];
      rowsums[current_first]=max_copy;

      max_copy=order[index_of_max];
      order[index_of_max]=order[current_first];
      order[current_first]=max_copy;
   }
} 

/*** get_max_range() is called only by sissort() ***/
int get_max_range(int *list, int first, int last)
{
   int subscript;    
   int index_max_so_far=first;         

   for(subscript=first+1; subscript<=last; subscript++) 
      if (list[subscript]>list[index_max_so_far])
         index_max_so_far=subscript;

   return index_max_so_far;
}  



/********************************************
 * siconj() computes the conjugate sequence *
 * before sampling each column.             *
 ********************************************/
int sisconj(int *colsums, int *conjseq, int nrow, int ncol, int i)
{
   int j, k;

   for(j=i+1; j<ncol; j++){
      for(k=0; k<colsums[j]; k++){
         conjseq[k]++;
      }
   }
}


/**********************************************
 * sisknots() determines the knots and        *
 * corresponding restrictions before sampling *
 * each column.                               *
 **********************************************/
void sisknots (int *rowsums, int *conjseq, int *k, int *v, int nrow)
{
   int rowpartial=0,
       conjpartial=0;

   int i,j;

   for(i=0; i<nrow; i++)
   {
      k[i]=0;
      v[i]=0;

      rowpartial+=rowsums[i];
      conjpartial+=conjseq[i];

      if(rowpartial>conjpartial)
      {
         k[i]=i+1;
         v[i]=rowpartial-conjpartial;
      }
   }

   for(i=0; i<nrow; i++)
      if(k[i]>0)
      {
         for(j=i+1; j<nrow; j++)
            if(v[j]<=v[i])
               k[j]=0;

         for(j=0; j<i; j++)
            if((v[i]-v[j]) >= (k[i]-k[j])) 
               k[j]=0;
      }

   removeknots(k, v, nrow);
} 

/*** removeknots() is called only by sisknots() ***/
void removeknots(int *k, int *v, int nknots)
{
   int current_first,
       index_of_min;   

   double min_copy;

   for(current_first=0; current_first<(nknots-1); current_first++) 
   {
      index_of_min=get_min_range(k, current_first, nknots-1);

      min_copy=k[index_of_min];
      k[index_of_min]=k[current_first];
      k[current_first]=min_copy;

      min_copy=v[index_of_min];
      v[index_of_min]=v[current_first];
      v[current_first]=min_copy;
   }
} 

/*** get_min_range() is called only by removeknots() ***/
int get_min_range(int *list, int first, int last)
{
   int subscript;    
   int index_min_so_far=first; 

   for(subscript=first+1; subscript<=last; subscript++)
      if(list[subscript]>0)     
         if ((list[subscript] < list[index_min_so_far]) || (list[index_min_so_far]==0))
            index_min_so_far = subscript;

   return index_min_so_far;
}  



/***************************************
 * sissamp() samples 1's and 0's for a *
 * single column.                      *
 ***************************************/
void sissamp(int *k, int *v, int *rowsums, int colsum, int nrow, int ncol, int *sampled, double *prob)
{
   int i, j; 
   int nsamp, 
       total_sampled=0;
   int nknots=0;
   
   int sample[colsum];
   int ones[colsum];

   double unifprob; 
   double prob_copy;

   double weights[nrow];


   unifprob=(double) 1/(sismin(k[0],colsum)-v[0]+1);
   prob_copy=log(unifprob);

   nsamp=v[0]+(int)(unif_rand()/unifprob);

   if(k[0]==nsamp)
   {
      for(i=0; i<k[0]; i++)
      {
         ones[i]=i+1;

         if(k[i]>0)
            nknots=i+1;
      }         

      for(i=k[0]; i<nrow; i++)
         if(k[i]>0)
            nknots=i+1;

      total_sampled+=nsamp;
   }

   if(k[0]>nsamp)
   {
      for(i=0; i<k[0]; i++)
      {
         weights[i]=(double) rowsums[i]/(ncol-rowsums[i]);

         if(k[i]>0)
            nknots=i+1;        
      }

      for(i=k[0]; i<nrow; i++)
         if(k[i]>0)
            nknots=i+1;

      SampleWithoutReplace(0, k[0], weights, nsamp, sample, nrow, prob);
      prob_copy+=*prob;

      for(i=total_sampled; i<total_sampled+nsamp; i++)
         ones[i]=sample[i];

      total_sampled+=nsamp;
   }

   //Second stage: Go through remaining knots
   for(i=1; i<nknots; i++)
   {
      unifprob=(double) 1/(sismin(k[i]-k[i-1], colsum-total_sampled)-sismax(v[i]-total_sampled,0)+1);
      prob_copy+=log(unifprob);

      nsamp=v[i]-total_sampled+(int)(unif_rand()/unifprob);

      if(k[i]-k[i-1]==nsamp)
      {
         for(j=total_sampled; j< total_sampled+nsamp; j++)
            ones[j]=k[i-1]+j-total_sampled+1;

         total_sampled+=nsamp;
      }

      if(k[i]-k[i-1]>nsamp)
      {
         for(j=k[i-1]; j<k[i]; j++)
            weights[j]=(double) rowsums[j]/(ncol-rowsums[j]);

         SampleWithoutReplace(k[i-1], k[i], weights, nsamp, sample, nrow, prob);
 
         prob_copy+=*prob;

         for(j=total_sampled; j<total_sampled+nsamp; j++)
            ones[j]=sample[j-total_sampled]+k[i-1];

         total_sampled+=nsamp;
      }
   }

   for(i=0; i<colsum; i++)
      sampled[ones[i]-1]=1;

   *prob=prob_copy;
}

/*** sismin() is called only by sissamp() ***/
int sismin(int x, int y)
{
   if (x<=y)
      return x;

   return y;
}


/*** sismax() is called only by sissamp() ***/
int sismax(int x, int y)
{
   if (x>=y)
      return x;

   return y;
}

/*** SampleWithoutReplace() is called only by sissamp() ***/
void SampleWithoutReplace(int start, int stop, double *weights, int nsample, int *sample, int nrow, double *prob)
{   
   double reduced_weights[stop-start];
   double p[stop-start];
   double log_p[stop-start];

   int reduced_perm[stop-start];

   double rT, mass;
   double prob_copy=0;

   int i, j, k, n, nknots;

   for(i=start; i<stop; i++)
   {
      reduced_weights[i-start]=weights[i];
      reduced_perm[i-start]=i-start+1;
   }   
   for(i=0, nknots=stop-start-1, n=nsample; i<n; i++, nknots--, nsample--) 
   {
      for(j=0; j<nknots+1; j++)
      {
         ComputeProb(reduced_weights, nknots+1, nsample, j+1, prob);
         log_p[j]=*prob;
         p[j]=exp(*prob)/nsample;
      }

      rT=unif_rand();
      mass=0;

      for(j=0; j<nknots; j++) 
      {
         mass+=p[j];

	 if(rT<=mass)
	    break;
      }

      sample[i]=reduced_perm[j];
      prob_copy+=log_p[j];

      for(k=j; k<nknots; k++) 
      {
         reduced_weights[k]=reduced_weights[k+1];
	 reduced_perm[k]=reduced_perm[k+1];
      }
   }

   *prob=prob_copy;
}

/*** ComputeProb() is called only by SampleWithoutReplace() ***/
void ComputeProb(double *weights, int nweights, int nsample, int weight_locat, double *prob)
{
   int i;

   double numerator,
          denominator;

   double reduced_weights[nweights-1];

   for(i=0; i<weight_locat-1; i++)
      reduced_weights[i]=weights[i];

   for(i=weight_locat-1; i<nweights-1; i++)
      reduced_weights[i]=weights[i+1];

   RecursiveProb(reduced_weights, nweights-1, nsample-1, prob);
   numerator=*prob;

   RecursiveProb(weights, nweights, nsample, prob);
   denominator=*prob;

   *prob=log(weights[weight_locat-1])+numerator-denominator;

   if(*prob>0)
     *prob=0;
}

/*** RecursiveProb() is called only by ComputeProb() ***/
void RecursiveProb(double *weights, int nweights, int nsample, double *prob)
{
  int k, j, minkm;
  double temp0, temp;
  double *current, *previous;

  current  = (double*) malloc(sizeof(double) * (nsample+1));
  previous = (double*) malloc(sizeof(double) * (nsample+1));

  *prob=0.0;
  if (nsample==0)
  {
     free(current);
     free(previous);

     return;
  }
  if (nsample==1)
  {
    for(k=0; k<nweights; k++)
    {
       *prob += weights[k];
    }

    *prob=log(*prob);

     free(current);
     free(previous);

    return;
  }
  current[0] = 0.0;
  current[1] = log(weights[0]);
  for(k=1; k<nweights; k++){
    for(j=0; j<nsample+1; j++){
     previous[j] = current[j];
    }
    minkm = k+1;
    if(minkm > nsample){minkm=nsample;}
    for(j=1; j<minkm; j++){
      temp0 = log(weights[k])+previous[j-1];
      temp = previous[j];
      if(temp0 > temp){ 
        current[j] = log(1.0 + exp(temp-temp0)) + temp0;
      }else{
        current[j] = log(exp(temp0-temp) + 1.0) + temp;
      }
    }
    if(k < nsample){
     current[minkm] = log(weights[k])+previous[minkm-1];
    }else{
      temp0 = log(weights[k])+previous[minkm-1];
      temp = previous[minkm];
      if(temp0 > temp){ 
        current[minkm] = log(1.0 + exp(temp-temp0)) + temp0;
      }else{
        current[minkm] = log(exp(temp0-temp) + 1.0) + temp;
      }
    }
  }
  *prob = current[nsample];
  free(current);
  free(previous);
}



/********************************************
 * sisufun() updates row sums and reorders  *
 * the sequence of 1's and 0's according to *
 * the switch that was made by sissort()    *
 ********************************************/
void sisufun(int *rowsums, int *sequence, int *order, int nrow)
{
   int i;

   int newseq[nrow];
   int newrowsums[nrow];
   int neworder[nrow];

   for(i=0; i<nrow; i++)
   {
      rowsums[i]-=sequence[i];
      newseq[order[i]-1]=sequence[i];
      newrowsums[order[i]-1]=rowsums[i];
      neworder[order[i]-1]=order[i];
   }

   for(i=0; i<nrow; i++)
   {
      sequence[i]=newseq[i];
      rowsums[i]=newrowsums[i];
      order[i]=neworder[i];
   }
}

void sissim(int *rowsums, int *colsums, int *newmat,
	   int *nrow, int *ncol, 
	   double *samplesize, 
	   int *verb, double *prob, double *probvec)
{
   int i, j, isamp;
   int nrowd, ncold;
   int conjugatevec[*nrow];
   int new_sample[*nrow];
   int v[*nrow];
   int k[*nrow];
   int rowsumd[*nrow];
   int colsumd[*ncol];
   int order[*nrow];
   long int nsamples;
   double prob_copy;

   nsamples=(long int)*samplesize;
   
   GetRNGstate();  /* R function enabling uniform RNG */

   for(isamp=0; isamp<nsamples; isamp++)
   {
     nrowd=*nrow;
     ncold=*ncol;

     /*** Initialize row sum and order vectors ***/
     for(i=0; i<*nrow; i++)
     {
        order[i]=i+1;
        rowsumd[i]=rowsums[i];
     }

     for(i=0; i<*ncol; i++)
        colsumd[i]=colsums[i];

     prob_copy=0;

     for(i=0; i<*ncol; i++)
     {    
	if(colsumd[i]>0)
        {
	   if(colsumd[i]<nrowd)
           {
              for(j=0; j<nrowd; j++)
              {
                 new_sample[j]=0;
                 conjugatevec[j]=0;
              }

              /*** 1. Sort rows ***/
              sissort(rowsumd, order, nrowd);

              /*** 2. Compute conjugate sequence ***/
              sisconj(colsumd, conjugatevec, nrowd, ncold, i);

              /*** 3. Compute Knots ***/
              sisknots(rowsumd, conjugatevec, k, v, nrowd);

              /*** 4. Sample ***/
              sissamp(k, v, rowsumd, colsumd[i], nrowd, ncold-i, new_sample, prob);

              prob_copy+=*prob;

              /*** 5. Update ***/ 
              sisufun(rowsumd, new_sample, order, nrowd);

              for(j=0; j<nrowd; j++)
              {
                 newmat[i*nrowd+j]=new_sample[j];
              }

           }

           if(colsumd[i]==nrowd)
           {
              for(j=0; j<nrowd; j++)
              {
		 new_sample[j]=1;
                 newmat[i*nrowd+j]=1;
              }

              sisufun(rowsumd, new_sample, order, nrowd);
           }
        }

	if(colsumd[i]==0)
        {
           for(j=0; j<nrowd; j++)
           {
              new_sample[j]=0;
              newmat[i*nrowd+j]=0;
           }
        }
     }

     probvec[isamp]=prob_copy;

   }
   PutRNGstate();  /* Disable RNG before returning */

}
