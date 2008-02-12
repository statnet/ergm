#include "changestats.h"


/********************  changestats:  A    ***********/
/*****************                       
 changestat: d_absdiff
*****************/
CHANGESTAT_FN(d_absdiff) { 
  double change;
  Vertex h, t;
  int i;

  CHANGE_STAT[0] = 0.0;
  for (i=0; i<ntoggles; i++) {
    h = heads[i]; 
    t = tails[i];
    change = fabs(INPUT_ATTRIB[h-1] - INPUT_ATTRIB[t-1]);
    CHANGE_STAT[0] += IS_OUTEDGE(h,t) ? -change : change;
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_absdiffcat
*****************/
CHANGESTAT_FN(d_absdiffcat) { 
  double change, absdiff, NAsubstitute, hval, tval;
  Vertex h, t, ninputs;
  int i, j;
  
  ninputs = N_INPUT_PARAMS - N_NODES;
  NAsubstitute = INPUT_PARAM[ninputs-1];
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    change = IS_OUTEDGE(h=heads[i], t=tails[i]) ? -1.0 : 1.0;
    hval = INPUT_ATTRIB[h-1];
    tval = INPUT_ATTRIB[t-1];
    if (hval == NAsubstitute ||  tval == NAsubstitute) absdiff = NAsubstitute;
    else absdiff = fabs(hval - tval);
	  if (absdiff>0) {
      for (j=0; j<ninputs; j++) {
        CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? change : 0.0;
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_altkstar
*****************/
CHANGESTAT_FN(d_altkstar) { 
  int i, isedge;
  double lambda, oneexpl, change;
  Vertex h, t, hd, td=0;
  
  change = 0.0;
  lambda = INPUT_PARAM[0];
  oneexpl = 1.0-1.0/lambda;
  
  for (i=0; i<ntoggles; i++) {
    isedge = IS_OUTEDGE(h=heads[i], t=tails[i]);
    hd = OUT_DEG[h] + IN_DEG[h] - isedge;
    td = OUT_DEG[t] + IN_DEG[t] - isedge;
    if(hd!=0){
      change += (1-2*isedge)*(1.0-pow(oneexpl,(double)hd));
    }
    if(td!=0){
      change += (1-2*isedge)*(1.0-pow(oneexpl,(double)td));
    }
    
    if (i+1 < ntoggles)
      TOGGLE(h, t);  /* Toggle this edge if more to come */
  }
  CHANGE_STAT[0] = change*lambda;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_asymmetric
*****************/
CHANGESTAT_FN(d_asymmetric) { 
  Vertex h, t;
  int i, edgeflag, refedgeflag;
  
  CHANGE_STAT[0] = 0.0;
  
  for (i=0; i<ntoggles; i++) {
    edgeflag = IS_OUTEDGE(h=heads[i], t=tails[i]);
    refedgeflag = IS_OUTEDGE(t,h);
    CHANGE_STAT[0] += ((edgeflag==refedgeflag) ? 1.0 : -1.0); 

    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  B    ***********/
/*****************
 changestat: d_b1concurrent
*****************/
CHANGESTAT_FN(d_b1concurrent) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  Vertex b1, b2, actdeg, *od;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  od=nwp->outdegree;
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    actdeg = od[b1];
    *(mtp->dstats) += (actdeg + echange > 1) - (actdeg > 1);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b1concurrent_by_attr
*****************/
CHANGESTAT_FN(d_b1concurrent_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, b1attr;
  Vertex b1, b2, b1deg, *od;
  TreeNode *oe;

  oe=nwp->outedges;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b1deg = od[b1];
    b1attr = mtp->inputparams[N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b1attr == mtp->inputparams[j]) { /* we have attr match */
        mtp->dstats[j] += (b1deg + echange > 1) - (b1deg > 1);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b1factor
*****************/
CHANGESTAT_FN(d_b1factor) { 
  double s, factorval;
  Vertex h, t;
  int i, j;
  
  for (i=0; i < N_CHANGE_STATS; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[h-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b1degree
*****************/
CHANGESTAT_FN(d_b1degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, j, echange;
  Vertex b1, b2, actdeg, d, *od;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    actdeg = od[b1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(mtp->inputparams[j]);
      mtp->dstats[j] += (actdeg + echange == d) - (actdeg == d);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b1degree_by_attr
*****************/
CHANGESTAT_FN(d_b1degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, b1attr;
  Vertex b1, b2, b1deg, d, *od;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b1deg = od[b1];
    b1attr = mtp->inputparams[2*N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b1attr == mtp->inputparams[2*j+1]) { /* we have attr match */
        d = (Vertex)mtp->inputparams[2*j];
        mtp->dstats[j] += (b1deg + echange == d) - (b1deg == d);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_balance
*****************/
CHANGESTAT_FN(d_balance) { 
  int i, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex node3, h, t;

  *(mtp->dstats) = 0.0;

  if (nwp->directed_flag) {
/* directed version */
   for (i=0; i<ntoggles; i++) 
    {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++)
             { 
             if (node3 != h && node3 != t)
               {   
               a = (EdgetreeSearch(t, h, nwp->outedges) != 0); 
               b = (EdgetreeSearch(t, node3, nwp->outedges) != 0);
               c = (EdgetreeSearch(node3, t, nwp->outedges) != 0);
               d = (EdgetreeSearch(node3, h, nwp->outedges) != 0);
               e = (EdgetreeSearch(h, node3, nwp->outedges) != 0);
               edgecount = (a + b + c + d + e);

               switch(edgecount)
               {  
               case 0:   /* 012 */
                 ++t012;

               case 1:   /* 021C, 021U, 021D, 102 */
                 {
                     if ((b == 1) || (d == 1))
                       ++t021C;
                     if (c == 1)
                       ++t021U;
                     if (e == 1)
                       ++t021D;
                     if (a == 1)
                       ++t102;
                 }
                 break;
          
               case 2:   /* 030C, 030T, 111U, 111D */
                 {
                     if ((b + d) == 2)       
                       ++t030C;
                     if (((b + e) == 2) || ((c + d) == 2) || ((c + e) == 2))
                       ++t030T;
                     if (((a + b) == 2) || ((a + e) == 2) || ((d + e) == 2))
                       ++t111U;
                     if (((a + c) == 2) || ((a + d) == 2) || ((b + c) == 2))
                       ++t111D;
                 }
                 break;
            
               case 3:   /* 120C, 120U, 120D, 201 */
                 {
                     if (a == 1)
                         {
                            if (((b + d) == 2) || ((c + e) == 2))
                              ++t120C;
                            if ((b + e) == 2)
                              ++t120U;
                            if ((c + d) == 2)
                              ++t120D;
                            if (((b + c) == 2) || ((d + e) == 2))
                              ++t201;
                         }
                     else 
                         {
                            if (b == 1)
                                {
                                    if (((c + d) == 2) || ((d + e) == 2))
                                      ++t120C;
                                    if ((c + e) == 2)
                                      ++t120D;
                                } 
                            else 
                                {
                                    ++t120U;
                                }
                         } 
                 }
                 break;
             
               case 4:   /* 210 */
                 ++t210;
                 break;
             
               case 5:   /* 300 */            
                 ++t300;
                 break;
               }

               switch(edgecount)
               {  
               case 1:   /* 102, 021D, 021U, 021C */
                 --t012;
                 break;
           
               case 2:   /* 030C, 030T, 111U, 111D */ 
                 {
                     if (((a + c) == 2) || ((a + e) == 2) || ((b + d) == 2) || 
                        ((c + e) == 2)) 
                       --t021C;
                     if (((a + d) == 2) || ((b + e) == 2))
                       --t021U;
                     if (((a + b) == 2) || ((c + d) == 2))
                       --t021D;
                     if (((b + c) == 2) || ((d + e) == 2))
                       --t102;
                 } 
                 break;

               case 3:   /* 201, 120D, 120U, 120C */
                 {
                     if (a == 1)
	    	           {
	    	           if ((c + e) == 2)       
                             --t030C;
                           if (((c + d) == 2) || ((b + e) == 2) || ((b + d) == 2))
                             --t030T;
                           if ((b + c) == 2)
                             --t111U;
                           if ((d + e) == 2)
	                      --t111D;
	 	               }
 		             else
 		               {
 		                   if (b == 1)
 		 	                 {
 		 	                   if ((c + d) == 2)
 			                     --t111U;
 			                   if (((c + e) == 2) || ((d + e) == 2))
			                     --t111D;
 		 	                 }
    	    	                   else
		    	                 --t111U;
    		               }
                   }
                   break;
          
               case 4:   /* 210 */
                 {
                     if (a == 1)
                         {
                            if (((b + c + e) == 3) || ((c + d + e) == 3))
                              --t120C;
                            if ((b + c + d) == 3)
                              --t120U;
                            if ((b + d + e) == 3)
                              --t120D;
                         }
                     else 
                         {
                            if ((b + c + d + e) == 4)
                              --t201;
                         } 
                 }
                 break;
              
               case 5:   /* 300 */            
                 --t210;
                 break;
               }
               }
             }    /* ******  move to next node3 ******** */
        }
        else 
          t012 = t012 + (N_NODES - 2);  

/*        t003 = (t300+t210+t120C+t120U+t120D+t201+t030C+t030T); 
        t003 = t003+(t111U+t111D+t021C+t021U+t021D+t102+t012); */
	b = t102 + t300; 
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;

      if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */ 
    }
    }else{
/*  undirected */
    for (i=0; i<ntoggles; i++) 
     {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++)
             { 
             if (node3 != h && node3 != t)
               {   
               a = (EdgetreeSearch(MIN(node3,t), MAX(node3,t), nwp->outedges) != 0);
               b = (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0);
               edgecount = (a + b);

               switch(edgecount)
               {  
               case 0:   /* 012 */
		{
                 ++t102;
                 --t012;
		}
                break;

               case 1:   /* 021C, 021U, 021D, 102 */
		{
                 ++t201;
                 --t102;
		}
                break;
          
               case 2:   /* 030C, 030T, 111U, 111D */
		{
                 ++t300;
                 --t201;
		}
                break;
            
               }
	       }

             }    /* ******  move to next node3 ******** */
        }
        else 
          t102 = t102 + (N_NODES - 2);  

        t003 = (t102+t201+t300);
	b = t102 + t300; 
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;
  
     if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */ 
     } /* i loop */
    }
    
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundeddegree
*****************/
CHANGESTAT_FN(d_boundeddegree) { 
  int i, j, echange;
  Vertex h, t, hd, td=0, deg, *id, *od;
  TreeNode *oe=nwp->outedges;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
      hd = od[h] + id[h];
      td = od[t] + id[t];
      for(j = 0; j+1 < nstats; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
          mtp->dstats[j] += (td + echange == deg) - (td == deg);
	}
      mtp->dstats[nstats-1] += (hd + echange >= bound) - (hd >= bound);
      mtp->dstats[nstats-1] += (td + echange >= bound) - (td >= bound);
      
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundedidegree
*****************/
CHANGESTAT_FN(d_boundedidegree) { 
  int i, j, echange;
  Vertex h, t, hd=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = nwp->indegree[h];
      for(j = 0; j < N_CHANGE_STATS; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
	}
      mtp->dstats[nstats-1] += (hd + echange >= bound) - (hd >= bound);
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 my_choose:
 Simple routine to return simple binomial coefficients quickly, 
 avoiding costly call to choose() function.  Note:  my_choose is
 usually not called directly; use CHOOSE macro instead.
*****************/
#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 
double my_choose(double n, int r) {
  const double recip_factorial[21] = {1.0, 1.0, 0.5,
	 1.66666666666667e-01, 4.16666666666667e-02, 8.33333333333333e-03,
         1.38888888888889e-03, 1.98412698412698e-04, 2.48015873015873e-05,
         2.75573192239859e-06, 2.75573192239859e-07, 2.50521083854417e-08,
         2.08767569878681e-09, 1.60590438368216e-10, 1.14707455977297e-11,
         7.64716373181982e-13, 4.77947733238739e-14, 2.81145725434552e-15,
         1.56192069685862e-16, 8.22063524662433e-18, 4.11031762331216e-19};
  double ans;

  if (r>20)
    return choose (n, (double)r); /* Use complicated function for large r */
  for(ans=recip_factorial[r]; r>0; r--)
    ans*=(n--);
  return ans;
}

/*****************
 changestat: d_boundedistar
*****************/
CHANGESTAT_FN(d_boundedistar) { 
  double change, tod;
  double newtod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex h, t;
  
  for (i=0; i < p; i++) 
    mtp->dstats[i] = 0.0;

  for (i=0; i < ntoggles; i++)
    {
      /* is there an edge for this toggle */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);

      tod = nwp->indegree[t];
      newtod = tod + (edgeflag ? -1 : 1);
      
      for(j=0; j < p; j++) 
	{
	  k =  ((int)mtp->inputparams[j]);
	  bound = (int)mtp->inputparams[j+p];
	  change = MIN(bound,CHOOSE(newtod, k))-MIN(bound,CHOOSE(tod, k));
            	  
	  mtp->dstats[j] += change;
      	}
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundedkstar
*****************/
CHANGESTAT_FN(d_boundedkstar) { 
  double change, hod, tod;
  double newhod, newtod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex h, t;
  
  for (i=0; i < p; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      /* is there an edge for this toggle */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);

      hod = nwp->outdegree[h] + nwp->indegree[h];
      newhod = hod + (edgeflag ? -1 : 1);
      tod = nwp->outdegree[t] + nwp->indegree[t];
      newtod = tod + (edgeflag ? -1 : 1);
      
      for(j=0; j < p; j++) 
	{
	  k =  ((int)mtp->inputparams[j]);
	  bound = (int)mtp->inputparams[j+p];
	  change = (MIN(bound,CHOOSE(newhod, k))-MIN(bound,CHOOSE(hod, k))) 
	          +(MIN(bound,CHOOSE(newtod, k))-MIN(bound,CHOOSE(tod, k)));
            	  
	  mtp->dstats[j] += change; /* (edgeflag ? - change : change); */
      	}
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundedodegree
*****************/
CHANGESTAT_FN(d_boundedodegree) { 
  int i, j, echange;
  Vertex h, t, hd=0, deg;
  int nstats = (int)N_CHANGE_STATS;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = nwp->outdegree[h];
      for(j = 0; j < N_CHANGE_STATS; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
	}
      mtp->dstats[nstats-1] += (hd + echange >= bound) - (hd >= bound);
      
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundedostar
*****************/
CHANGESTAT_FN(d_boundedostar) { 
  double change, hod;
  double newhod;
  int edgeflag, i, j, k, bound;
  int p = N_CHANGE_STATS;
  Vertex h, t;
  
  for (i=0; i < p; i++) 
    mtp->dstats[i] = 0.0;

  for (i=0; i < ntoggles; i++)
    {
      /* is there an edge for this toggle */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);

      hod = nwp->outdegree[h];
      newhod = hod + (edgeflag ? -1 : 1);
      
      for(j=0; j < p; j++) 
	{
	  k =  ((int)mtp->inputparams[j]);
	  bound = (int)mtp->inputparams[j+p];
	  change = MIN(bound,CHOOSE(newhod, k))-MIN(bound,CHOOSE(hod, k));
            	  
	  mtp->dstats[j] += change;
      	}
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_boundedtriangle
*****************/
CHANGESTAT_FN(d_boundedtriangle) { 
  Edge e;
  Vertex h, t, node3;
  Vertex change;
  double boundedchange, htcount;
  Vertex htri, ttri;
  int edgeflag, i;
  int bound = (int)mtp->inputparams[0];

  *(mtp->dstats) = 0.0;

  for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);

      change=0;
      htri=0;
      ttri=0;
      
      for(e = EdgetreeMinimum(nwp->outedges, h);
	 (node3 = nwp->outedges[e].value) != 0;
	  e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
	{
	  htri += CountTriangles(h, node3, 1, 1, nwp);
	}
      
      for(e=EdgetreeMinimum(nwp->inedges, h); 
	 (node3=nwp->inedges[e].value)!=0;
	  e=EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
	  {
	   htri += CountTriangles(h, node3, 1, 1, nwp);
	  }
      
      for(e = EdgetreeMinimum(nwp->outedges, t);
	 (node3 = nwp->outedges[e].value) != 0;
	  e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	{
	  ttri += CountTriangles(t, node3, 1, 1, nwp);
	}
      
      for(e=EdgetreeMinimum(nwp->inedges, t); 
	 (node3=nwp->inedges[e].value)!=0;
	  e=EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	  {
	   ttri += CountTriangles(t, node3, 1, 1, nwp);
	  }
      htri = htri/2;
      ttri = ttri/2;
      htcount = CountTriangles(h, t, 1, 1, nwp);
      boundedchange = (MIN(ttri+(edgeflag ? -1:1)*htcount,bound)-MIN(ttri,bound)+
	               MIN(htri+(edgeflag ? -1:1)*htcount,bound)-MIN(htri,bound));
      *(mtp->dstats) += boundedchange;

      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]);
}

/*****************
 CountTriangles: called by d_boundedtriangle
*****************/
Vertex CountTriangles (Vertex h, Vertex t, int outcount, int incount, 
		       Network *nwp) {
  Edge e;
  Vertex change;
  Vertex k;
  
  change=0;
  if(outcount){
    for(e = EdgetreeMinimum(nwp->outedges, t);
	(k = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
	if (EdgetreeSearch(MIN(k,h), MAX(k,h), nwp->outedges) != 0)
	  ++change;
      }
  }
  
  if(incount){
    for(e = EdgetreeMinimum(nwp->inedges, t); 
	(k = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
	if (EdgetreeSearch(MIN(k,h), MAX(k,h), nwp->outedges) != 0)
	  ++change;
      }
  }
  return(change);
}



/********************  changestats:  C    ***********/
/*****************
 changestat: d_concurrent
*****************/
CHANGESTAT_FN(d_concurrent) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  Vertex b1, b2, actdeg, *od, *id;
  TreeNode *oe=nwp->outedges;

  od=nwp->outdegree;
  id=nwp->indegree;
  *(mtp->dstats) = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    actdeg = od[b1];
    if(!nwp[0].directed_flag){
      actdeg += id[b1];
    }
    *(mtp->dstats) += (actdeg + echange > 1) - (actdeg > 1);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_concurrent_by_attr
*****************/
CHANGESTAT_FN(d_concurrent_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, b1attr;
  Vertex b1, b2, b1deg, *od, *id;
  TreeNode *oe=nwp->outedges;

  od=nwp->outdegree;
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b1deg = od[b1];
    if(!nwp[0].directed_flag){
      b1deg += id[b1];
    }
    b1attr = mtp->inputparams[N_CHANGE_STATS + b1 - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
/*Rprintf("j %d b1deg %d b1attr %d inp %d\n",j,b1deg,b1attr,mtp->inputparams[j]);*/
      if (b1attr == mtp->inputparams[j]) { /* we have attr match */
        mtp->dstats[j] += (b1deg + echange > 1) - (b1deg > 1);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_ctriple
*****************/
CHANGESTAT_FN(d_ctriple) { 
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr;
  
  ninputs = mtp->ninputparams;
  nstats = N_CHANGE_STATS;
  
  if(ninputs>0){
    /* match on attributes */
    if(nstats>1){
      for (j=0; j<nstats; j++){
        mtp->dstats[j] = 0.0;
      }
    }else{
      *(mtp->dstats) = 0.0;
    }
    for (i=0; i<ntoggles; i++) 
    {
      hattr = mtp->attrib[heads[i]-1];
      if(hattr == mtp->attrib[tails[i]-1]){
        edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
        change = 0;
        
        for(e = EdgetreeMinimum(nwp->outedges, t);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
        {
          if(hattr == mtp->attrib[node3-1]){
            if (EdgetreeSearch(node3, h, nwp->outedges) != 0) ++change;
          }
        }
        
        if(nstats>1){
          for (j=0; j<nstats; j++){
            mtp->dstats[j] += (hattr==mtp->inputparams[j]) ? 
            (edgeflag ? -(double)change : change) : 0.0;
          }
        }else{
          *(mtp->dstats) += edgeflag ? -(double)change : change;
        }
        
      }
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    /* no attribute matching */
    *(mtp->dstats) = 0.0;
    for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      
      change = 0;
      
      for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        if (EdgetreeSearch(node3, h, nwp->outedges) != 0) ++change;
      }
      
      *(mtp->dstats) += edgeflag ? -(double)change : change;
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_cycle
*****************/
CHANGESTAT_FN(d_cycle) { 
  int edgeflag,i,j,k,nstats,directed;
  Vertex h, t;
  long int maxlen;
  double *countv,emult;
  
  /*Perform initial setup*/
  directed=(int)(mtp->inputparams[0]);
  maxlen=(long int)(mtp->inputparams[1]);
  nstats=(int)N_CHANGE_STATS;
  countv=(double *)R_alloc(sizeof(double),maxlen-1);
  for(i=0;i<nstats;i++){
    mtp->dstats[i] = 0.0;
  }
      
  for (i=0; i < ntoggles; i++)
    {
      for(j=0;j<maxlen-1;j++)  /*Clear out the count vector*/
        countv[j]=0.0;

      /*Count the cycles associated with this edge*/
      /*Note: the ergm toggle system gets heads and tails reversed!*/
/*      edgewise_cycle_census(g,tails[i],heads[i],countv,maxlen,directed);*/
      edgewise_cycle_census(nwp,heads[i],tails[i],countv,maxlen,directed);

      /*Make the change, as needed*/
/*      edgeflag = (EdgetreeSearch(t=tails[i], h=heads[i], g.outedges) != 0);*/
      if((!directed)&&(heads[i]>tails[i]))
        edgeflag = (EdgetreeSearch(t=tails[i],h=heads[i], nwp->outedges) != 0);
      else
        edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      emult = edgeflag ? -1.0 : 1.0;    /*Set decrement or increment*/
      k=0;
      for(j=0;j<maxlen-1;j++)
        if(mtp->inputparams[2+j]>0.0)
          mtp->dstats[k++]+=emult*countv[j];
      if (i+1 < ntoggles)
        ToggleEdge(h, t, nwp);  /* Toggle this edge if more to come */
    }
  i--;
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]);
}

/*****************
 edgewise_path_recurse:  Called by d_cycle
*****************/
void edgewise_path_recurse(Network *g, Vertex dest, Vertex curnode, Vertex *availnodes, long int availcount, long int curlen, double *countv, long int maxlen, int directed)
{
  Vertex *newavail,i,j;
  long int newavailcount;
  
  /*If we've found a path to the destination, increment the census vector*/ 
  if(directed||(curnode<dest)){
    if(EdgetreeSearch(curnode,dest,g->outedges) != 0)
      countv[curlen]++;
  }else{
    if(EdgetreeSearch(dest,curnode,g->outedges) != 0)
      countv[curlen]++;
  }
  
  /*If possible, keep searching for novel paths*/
  if((availcount>0)&&(curlen<maxlen-2)){
    if(availcount>1){    /*Remove the current node from the available list*/
      if((newavail=(Vertex *)malloc(sizeof(Vertex)*(availcount-1)))==NULL){
        Rprintf("Unable to allocate %d bytes for available node list in edgewise_path_recurse.  Trying to terminate recursion gracefully, but your path count is probably wrong.\n",sizeof(Vertex)*(availcount-1));
        return;
      }
      j=0;
      for(i=0;i<availcount;i++)      /*Create the reduced list, fur passin'*/
        if(availnodes[i]!=curnode)
          newavail[j++]=availnodes[i];
    }else
      newavail=NULL;                 /*Set to NULL if we're out of nodes*/
    newavailcount=availcount-1;      /*Decrement the available count*/
    /*Recurse on all available nodes*/
    for(i=0;i<newavailcount;i++)
      if(directed||(curnode<newavail[i])){
        if(EdgetreeSearch(curnode,newavail[i],g->outedges) != 0)
          edgewise_path_recurse(g,dest,newavail[i],newavail,newavailcount,
            curlen+1,countv,maxlen,directed);
      }else{
        if(EdgetreeSearch(newavail[i],curnode,g->outedges) != 0)
          edgewise_path_recurse(g,dest,newavail[i],newavail,newavailcount,
            curlen+1,countv,maxlen,directed);
      }
    /*Free the available node list*/
    if(newavail!=NULL)
      free((void *)newavail);
  }

  /*Check for interrupts (if recursion is taking way too long...)*/
  R_CheckUserInterrupt();
}

/*****************
 edgewise_cycle_census:  Called by d_cycle
*****************/
void edgewise_cycle_census(Network *g, Vertex t, Vertex h, double *countv, long int maxlen, int directed)
{
  long int n,i,j;
  Vertex *availnodes;

  /*Set things up*/
  n=g->nnodes;

  /*First, check for a 2-cycle (but only if directed)*/
  if(directed&&(EdgetreeSearch(h,t,g->outedges) != 0))
    countv[0]++;
  if(n==2)
    return;                 /*Failsafe for graphs of order 2*/
  
  /*Perform the recursive path count*/
  if((availnodes=(Vertex *)malloc(sizeof(Vertex)*(n-2)))==NULL){
    Rprintf("Unable to allocate %d bytes for available node list in edgewise_cycle_census.  Exiting.\n",sizeof(Vertex)*(n-2));
    return;
  }
  j=0;                             /*Initialize the list of available nodes*/
  for(i=1;i<=n;i++)
    if((i!=h)&&(i!=t))
      availnodes[j++]=i;
  for(i=0;i<n-2;i++)               /*Recurse on each available vertex*/
    if(directed||(h<availnodes[i])){
      if(EdgetreeSearch(h,availnodes[i],g->outedges) != 0)
        edgewise_path_recurse(g,t,availnodes[i],availnodes,n-2,1,countv,maxlen,
          directed);
    }else{
      if(EdgetreeSearch(availnodes[i],h,g->outedges) != 0)
        edgewise_path_recurse(g,t,availnodes[i],availnodes,n-2,1,countv,maxlen,
          directed);
    }
  free((void *)availnodes);  /*Free the available node list*/
}

/********************  changestats:  D    ***********/
/*****************
 changestat: d_degree
*****************/
CHANGESTAT_FN(d_degree) { 
  int i, j, echange;
  Vertex head, tail, headdeg, taildeg, deg, *id, *od;
  TreeNode *oe=nwp->outedges;

  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(head=heads[i], tail=tails[i], oe)==0)? 1:-1;
    headdeg = od[head] + id[head];
    taildeg = od[tail] + id[tail];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      deg = (Vertex)mtp->inputparams[j];
      mtp->dstats[j] += (headdeg + echange == deg) - (headdeg == deg);
      mtp->dstats[j] += (taildeg + echange == deg) - (taildeg == deg);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_degree_by_attr
*****************/
CHANGESTAT_FN(d_degree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr, testattr;
  Vertex head, tail, headdeg, taildeg, d, *id, *od;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    echange=(EdgetreeSearch(head=heads[i], tail=tails[i], oe)==0)? 1:-1;
    headdeg = od[head] + id[head];
    taildeg = od[tail] + id[tail];
    headattr = mtp->inputparams[2*N_CHANGE_STATS + head - 1]; 
    tailattr = mtp->inputparams[2*N_CHANGE_STATS + tail - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)mtp->inputparams[2*j];
      testattr = mtp->inputparams[2*j + 1]; 
      if (headattr == testattr)  /* we have head attr match */
        mtp->dstats[j] += (headdeg + echange == d) - (headdeg == d);
      if (tailattr == testattr)  /* we have tail attr match */
        mtp->dstats[j] += (taildeg + echange == d) - (taildeg == d);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_degree_w_homophily
*****************/
CHANGESTAT_FN(d_degree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr;
  Vertex head, tail, headdeg, taildeg, deg, tmp;
  TreeNode *ie=nwp->inedges, *oe=nwp->outedges;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    head=heads[i];
    tail=tails[i];
    headattr = (int)nodeattr[head];
    tailattr = (int)nodeattr[tail];    
    if (headattr == tailattr) { /* They match; otherwise don't bother */
      echange=(EdgetreeSearch(head, tail, oe)==0)? 1:-1;
      headdeg=taildeg=0;
      for(e = EdgetreeMinimum(oe, head);
      (tmp = oe[e].value) != 0;
      e = EdgetreeSuccessor(oe, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(e = EdgetreeMinimum(ie, head);
      (tmp = ie[e].value) != 0;
      e = EdgetreeSuccessor(ie, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      }
      for(e = EdgetreeMinimum(oe, tail);
      (tmp = oe[e].value) != 0;
      e = EdgetreeSuccessor(oe, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(e = EdgetreeMinimum(ie, tail);
      (tmp = ie[e].value) != 0;
      e = EdgetreeSuccessor(ie, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)mtp->inputparams[j];
        mtp->dstats[j] += (headdeg + echange == deg) - (headdeg == deg);
	mtp->dstats[j] += (taildeg + echange == deg) - (taildeg == deg);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}                                        

/*****************
 changestat: d_density
*****************/
CHANGESTAT_FN(d_density) {
  int edgeflag, i;
  Vertex h, t, ndyads;
  
  ndyads = (N_NODES)*(N_NODES-1);
  if(!nwp->directed_flag){
    ndyads = ndyads / 2;
  }

  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    *(mtp->dstats) += edgeflag ? - 1 : 1;
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = *(mtp->dstats) / ndyads;
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_dsp
*****************/
CHANGESTAT_FN(d_dsp) { 
  Edge e, f;
  int i, j, echange;
  int L2hu, L2ut;
  Vertex deg;
  Vertex h, t, u, v;
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++){      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(nwp->outedges, t);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != h){
        L2hu=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2hu + echange == deg)
          - (L2hu == deg));
        }
      }
    }
    /* step through inedges of t */
    for(e = EdgetreeMinimum(nwp->inedges, t);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != h){
        L2hu=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2hu + echange == deg)
          - (L2hu == deg));
        }
      }
    }
    /* step through outedges of h */
    for(e = EdgetreeMinimum(nwp->outedges, h);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != t){
        L2ut=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2ut + echange == deg)
          - (L2ut == deg));
        }
      }
    }
    /* step through inedges of h */
    for(e = EdgetreeMinimum(nwp->inedges, h);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != t){
        L2ut=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2ut + echange == deg)
          - (L2ut == deg));
        }
      }
    }
    
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_dyadcov
*****************/
CHANGESTAT_FN(d_dyadcov) { 
  double val;
  Vertex h, t;
  int i, edgeflag, refedgeflag;
  long int nrow, noffset, index;
  
  noffset = nwp->bipartite;
  if(noffset > 0){
   nrow = (N_NODES)-(long int)(mtp->inputparams[0]);
  }else{
   nrow = (long int)(mtp->inputparams[0]);
  }
  
/*  Rprintf("nrow %d noffset %d\n",nrow, noffset);
  Rprintf("attrib: ");
  for(i=0;i<1000;i++)
   Rprintf("%1.0f",mtp->attrib[i]);

  Rprintf("\n;"); */

  if(nwp->directed_flag){
  /* directed version */

  for(i=0;i<3;i++)
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial state of the edge and its reflection*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      refedgeflag = (EdgetreeSearch(t, h, nwp->outedges) != 0);
      
      /*Get the dyadic covariate*/
/*    val = mtp->attrib[(t-1-nrow)+(h-1)*ncols]; */
      index = (t-1-noffset)*nrow+(h-1);
      if(index >= 0 && index <= nrow*nrow){
       val = mtp->attrib[(t-1-noffset)*nrow+(h-1)];
/*  Rprintf("h %d t %d nrow %d ncols %d val %f\n",h, t, nrow, ncols, val); */
      
       /*Update the change statistics, as appropriate*/
       if(refedgeflag){      /* Reflected edge is present */
         if(edgeflag){         /* Toggled edge _was_ present */
	  if(t>h){              /* Mut to low->high */
	    mtp->dstats[0] -= val;
	    mtp->dstats[1] += val;
	  }else{                /* Mut to high->low */
	    mtp->dstats[0] -= val;
	    mtp->dstats[2] += val;
	  }
         }else{                /* Toggled edge _was not_ present */
	  if(t>h){              /* Low->high to mut */
	    mtp->dstats[1] -= val;
	    mtp->dstats[0] += val;
	  }else{                /* High->low to mut */
	    mtp->dstats[2] -= val;
	    mtp->dstats[0] += val;
	  }
	}
       }else{                /* Reflected edge is absent */
        if(edgeflag){         /* Toggled edge _was_ present */
	  if(t>h){              /* High->low to null */
	    mtp->dstats[2] -= val;
	  }else{                /* Low->high to null */
	    mtp->dstats[1] -= val;
	  }
        }else{                /* Toggled edge _was not_ present */
	  if(t>h){              /* Null to high->low */
	    mtp->dstats[2] += val;
	  }else{                /* Null to low->high */
	    mtp->dstats[1] += val;
	  }
	}
       }
      }
      
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
/* undirected case (including bipartite) */
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Get the covariate value*/
/*    val = mtp->attrib[(t-1-nrow)+(h-1)*ncols]; */
      index = (t-1-noffset)*nrow+(h-1);
      if(index >= 0 && index <= nrow*((long int)(mtp->inputparams[0]))){
       val = mtp->attrib[(t-1-noffset)*nrow+(h-1)];
      /*Update the change statistic, based on the toggle type*/
/*  Rprintf("h %d t %d nrow %d noffset %d val %f\n",h, t, nrow, noffset, val); */
      /*Update the change statistic, based on the toggle type*/
       *(mtp->dstats) += edgeflag ? -val : val;
      }
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}


/********************  changestats:  E    ***********/
/*****************
 changestat: d_b2concurrent
*****************/
CHANGESTAT_FN(d_b2concurrent) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  Vertex b1, b2, b2deg, *id;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  id=nwp->indegree;
  *(mtp->dstats) = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2];
    *(mtp->dstats) += (b2deg + echange > 1) - (b2deg > 1);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b2concurrent_by_attr
*****************/
CHANGESTAT_FN(d_b2concurrent_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, b2attr;
  Vertex b1, b2, b2deg, *id;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2];
    b2attr = mtp->inputparams[N_CHANGE_STATS + b2 - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b2attr == mtp->inputparams[j]) { /* we have attr match */
        mtp->dstats[j] += (b2deg + echange > 1) - (b2deg > 1);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b2degree
*****************/
CHANGESTAT_FN(d_b2degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, j, echange;
  Vertex b1, b2, b2deg, d, *id;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)(mtp->inputparams[j]);
      mtp->dstats[j] += (b2deg + echange == d) - (b2deg == d);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b2degree_by_attr
*****************/
CHANGESTAT_FN(d_b2degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first 2*nstats values are in pairs:  (degree, attrvalue)
    The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, b2attr;
  Vertex b1, b2, b2deg, d, *id;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2];
    b2attr = mtp->inputparams[2*N_CHANGE_STATS + b2 - 1];
    for(j = 0; j < N_CHANGE_STATS; j++) {
      if (b2attr == mtp->inputparams[2*j+1]) { /* we have attr match */
        d = (Vertex)mtp->inputparams[2*j];
        mtp->dstats[j] += (b2deg + echange == d) - (b2deg == d);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_edgecov
*****************/
CHANGESTAT_FN(d_edgecov) { 
  double val;
  Vertex h, t;
  long int nrow, noffset;
  int i, edgeflag;
  
  noffset = nwp->bipartite;
  if(noffset > 0){
/*   nrow = (N_NODES)-(long int)(mtp->inputparams[0]); */
    nrow = noffset;
  }else{
   nrow = (long int)(mtp->inputparams[0]);
  }
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Get the covariate value*/
/*    val = mtp->attrib[(t-1-nrow)+(h-1)*ncols]; */
      val = mtp->attrib[(t-1-noffset)*nrow+(h-1)];  /*Note: h/t are backwards!*/
/*  Rprintf("h %d t %d nrow %d val %f\n",h, t, nrow, val); */
      /*Update the change statistic, based on the toggle type*/
      *(mtp->dstats) += edgeflag ? -val : val;
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_edges
*****************/
CHANGESTAT_FN(d_edges) {
  int edgeflag, i;
  Vertex h, t;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      *(mtp->dstats) += edgeflag ? - 1 : 1;
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_esp
*****************/
CHANGESTAT_FN(d_esp) { 
  Edge e, f;
  int i, j, echange;
  int L2ht, L2hu, L2ut;
  Vertex deg;
  Vertex h, t, u, v;
  TreeNode *oe=nwp->outedges, *ie=nwp->inedges;

  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++){      
    L2ht=0;
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(oe, t); 
    (u = oe[e].value) != 0; e = EdgetreeSuccessor(oe, e)) {
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), oe) != 0){
        L2ht++;
        L2hu=0;
        L2ut=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(oe, u); (v = oe[f].value) != 0;
        f = EdgetreeSuccessor(oe, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),oe)!= 0) L2ut++;
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),oe)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(ie, u); (v = ie[f].value) != 0;
        f = EdgetreeSuccessor(ie, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),oe)!= 0) L2ut++;
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),oe)!= 0) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2hu + echange == deg)
          - (L2hu == deg));
          mtp->dstats[j] += ((L2ut + echange == deg)
          - (L2ut == deg));
        }
      }
    }
    /* step through inedges of t */
    for (e = EdgetreeMinimum(ie, t); (u = ie[e].value) != 0;
    e = EdgetreeSuccessor(ie, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), oe) != 0){
        L2ht++;
        L2hu=0;
        L2ut=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(oe, u); (v = oe[f].value) != 0;
        f = EdgetreeSuccessor(oe, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),oe)!= 0) L2ut++;
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),oe)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(ie, u); (v = ie[f].value) != 0;
        f = EdgetreeSuccessor(ie, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),oe)!= 0) L2ut++;
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),oe)!= 0) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2hu + echange == deg)
          - (L2hu == deg));
          mtp->dstats[j] += ((L2ut + echange == deg)
          - (L2ut == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)mtp->inputparams[j];
/*      mtp->dstats[j] += echange*((L2ht == deg) - (0 == deg)); */
      mtp->dstats[j] += echange*(L2ht == deg);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_b2factor
*****************/
CHANGESTAT_FN(d_b2factor) { 
  double s, factorval;
  Vertex h, t;
  int i, j;
  long int nb1, nb2;
  
  nb1 = nwp[0].bipartite;
  nb2 = nwp[0].nnodes - nb1;

  for (i=0; i < N_CHANGE_STATS; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[t-nb1-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}
/********************  changestats:  F    ***********/

/********************  changestats:  G    ***********/
/*****************
 changestat: d_gwb1degree
*****************/
CHANGESTAT_FN(d_gwb1degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  double decay, oneexpd;
  Vertex b1, b1deg, *od;
  TreeNode *oe;  
  
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  oe=nwp->outedges;
  od=nwp->outdegree;
  mtp->dstats[0] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], tails[i], oe)==0) ? 1 : -1;
    b1deg = od[b1]+(echange-1)/2;
    mtp->dstats[0] += echange*pow(oneexpd,(double)b1deg);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwb1degree_by_attr
*****************/
CHANGESTAT_FN(d_gwb1degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first value is theta (as in Hunter et al, JASA 200?), controlling decay
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, echange, b1attr;
  double decay, oneexpd;
  Vertex b1, b1deg, *od;
  TreeNode *oe;
  
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  oe=nwp->outedges;
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(b1=heads[i], tails[i], oe)==0) ? 1 : -1;
    b1deg = od[b1]+(echange-1)/2;
    b1attr = mtp->inputparams[b1]; 
/*  Rprintf("b1 %d tails %d b1deg %d b1attr %d echange %d\n",b1, tails[i], b1deg, b1attr, echange); */
    mtp->dstats[b1attr-1] += echange * pow(oneexpd,(double)b1deg);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwdegree
*****************/
CHANGESTAT_FN(d_gwdegree) { 
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
    hd = od[h] + id[h] + (echange - 1)/2;
    td = od[t] + id[t] + (echange - 1)/2;
    change += echange*(pow(oneexpd,(double)hd)+pow(oneexpd,(double)td));
      
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwdegree_by_attr
*****************/
CHANGESTAT_FN(d_gwdegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, hattr, tattr, echange=0;
  double decay, oneexpd;
  Vertex h, t, hd, td=0, *id, *od;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    hd = od[h] + id[h] + (echange - 1)/2;
    hattr = mtp->inputparams[h]; 
    mtp->dstats[hattr-1] += echange*(pow(oneexpd,(double)hd));
    
    td = od[t] + id[t] + (echange - 1)/2;
    tattr = mtp->inputparams[t]; 
    mtp->dstats[tattr-1] += echange*(pow(oneexpd,(double)td));
      
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwdsp
****************/
CHANGESTAT_FN(d_gwdsp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2hu, L2ut;
  Vertex h, t, u, v;
  double alpha, oneexpa, cumchange;
  
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++){      
    cumchange=0.0;
    ochange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(nwp->outedges, t);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != h){
        L2hu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        cumchange += pow(oneexpa,(double)L2hu);
      }
    }
    /* step through inedges of t */
    for(e = EdgetreeMinimum(nwp->inedges, t);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != h){
        L2hu=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u);
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
        }
        cumchange += pow(oneexpa,(double)L2hu);
      }
    }
    
    /* step through outedges of h  */
    for(e = EdgetreeMinimum(nwp->outedges, h);
    (u = nwp->outedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->outedges, e)){
      if (u != t){
        L2ut=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        cumchange += pow(oneexpa,(double)L2ut);
      }
    }
    /* step through inedges of h */
    for(e = EdgetreeMinimum(nwp->inedges, h);
    (u = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)){
      if (u != t){
        L2ut=ochange;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(nwp->outedges, u);
        (v = nwp->outedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->outedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        /* step through inedges of u */
        for(f = EdgetreeMinimum(nwp->inedges, u); 
        (v = nwp->inedges[f].value) != 0;
        f = EdgetreeSuccessor(nwp->inedges, f)){
          if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
        }
        cumchange += pow(oneexpa,(double)L2ut);
      }
    }
    
    cumchange  = echange*cumchange;
    (*(mtp->dstats)) += cumchange;
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwb2degree
*****************/
CHANGESTAT_FN(d_gwb2degree) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  */
  int i, echange;
  double decay, oneexpd;
  Vertex b2, b2deg, *id;
  TreeNode *oe;  
  
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  oe=nwp->outedges;
  id=nwp->indegree;
  mtp->dstats[0] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2]+(echange-1)/2;
    mtp->dstats[0] += echange*pow(oneexpd,(double)b2deg);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwb2degree_by_attr
*****************/
CHANGESTAT_FN(d_gwb2degree_by_attr) { 
  /* It is assumed that in this bipartite network, the only edges are
  of the form (b1, b2), where b1 is always strictly less
  than b2.  In other words, the degree of a b1 is equivalent
  to its outdegree and the degree of a b2 is equivalent to its
  indegree.
  The inputparams are assumed to be set up as follows:
    The first value is theta (as in Hunter et al, JASA 200?), controlling decay
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, echange, b2attr;
  double decay, oneexpd;
  Vertex b2, b2deg, *id;
  TreeNode *oe;
  
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  oe=nwp->outedges;
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(heads[i], b2=tails[i], oe)==0) ? 1 : -1;
    b2deg = id[b2]+(echange-1)/2;
    b2attr = mtp->inputparams[b2]; 
/*  Rprintf("h %d b2 %d b2deg %d b2attr %d echange %d\n",heads[i], b2, b2deg, b2attr, echange); */
    mtp->dstats[b2attr-1] += echange * pow(oneexpd,(double)b2deg);
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwesp
*****************/
CHANGESTAT_FN(d_gwesp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2ht, L2hu, L2ut;
  Vertex h, t, u, v;
  double alpha, oneexpa, cumchange;
  
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++){      
    cumchange=0.0;
    L2ht=0;
    ochange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of t  */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=ochange;
	L2ut=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
	  if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
	  if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
	}
	cumchange += pow(oneexpa,(double)L2hu) +
	  pow(oneexpa,(double)L2ut) ;
      }
    }
    /* step through inedges of t */
    
    for(e = EdgetreeMinimum(nwp->inedges, t);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=ochange;
	L2ut=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
	  if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
	}
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(MIN(v,t),MAX(v,t),nwp->outedges)!= 0) L2ut++;
	  if(EdgetreeSearch(MIN(v,h),MAX(v,h),nwp->outedges)!= 0) L2hu++;
	}
	cumchange += pow(oneexpa,(double)L2hu) +
	  pow(oneexpa,(double)L2ut) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2ht)) ;
    }else{
      cumchange += (double)L2ht;
    }
    cumchange  = echange*cumchange;
    (*(mtp->dstats)) += cumchange;
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwidegree
*****************/
CHANGESTAT_FN(d_gwidegree) { 
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex t, td=0, *id;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    td = id[t] + (echange - 1)/2;
    change += echange * pow(oneexpd,(double)td);
    
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwidegree_by_attr
*****************/
CHANGESTAT_FN(d_gwidegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, tattr, echange=0;
  double decay, oneexpd;
  Vertex t, td=0, *id;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    td = id[t] + (echange - 1)/2;
    tattr = mtp->inputparams[t]; 
    mtp->dstats[tattr-1] += echange*(pow(oneexpd,(double)td));
      
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwodegree
*****************/
CHANGESTAT_FN(d_gwodegree) { 
  int i, echange=0;
  double decay, oneexpd;
  Vertex h, hd=0, *od;
  TreeNode *oe=nwp->outedges;
  
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], tails[i], oe) == 0) ? 1 : -1;
    hd = od[h] + (echange - 1)/2;
    mtp->dstats[0] += echange * pow(oneexpd,(double)hd);
    
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwodegree_by_attr
*****************/
CHANGESTAT_FN(d_gwodegree_by_attr) { 
  /*The inputparams are assumed to be set up as follows:
    The first value is the decay parameter (as in Hunter et al, JASA 200?)
    The next sequence of values is the nodal attributes, coded as integers
         from 1 through N_CHANGE_STATS
  */
  int i, hattr, echange=0;
  double decay, oneexpd;
  Vertex h, hd, *od;
  TreeNode *oe=nwp->outedges;
  
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], tails[i], oe) == 0) ? 1 : -1;
    hd = od[h] + (echange - 1)/2;
    hattr = mtp->inputparams[h]; 
    mtp->dstats[hattr-1] += echange*(pow(oneexpd,(double)hd));
      
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwtdsp
****************/
CHANGESTAT_FN(d_gwtdsp) {
  Edge e, f;
  int i, echange, ochange, L2hu, L2ut;
  Vertex h, t, u, v;
  double alpha, oneexpa, cumchange;
  
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++){
    h=heads[i]; t=tails[i];
    cumchange=0.0;
    ochange = -IS_OUTEDGE(h,t);
    echange = 2*ochange + 1;
    /* step through outedges of t */
    for(e = MIN_OUTEDGE(t); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != h){
        L2hu=ochange; /* L2hu will be # shrd prtnrs of (h,u) not incl. t */
        /* step through inedges of u, incl. (t,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(h,v)) L2hu++;
        }
        cumchange += pow(oneexpa,(double)L2hu); /* sign corrected below */
      }
    }
    /* step through inedges of h */
    for(e = MIN_INEDGE(h); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != t){
        L2ut=ochange; /* L2ut will be # shrd prtnrs of (u,t) not incl. h */
        /* step through outedges of u , incl. (u,h) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,t)) L2ut++;
        }
        cumchange += pow(oneexpa,(double)L2ut); /* sign corrected below */
      }
    }
    CHANGE_STAT[0] += echange * cumchange;
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_gwtesp
*****************/
CHANGESTAT_FN(d_gwtesp) { 
  Edge e, f;
  int i, echange, ochange;
  int L2ht, L2hu, L2ut;
  Vertex h, t, u, v;
  double alpha, oneexpa, cumchange;
  
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++){      
    cumchange=0.0;
    L2ht=0;
    ochange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 0 : -1;
    echange = 2*ochange + 1;
    /* step through outedges of t  */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(h, u, nwp->outedges) != 0){
	L2hu=ochange;
	/* step through inedges of u */
	for(f = EdgetreeMinimum(nwp->inedges, u); 
	    (v = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)){
	  if(EdgetreeSearch(h,v,nwp->outedges)!= 0) L2hu++;
	}
	cumchange += pow(oneexpa,(double)L2hu);
      }
    }
    /* step through inedges of t */
    
    for(e = EdgetreeMinimum(nwp->inedges, t);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(h, u, nwp->outedges) != 0){
	L2ht++;
      }
      if (EdgetreeSearch(u, h, nwp->outedges) != 0){
	L2ut=ochange;
	/* step through outedges of u */
	for(f = EdgetreeMinimum(nwp->outedges, u);
	    (v = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)){
	  if(EdgetreeSearch(v,t,nwp->outedges)!= 0) L2ut++;
	}
	cumchange += pow(oneexpa,(double)L2ut) ;
      }
    }
    
    if(alpha < 100.0){
      cumchange += exp(alpha)*(1.0-pow(oneexpa,(double)L2ht)) ;
    }else{
      cumchange += (double)L2ht;
    }
    cumchange  = echange*cumchange;
    (*(mtp->dstats)) += cumchange;
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  H    ***********/
/*****************
 changestat: d_hamming
*****************/
CHANGESTAT_FN(d_hamming) { 
  Vertex h, t;
  int i, nhedge, discord;
  
  nhedge = nwp[1].nedges;
/*Rprintf("nhedge %d\n",nhedge); */
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    /*Get the initial state of the edge and its alter in x0*/
/*  edgeflag =(EdgetreeSearch(h=heads[i], t=tails[i], nwp[0].outedges) != 0); */
    discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
/*    if(!nwp[0].directed_flag && h < t){
      hh = t;
      ht = h;
    }else{
      hh = h;
      ht = t;
    }
     if we will dissolve an edge discord=-1
     discord = edgeflag ? -1 : 1;
    

  so moving away one step
    discord = (edgeflag0!=edgeflag) ? -1 : 1;

Rprintf("h %d t %d discord %d\n",h, t, discord);
  if(nhedge>0)
  Rprintf("h %d t %d discord %d nhedge %d\n",h, t, discord, nhedge); */

    /*Update the change statistics, as appropriate*/
/*    *(mtp->dstats) += ((edgeflag0!=edgeflag) ? -1.0 : 1.0); */

    *(mtp->dstats) += (discord ? -1.0 : 1.0);

    if (i+1 < ntoggles){
      ToggleEdge(heads[i], tails[i], &nwp[0]);  /* Toggle this edge if more to come */
      ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
    }
  }
  i--; 
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], &nwp[0]); 
    ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
  }
}

/*****************
 changestat: d_hamming_weighted
*****************/
CHANGESTAT_FN(d_hamming_weighted) { 
  Vertex h, t;
  double val;
  long int nnodes, nb1, nb2, n0edge;
  int i, discord;

  n0edge =  mtp->inputparams[0];
  nnodes = nwp[0].nnodes;
  nb1 = nwp[0].bipartite;
  nb2 = nwp[0].nnodes - nb1;
/*  Rprintf("nb1 %d i0 %f i1 %f i2 %f i3 %f\n", nb1,
                                 mtp->inputparams[0],
                                 mtp->inputparams[1],
                                 mtp->inputparams[2],
                                 mtp->inputparams[3]
		  );
  for (i=0; i<1000; i++) {
  Rprintf("i %d inp %f\n", i, mtp->inputparams[i]);
  } */

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial discord state*/
/*    edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[0].outedges) != 0); */
      discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
      /*Get the covariate value*/
      val = mtp->inputparams[1+(t-nb1-1)*nb1+(h-1)+2*n0edge];
      /*Update the change statistic, based on the toggle type*/
      *(mtp->dstats) += discord ? -val : val;
 /* Rprintf("nnodes %d n0edge %d h %d t %d discord %d val %f\n",nnodes, n0edge, h, t-nb1, discord, val); */
      if (i+1 < ntoggles){
        ToggleEdge(heads[i], tails[i], &nwp[0]);  /* Toggle this edge if more to come */
        ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
      }
  }
  i--; 
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], &nwp[0]); 
    ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
  }
}

/*****************
 changestat: d_hammingmix_constant
*****************/
CHANGESTAT_FN(d_hammingmix_constant) { 
  Vertex h, t;
  int i, nhedge, discord;
  int matchvalh, matchvalt;
  
  nhedge = mtp->inputparams[0];
/*  Rprintf("nhedge %d\n", nhedge); */
  if(ntoggles==2){
   matchvalh = mtp->inputparams[heads[0]+2*nhedge];
   matchvalt = mtp->inputparams[tails[0]+2*nhedge];
   if(matchvalh != mtp->inputparams[heads[1]+2*nhedge] ||
      matchvalt != mtp->inputparams[tails[1]+2*nhedge]){
      *(mtp->dstats) = 10000.0;
      return;
   }
  }
  *(mtp->dstats) = 0.0;
/*  Rprintf("Warning: hammingconstantmix can only be used with ConstantEdges terms.\n");
  Rprintf("nhedge %d i0 %f i1 %f i2 %f i3 %f\n", nhedge, mtp->inputparams[0],
                                 mtp->inputparams[1],
                                 mtp->inputparams[2],
                                 mtp->inputparams[3]
		  ); */
     
  for (i=0; i<ntoggles; i++)
    {
      discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
      *(mtp->dstats) += (discord ? -1.0 : 1.0);

    if (i+1 < ntoggles){
      ToggleEdge(heads[i], tails[i], &nwp[0]);  /* Toggle this edge if more to come */
      ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
    }
  }
  i--; 
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], &nwp[0]); 
    ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
  }
}

/*****************
 changestat: d_hammingmix
*****************/
CHANGESTAT_FN(d_hammingmix) { 
  Vertex h, t;
  int i, j, nhedge, edgeflag, discord;
  int matchvalh, matchvalt;
  int nstats;
  
  nhedge =  mtp->inputparams[0];
  nstats = N_CHANGE_STATS;
/*  Rprintf("nstats %d nhedge %d i0 %f i1 %f i2 %f i3 %f\n",nstats, nhedge, mtp->inputparams[0],
                                 mtp->inputparams[1],
                                 mtp->inputparams[2],
                                 mtp->inputparams[3]
		  ); */
  for (i=0; i < N_CHANGE_STATS; i++)
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++)
    {
      h=heads[i];
      t=tails[i];
      matchvalh = mtp->inputparams[h+2*nstats+2*nhedge];
      matchvalt = mtp->inputparams[t+2*nstats+2*nhedge];
      edgeflag=(EdgetreeSearch(h, t, nwp[0].outedges) != 0); /*Get edge state*/
      discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
      for (j=0; j<nstats; j++) 
	  {
/*   Rprintf("h %d t %d matchvalh %d matchvalt %d edgeflag %d discord %d j %d p0 %f p1 %f\n",h,t,matchvalh,matchvalt,edgeflag,discord,j,mtp->inputparams[2*nhedge+  j], mtp->inputparams[2*nhedge+ nstats+j]); */
           if(matchvalh==mtp->inputparams[2*nhedge+1+       j] &&
	      matchvalt==mtp->inputparams[2*nhedge+1+nstats+j]
	     ){
		mtp->dstats[j] += (discord ? -1.0 : 1.0);
	      }
	  }

    if (i+1 < ntoggles){
      ToggleEdge(heads[i], tails[i], &nwp[0]);  /* Toggle this edge if more to come */
      ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
    }
  }
  i--; 
  while (--i>=0){  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], &nwp[0]); 
    ToggleEdge(heads[i], tails[i], &nwp[1]);  /* Toggle the discord for this edge */
  }
}

/********************  changestats:  I    ***********/
/*****************
 changestat: d_idegree
*****************/
CHANGESTAT_FN(d_idegree) { 
  int echange, i, j;
  Edge e;
  Vertex h, t, node3, deg, td=0;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	hattr = mtp->attrib[h-1];
	if(hattr == mtp->attrib[t-1]){
	  td = 0;
	  for(e = EdgetreeMinimum(nwp->inedges, t);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	    {
	      if(hattr == mtp->attrib[node3-1]){++td;}
	    }
	  
	  for(j=0; j < N_CHANGE_STATS; j++) 
	    {
	      deg = (Vertex)mtp->inputparams[j];
	      mtp->dstats[j] += (td + echange == deg) - (td == deg);
	    }
	}
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	td = nwp->indegree[t];
	
	for(j=0; j < N_CHANGE_STATS; j++) 
	  {
	    deg = (Vertex)mtp->inputparams[j];
	    mtp->dstats[j] += (td + echange == deg) - (td == deg);
	  }
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_idegree_by_attr
*****************/
CHANGESTAT_FN(d_idegree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, tailattr, testattr;
  Vertex tail, taildeg, d, *id;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    echange=(EdgetreeSearch(heads[i], tail=tails[i], oe)==0)? 1:-1;
    taildeg = id[tail];
    tailattr = mtp->inputparams[2*N_CHANGE_STATS + tail - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)mtp->inputparams[2*j];
      testattr = mtp->inputparams[2*j + 1]; 
      if (tailattr == testattr)  /* we have tail attr match */
        mtp->dstats[j] += (taildeg + echange == d) - (taildeg == d);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_idegree_w_homophily
*****************/
CHANGESTAT_FN(d_idegree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr;
  Vertex head, tail, taildeg, deg, tmp;
  TreeNode *ie=nwp->inedges, *oe=nwp->outedges;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    head=heads[i];
    tail=tails[i];
    headattr = (int)nodeattr[head];
    tailattr = (int)nodeattr[tail];    
    if (headattr == tailattr) { /* They match; otherwise don't bother */
      echange=(EdgetreeSearch(head, tail, oe)==0)? 1:-1;
      taildeg=0;
/*      for(e = EdgetreeMinimum(oe, tail);
      (tmp = oe[e].value) != 0;
      e = EdgetreeSuccessor(oe, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      } */
      for(e = EdgetreeMinimum(ie, tail);
      (tmp = ie[e].value) != 0;
      e = EdgetreeSuccessor(ie, e)) {
        taildeg += (nodeattr[tmp]==tailattr);
      }
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)mtp->inputparams[j];
        mtp->dstats[j] += (taildeg + echange == deg) - (taildeg == deg);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_intransitive
*****************/
CHANGESTAT_FN(d_intransitive) { 
  Edge e;
  Vertex h, t, node2;
  double change;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
/*           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag); */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node2 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (node2 != h){
       if (EdgetreeSearch(h, node2, nwp->outedges) == 0){
        change = change + 1.0;
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, t);
	    (node2 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      if (node2 != h){
       if (EdgetreeSearch(h, node2, nwp->outedges) != 0){
        change = change - 1.0;
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, h);
	    (node2 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of hail */
    {
      if (node2 != t){
       if (EdgetreeSearch(node2, t, nwp->outedges) == 0){
        change = change + 1.0;
       }
      }
    }
    
    *(mtp->dstats) += edgeflag ? -change : change;
/*  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change); */

    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_isolates
*****************/
CHANGESTAT_FN(d_isolates) { 
  int i, echange;
  Vertex h, t, hd, td=0, *id, *od;
  TreeNode *oe=nwp->outedges;

  id=nwp->indegree;
  od=nwp->outdegree;
  
  *(mtp->dstats) = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
      hd = od[h] + id[h];
      td = od[t] + id[t];
      *(mtp->dstats) += (hd + echange == 0) - (hd == 0);
      *(mtp->dstats) += (td + echange == 0) - (td == 0);
      
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_istar
*****************/
CHANGESTAT_FN(d_istar) { 
  double change, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hattr = mtp->attrib[h-1];
      if(hattr == mtp->attrib[t-1]){
        td = - edgeflag;
        for(e = EdgetreeMinimum(nwp->inedges, t);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) {/* step through inedges of tail */
          if(hattr == mtp->attrib[node3-1]){++td;}
        }	  
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)mtp->inputparams[j]) - 1;
          change = CHOOSE(td, kmo); 
          mtp->dstats[j] += (edgeflag ? - change : change); 
        }
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      td = nwp->indegree[t] - edgeflag;	
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)mtp->inputparams[j]) - 1;
        change = CHOOSE(td, kmo); 
        mtp->dstats[j] += (edgeflag ? - change : change); 
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  K    ***********/
/*****************
 changestat: d_kstar
*****************/
CHANGESTAT_FN(d_kstar) { 
  double change, hd, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
    
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hattr = mtp->attrib[h-1];
      if(hattr == mtp->attrib[t-1]){
        hd = - edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, h);
        (node3 = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) {/* step through outedges of head */
          if(hattr == mtp->attrib[node3-1]){++hd;}
        }
        for(e = EdgetreeMinimum(nwp->inedges, h);
        (node3 = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) { /* step through inedges of head */
          if(hattr == mtp->attrib[node3-1]){++hd;}
        }
        td = - edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, t);
        (node3 = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) {/* step through outedges of tail */
          if(hattr == mtp->attrib[node3-1]){++td;}
        }
        for(e = EdgetreeMinimum(nwp->inedges, t);
        (node3 = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) {/* step through inedges of tail */
          if(hattr == mtp->attrib[node3-1]){++td;}
        }
        
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)mtp->inputparams[j]) - 1;
/*          if (kmo==0) {
            change=1;
          } else { */
            change = CHOOSE(hd, kmo) + CHOOSE(td, kmo); 
/*          } uncomment these few lines to define 1-stars as equivalent to 
              edges (currently, each edge is counted as two 1-stars) */
          mtp->dstats[j] += (edgeflag ? - change : change); 
        }
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    for (i=0; i < ntoggles; i++)
    {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hd = nwp->outdegree[h] + nwp->indegree[h] - edgeflag; 
      td = nwp->outdegree[t] + nwp->indegree[t] - edgeflag;
      for(j=0; j < N_CHANGE_STATS; j++) 
      {
        kmo = ((int)mtp->inputparams[j]) - 1;
/*        if (kmo==0) {
          change=1;
        } else { */
          change = CHOOSE(hd, kmo) + CHOOSE(td, kmo); 
/*      } uncomment these few lines to define 1-stars as equivalent to 
          edges (currently, each edge is counted as two 1-stars) */
        mtp->dstats[j] += (edgeflag ? - change : change); 
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  L    ***********/
/*****************
 changestat: d_localtriangle
*****************/
CHANGESTAT_FN(d_localtriangle) { 
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i;
  int ninputs, nstats;
  long int nmat;
  
  ninputs = mtp->ninputparams;
  nstats = N_CHANGE_STATS;
  nmat = (long int)(mtp->inputparams[0]);
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      change = 0;
      
      if(mtp->inputparams[1+(tails[i]-1)+(heads[i]-1)*nmat] != 0){
	for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	  {
	    if(mtp->inputparams[1+(node3-1)+(heads[i]-1)*nmat] != 0 && 
	       mtp->inputparams[1+(node3-1)+(tails[i]-1)*nmat] != 0 ){
	      if (nwp->directed_flag){
		if (EdgetreeSearch(node3, h, nwp->outedges) != 0) ++change;
		if (EdgetreeSearch(node3, h, nwp->inedges) != 0) ++change;
	      }else{
		if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		  ++change;
		}
	      }
	    }
	  }
	
	for(e = EdgetreeMinimum(nwp->inedges, t); 
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	  {
	    if(mtp->inputparams[1+(node3-1)+(heads[i]-1)*nmat] != 0 && 
	       mtp->inputparams[1+(node3-1)+(tails[i]-1)*nmat] != 0 ){
	      if (nwp->directed_flag)
		{
		  if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
		    ++change;
		  if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
		    ++change;
		}
	      else
		{
		  if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		    ++change;
		  }
		}
	    }
	  }
	
	*(mtp->dstats) += edgeflag ? -(double)change : change;
      }
      
      if (i+1 < ntoggles) 
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  M    ***********/
/*****************
 changestat: d_m2star
*****************/
CHANGESTAT_FN(d_m2star) {
  Vertex h, t;
  int hid, tod, change;
  int i, edgeflag, backedgeflag;
    
  *(mtp->dstats) = 0.0;

  for (i=0; i < ntoggles; i++)
    {
      /*  edgeflag is 1 if the edge from heads[i] to tails[i]  */
      /*   exists and will disappear */
      /*  edgeflag is 0 if the edge does not exist */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      backedgeflag = (EdgetreeSearch(t, h, nwp->outedges) != 0);

      hid = nwp->indegree[h]; 
      tod = nwp->outdegree[t];
      change = hid + tod - 2*backedgeflag; 
      *(mtp->dstats) += (edgeflag ? -change : change); 

      if (i+1 < ntoggles)
        ToggleEdge(h, t, nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_meandeg
*****************/
CHANGESTAT_FN(d_meandeg) {
  int edgeflag, i;
  Vertex h, t;
  
  *(mtp->dstats) = 0.0;

  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      *(mtp->dstats) += (edgeflag ? - 2 : 2);
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  
  *(mtp->dstats)/=N_NODES;

  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_mix
 This appears to be the version of nodemix used for 
 bipartite networks (only)
*****************/
CHANGESTAT_FN(d_mix) {
  Vertex h, t, tmpi;
  int matchvalh, matchvalt;
  int i, j, edgeflag=0, nstats;

  nstats = N_CHANGE_STATS;
  for (i=0; i < N_CHANGE_STATS; i++)
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    h=heads[i];
    t=tails[i];
    edgeflag=(EdgetreeSearch(h, t, nwp[0].outedges) != 0); /*Get edge state*/
    if (nwp->bipartite > 0 && h > t) { 
      tmpi = h; h = t; t = tmpi; /* swap h, t */
    }
    matchvalh = mtp->inputparams[h-1+2*nstats];
    matchvalt = mtp->inputparams[t-1+2*nstats];
    for (j=0; j<nstats; j++) {
      if(matchvalh==mtp->inputparams[       j] &&
	      matchvalt==mtp->inputparams[nstats+j]
      ){mtp->dstats[j] += edgeflag ? -1.0 : 1.0;}
	  }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], &nwp[0]);  /* Toggle this edge if more to come */
  }

  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], &nwp[0]);
}

/*****************
 changestat: d_mutual

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
CHANGESTAT_FN(d_mutual) { 
  Vertex h, t;
  int i, edgeflag, refedgeflag;
  
  *(mtp->dstats) = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      refedgeflag = (EdgetreeSearch(t, h, nwp->outedges) != 0);
      
      /* from (1,1) to anything = -1 */
      if (edgeflag == 1 && refedgeflag == 1)
	*(mtp->dstats) -= 1;
      
      /* from anything to (1,1) = +1 */
      if (edgeflag == 0 && refedgeflag == 1)
	*(mtp->dstats) += 1;
      
      if (i+1 < ntoggles) 
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  N    ***********/
/*****************
 changestat: d_nearsimmelian
*****************/
CHANGESTAT_FN(d_nearsimmelian) { 
  Vertex h, t, node3;
  double change;
  int edgeflag, i, edgeflagth, sc;

 *(mtp->dstats) = 0.0;

 for (i=0; i<ntoggles; i++) 
 {
  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
  edgeflagth = (EdgetreeSearch(t, h, nwp->outedges) == 0);
   
  for(node3=1;node3<=N_NODES;node3++){
    if((node3!=h)&&(node3!=t)){
     sc = edgeflagth + (EdgetreeSearch(node3, h, nwp->outedges) == 0);
     if(sc < 2){
      sc += (EdgetreeSearch(h, node3, nwp->outedges) == 0);
      if(sc < 2){
       sc += (EdgetreeSearch(node3, t, nwp->outedges) == 0);
       if(sc < 2){
        sc += (EdgetreeSearch(t, node3, nwp->outedges) == 0);
        if(sc < 2){
         change=0.0;
         if (sc == 0 && edgeflag == 0 ){--change;}
         if (sc == 0 && edgeflag == 1 ){++change;}
         if (sc == 1 && edgeflag == 0 ){++change;}
         if (sc == 1 && edgeflag == 1 ){--change;}
         *(mtp->dstats) += change;
	}
       }
      }
     }
    }
   }
   
   if (i+1 < ntoggles) 
     TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodecov
*****************/
CHANGESTAT_FN(d_nodecov) { 
  double sum;
  Vertex h, t;
  int i, edgeflag;

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      sum = mtp->attrib[h-1] + mtp->attrib[t-1];
      *(mtp->dstats) += edgeflag ? -sum : sum;
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodefactor
*****************/
CHANGESTAT_FN(d_nodefactor) { 
  double s, factorval;
  Vertex h, t;
  int i, j;
  
  for (i=0; i < N_CHANGE_STATS; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[h-1] != factorval) ? 0.0 : s)
	    + ((mtp->attrib[t-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodeicov
*****************/
CHANGESTAT_FN(d_nodeicov) { 
  double sum;
  Vertex h, t;
  int i, edgeflag;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      sum = mtp->attrib[t-1];
      *(mtp->dstats) += edgeflag ? -sum : sum;
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodeifactor
*****************/
CHANGESTAT_FN(d_nodeifactor) { 
  double s, factorval;
  Vertex h, t;
  int i, j;
  
  for (i=0; i < N_CHANGE_STATS; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[t-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodematch
*****************/
CHANGESTAT_FN(d_nodematch) { 
  double matchval/*, checksum=0.0*/;
  Vertex h, t, ninputs;
  int i, j, edgeflag=0, matchflag;
  
  ninputs = mtp->ninputparams - N_NODES;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    h=heads[i];
    t=tails[i];
    matchval = mtp->inputparams[h+ninputs-1];
    matchflag = (matchval == mtp->inputparams[t+ninputs-1]);
    if (matchflag) 
      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
    if (ninputs==0) {/* diff=F in network statistic specification */
      *(mtp->dstats) += matchflag ? (edgeflag ? -1.0 : 1.0) : 0.0;
    }
    else { /* diff=T (and more than one category?)  */
      for (/*checksum=0.0,*/ j=0; j<ninputs; j++) {
	      mtp->dstats[j] += (matchflag && matchval==mtp->inputparams[j]) ? 
        (edgeflag ? -1.0 : 1.0) : 0.0;
	      /*checksum += mtp->dstats[j];*/
	    }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodemix
 Update mixing matrix, non-bipartite networks only 
 (but see also d_mix)
*****************/
CHANGESTAT_FN(d_nodemix) {
  Vertex h, t, ninputs, ninputs2;
  int i, j, edgeflag=0, matchflag;
  double rtype, ctype, tmp;

  ninputs = mtp->ninputparams - N_NODES;
  ninputs2 = ninputs/2;

  for (i=0; i < N_CHANGE_STATS; i++)
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++)
    {
      h=heads[i];
      t=tails[i];
      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0); /*Get edge state*/
      matchflag=0;
      /*Find the node covariate values (types) for the head and tail*/
      rtype=mtp->inputparams[h+ninputs-1];
      ctype=mtp->inputparams[t+ninputs-1];
      if (!nwp->directed_flag && rtype > ctype)  {
        tmp = rtype; rtype = ctype; ctype = tmp; /* swap rtype, ctype */
      }
      /*Find the right statistic to update */
      for(j=0;(j<ninputs2)&&(!matchflag);j++){
        if((mtp->inputparams[j          ]==rtype)&&
           (mtp->inputparams[j+ninputs2]==ctype)){
            mtp->dstats[j] += (edgeflag ? -1.0 : 1.0);
            matchflag++;
        }
      } 
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]);
}

/*****************
 changestat: d_nodeocov
*****************/
CHANGESTAT_FN(d_nodeocov) { 
  double sum;
  Vertex h, t;
  int i, edgeflag;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      sum = mtp->attrib[h-1];
      *(mtp->dstats) += edgeflag ? -sum : sum;
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_nodeofactor
*****************/
CHANGESTAT_FN(d_nodeofactor) { 
  double s, factorval;
  Vertex h, t;
  int i, j;
  
  for (i=0; i < N_CHANGE_STATS; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(N_CHANGE_STATS); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[h-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  O    ***********/
/*****************
 changestat: d_odegree
*****************/
CHANGESTAT_FN(d_odegree) { 
  int echange, i, j;
  Edge e;
  Vertex h, t, node3, deg, td=0;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	hattr = mtp->attrib[t-1];
	if(hattr == mtp->attrib[h-1]){
	  td = 0;
	  for(e = EdgetreeMinimum(nwp->outedges, h);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
	    {
	      if(hattr == mtp->attrib[node3-1]){++td;}
	    }
	  
	  for(j=0; j < N_CHANGE_STATS; j++) 
	    {
	      deg = (Vertex)mtp->inputparams[j];
	      mtp->dstats[j] += (td + echange == deg) - (td == deg);
	    }
	}
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	td = nwp->outdegree[h];
	
	for(j=0; j < N_CHANGE_STATS; j++) 
	  {
	    deg = (Vertex)mtp->inputparams[j];
	    mtp->dstats[j] += (td + echange == deg) - (td == deg);
	  }
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]);    
}

/*****************
 changestat: d_odegree_by_attr
*****************/
CHANGESTAT_FN(d_odegree_by_attr) { 
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, testattr;
  Vertex head, headdeg, d, *od;
  TreeNode *oe=nwp->outedges;
  
  od=nwp->outdegree;
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    echange=(EdgetreeSearch(head=heads[i], tails[i], oe)==0)? 1:-1;
    headdeg = od[head];
    headattr = mtp->inputparams[2*N_CHANGE_STATS + head - 1]; 
    for(j = 0; j < N_CHANGE_STATS; j++) {
      d = (Vertex)mtp->inputparams[2*j];
      testattr = mtp->inputparams[2*j + 1]; 
      if (headattr == testattr) { /* we have head attr match */
        mtp->dstats[j] += (headdeg + echange == d) - (headdeg == d);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_odegree_w_homophily
*****************/
CHANGESTAT_FN(d_odegree_w_homophily) { 
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr;
  Vertex head, tail, headdeg, deg, tmp;
  TreeNode *oe=nwp->outedges;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + N_CHANGE_STATS - 1;  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    head=heads[i];
    tail=tails[i];
    headattr = (int)nodeattr[head];
    tailattr = (int)nodeattr[tail];
    if (headattr == tailattr) { /* They match; otherwise don't bother */
      echange=(EdgetreeSearch(head, tail, oe)==0)? 1:-1;
      headdeg=0;
      for(e = EdgetreeMinimum(oe, head);
      (tmp = oe[e].value) != 0;
      e = EdgetreeSuccessor(oe, e)) {
        headdeg += (nodeattr[tmp]==headattr);
      }
/*      for(e = EdgetreeMinimum(ie, head); */
/*      (tmp = ie[e].value) != 0; */
/*      e = EdgetreeSuccessor(ie, e)) { */
/*        headdeg += (nodeattr[tmp]==headattr); */
/*      } */
      for(j = 0; j < N_CHANGE_STATS; j++) {
        deg = (Vertex)mtp->inputparams[j];
        mtp->dstats[j] += (headdeg + echange == deg) - (headdeg == deg);
      }
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_ostar
*****************/
CHANGESTAT_FN(d_ostar) { 
  double change, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hattr = mtp->attrib[t-1];
      if(hattr == mtp->attrib[h-1]){
        td = - edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, h);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) { /* step through outedges of tail */
          if(hattr == mtp->attrib[node3-1]){++td;}
        }
        for(j=0; j < N_CHANGE_STATS; j++) {
          kmo = ((int)mtp->inputparams[j]) - 1;
          change = CHOOSE(td, kmo); 
          mtp->dstats[j] += (edgeflag ? - change : change); 
        }
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    for (i=0; i < ntoggles; i++) {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      td = nwp->outdegree[h] - edgeflag;      
      for(j=0; j < N_CHANGE_STATS; j++) {
        kmo = ((int)mtp->inputparams[j]) - 1;
        change = CHOOSE(td, kmo); 
        mtp->dstats[j] += (edgeflag ? - change : change); 
      }
      if (i+1 < ntoggles)
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  R    ***********/
/*****************
 changestat: d_receiver
*****************/
CHANGESTAT_FN(d_receiver) { 
  int i, j, echange;
  Vertex h, t, deg;
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
    if(t == 1){
      echange = -echange;
      for (j=0; j < N_CHANGE_STATS; j++){ 
        deg = (Vertex)mtp->inputparams[j];
        if(deg != 1)
          mtp->dstats[j] += echange;
      }
    }else{
      j=0;
      deg = (Vertex)mtp->inputparams[j];
      while(deg != t && j < N_CHANGE_STATS){
        j++;
        deg = (Vertex)mtp->inputparams[j];
      }
      if(j < N_CHANGE_STATS)
        mtp->dstats[j] += echange;
    }
    
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  S    ***********/
/*****************
 changestat: d_sender
*****************/
CHANGESTAT_FN(d_sender) { 
  int i, j, echange;
  Vertex h, t, deg;
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      if(h == 1){
       echange = -echange;
       for (j=0; j < N_CHANGE_STATS; j++){
         deg = (Vertex)mtp->inputparams[j];
         if(deg != 1){mtp->dstats[j] += echange;}
       }
      }else{
       j=0;
       deg = (Vertex)mtp->inputparams[j];
       while(deg != h && j < N_CHANGE_STATS){
	j++;
	deg = (Vertex)mtp->inputparams[j];
       }
       if(j < N_CHANGE_STATS){mtp->dstats[j] += echange;}
      }
      
      if (i+1 < ntoggles)
	TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_simmelian
*****************/
CHANGESTAT_FN(d_simmelian) { 
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i;
  
 *(mtp->dstats) = 0.0;
 for (i=0; i<ntoggles; i++) 
 {
  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
   
  if(EdgetreeSearch(t, h, nwp->outedges) != 0){
   change = 0;
   
   for(e = EdgetreeMinimum(nwp->outedges, t);
       (node3 = nwp->outedges[e].value) != 0;
       e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
   {
     if (node3 != h
      && EdgetreeSearch(node3, h, nwp->outedges) != 0 
      && EdgetreeSearch(h, node3, nwp->outedges) != 0 
      && EdgetreeSearch(node3, t, nwp->outedges) != 0 
        ){++change;}
   }
      
   *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
   }
   
   if (i+1 < ntoggles) 
     TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_simmelianties
*****************/
CHANGESTAT_FN(d_simmelianties) { 
  Edge e;
  Vertex h, t, change, node3, node4, first, firstht;
  int edgeflag, i;
  
 *(mtp->dstats) = 0.0;
 for (i=0; i<ntoggles; i++) 
 {
  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
   
  if(EdgetreeSearch(t, h, nwp->outedges) != 0){
   change = 0;
   firstht = 0;
   
   for(e = EdgetreeMinimum(nwp->outedges, t);
       (node3 = nwp->outedges[e].value) != 0;
       e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
   {
     if (node3 != h
      && EdgetreeSearch(node3, h, nwp->outedges) != 0 
      && EdgetreeSearch(h, node3, nwp->outedges) != 0 
      && EdgetreeSearch(node3, t, nwp->outedges) != 0 
        ){
          firstht = 1;
          ++change;
          first = 1;
          for(e = EdgetreeMinimum(nwp->outedges, h);
              (node4 = nwp->outedges[e].value) != 0;
              e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
          {
          if (node4 != t  && node4 != node3
           && EdgetreeSearch(node4, h, nwp->outedges) != 0 
           && EdgetreeSearch(node4, node3, nwp->outedges) != 0 
           && EdgetreeSearch(node3, node4, nwp->outedges) != 0 
             ){first = 0;}
	  }
	  if(first){++change;}

          first = 1;
          for(e = EdgetreeMinimum(nwp->outedges, t);
              (node4 = nwp->outedges[e].value) != 0;
              e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
          {
          if (node4 != h  && node4 != node3
           && EdgetreeSearch(node4, t, nwp->outedges) != 0 
           && EdgetreeSearch(node4, node3, nwp->outedges) != 0 
           && EdgetreeSearch(node3, node4, nwp->outedges) != 0 
             ){first = 0;}
	  }
	  if(first){++change;}
         }
   }
/*   if(firstht){++change;} */
      
   change = 2*change;
   *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
   }
   
   if (i+1 < ntoggles) 
     TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_smalldiff
*****************/
CHANGESTAT_FN(d_smalldiff) { 
  Vertex h, t;
  int i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) {
    h=heads[i];
    t=tails[i];
    *(mtp->dstats) += (fabs(mtp->attrib[h-1] - mtp->attrib[t-1])
    > *(mtp->inputparams)) ? 0.0 :
    ((EdgetreeSearch(h, t, nwp->outedges) != 0) ? -1.0 : 1.0); 
    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_sociality
*****************/
CHANGESTAT_FN(d_sociality) { 
  int i, j, echange;
  Vertex h, t, deg;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)N_CHANGE_STATS;
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i<ntoggles; i++)
      {      
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	hattr = mtp->attrib[h-1+nstats];
	if(hattr == mtp->attrib[t-1+nstats]){
	  j=0;
	  deg = (Vertex)mtp->inputparams[j];
	  while(deg != h && j < nstats){
	    j++;
	    deg = (Vertex)mtp->inputparams[j];
	  }
	  if(j < nstats){mtp->dstats[j] += echange;}
	  j=0;
	  deg = (Vertex)mtp->inputparams[j];
	  while(deg != t && j < nstats){
	    j++;
	    deg = (Vertex)mtp->inputparams[j];
	  }
	  if(j < nstats){mtp->dstats[j] += echange;}
	}
	
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i<ntoggles; i++)
      {      
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	j=0;
	deg = (Vertex)mtp->inputparams[j];
	while(deg != h && j < nstats){
	  j++;
	  deg = (Vertex)mtp->inputparams[j];
	}
	if(j < nstats){mtp->dstats[j] += echange;}
	j=0;
	deg = (Vertex)mtp->inputparams[j];
	while(deg != t && j < nstats){
	  j++;
	  deg = (Vertex)mtp->inputparams[j];
	}
	if(j < nstats){mtp->dstats[j] += echange;}
	
	if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/********************  changestats:  T    ***********/
/*****************
 changestat: d_tdsp
*****************/
CHANGESTAT_FN(d_tdsp) {
  Edge e, f;
  int i, j, echange, L2hu, L2ut;
  Vertex deg, h, t, u, v;
  
  for (i=0; i < N_CHANGE_STATS; i++) 
    CHANGE_STAT[i] = 0.0;

  for (i=0; i<ntoggles; i++){
    h = heads[i]; t=tails[i];
    echange = 1-2*IS_OUTEDGE(h,t);
    /* step through outedges of t */
    for(e = MIN_OUTEDGE(t); (u=OUTVAL(e))!=0; e=NEXT_OUTEDGE(e)) { 
      if (u != h){
        L2hu=0; /* This will be # of shared partners of (h,u) */
        /* step through inedges of u, incl. (t,u) itself */
        for(f = MIN_INEDGE(u); (v=INVAL(f))!=0; f=NEXT_INEDGE(f)) {
          if(IS_OUTEDGE(h,v)) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2hu + echange == deg) - (L2hu == deg));
        }
      }
    }
    /* step through inedges of h */
    for(e = MIN_INEDGE(h); (u=INVAL(e))!=0; e=NEXT_INEDGE(e)) {
      if (u != t){
        L2ut=0; /* This will be # of shared partners of (u,t) */
        /* step through outedges of u , incl. (u,h) itself */
        for(f = MIN_OUTEDGE(u);(v=OUTVAL(f))!=0; f=NEXT_OUTEDGE(f)){
          if(IS_OUTEDGE(v,t)) L2ut++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)INPUT_PARAM[j];
          CHANGE_STAT[j] += ((L2ut + echange == deg) - (L2ut == deg));
        }
      }
    }
    
    if (i+1 < ntoggles) TOGGLE(h,t);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_tesp
*****************/
CHANGESTAT_FN(d_tesp) { 
  Edge e, f;
  int i, j, echange;
  int L2ht, L2hu, L2ut;
  Vertex deg;
  Vertex h, t, u, v;
  TreeNode *oe=nwp->outedges, *ie=nwp->inedges;

  
  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++){      
    L2ht=0;
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(oe, t); 
    (u = oe[e].value) != 0; e = EdgetreeSuccessor(oe, e)) {
      if (EdgetreeSearch(h, u, oe) != 0){
        L2hu=0;
        /* step through inedges of u */
        for(f = EdgetreeMinimum(ie, u); (v = ie[f].value) != 0;
        f = EdgetreeSuccessor(ie, f)){
          if(EdgetreeSearch(h,v,oe)!= 0) L2hu++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2hu + echange == deg) - (L2hu == deg));
        }
      }
    }
    /* step through inedges of t */
    for (e = EdgetreeMinimum(ie, t); (u = ie[e].value) != 0;
    e = EdgetreeSuccessor(ie, e)){
      if (EdgetreeSearch(h, u, oe) != 0){
        L2ht++;
      }
      if (EdgetreeSearch(u, h, oe) != 0){
        L2ut=0;
        /* step through outedges of u */
        for(f = EdgetreeMinimum(oe, u); (v = oe[f].value) != 0;
        f = EdgetreeSuccessor(oe, f)){
          if(EdgetreeSearch(v, t,oe)!= 0) L2ut++;
        }
        for(j = 0; j < N_CHANGE_STATS; j++){
          deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += ((L2ut + echange == deg) - (L2ut == deg));
        }
      }
    }
    for(j = 0; j < N_CHANGE_STATS; j++){
      deg = (Vertex)mtp->inputparams[j];
      mtp->dstats[j] += echange*(L2ht == deg);
    }
    if (i+1 < ntoggles)
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_transitive
*****************/
CHANGESTAT_FN(d_transitive) { 
  Edge e;
  Vertex h, t, node2;
  double change;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
/*           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag); */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node2 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (node2 != h){
       if (EdgetreeSearch(h, node2, nwp->outedges) == 0){
        change = change - 1.0;
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, t);
	    (node2 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      if (node2 != h){
       if (EdgetreeSearch(h, node2, nwp->outedges) != 0){
        change = change + 1.0;
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, h);
	    (node2 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of hail */
    {
      if (node2 != t){
       if (EdgetreeSearch(node2, t, nwp->outedges) == 0){
        change = change - 1.0;
       }
      }
    }
    
    *(mtp->dstats) += edgeflag ? -change : change;
/*  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change); */

    if (i+1 < ntoggles) 
      TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_triadcensus
*****************/
CHANGESTAT_FN(d_triadcensus) { 
  int i, j, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex triadtype, node3, h, t;

  for (i=0; i < N_CHANGE_STATS; i++) 
    mtp->dstats[i] = 0.0;

  if (nwp->directed_flag) {
/* directed version */
   for (i=0; i<ntoggles; i++) 
    {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0;
      t210 = 0;
      t120C = 0;  t120U = 0;   t120D = 0;  t201 = 0;
      t030C = 0;  t030T = 0;   t111U = 0;  t111D = 0;
      t021C = 0;  t021U = 0;   t021D = 0;  t102 = 0;
      t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++)
             { 
             if (node3 != h && node3 != t)
               {   
               a = (EdgetreeSearch(t, h, nwp->outedges) != 0); 
               b = (EdgetreeSearch(t, node3, nwp->outedges) != 0);
               c = (EdgetreeSearch(node3, t, nwp->outedges) != 0);
               d = (EdgetreeSearch(node3, h, nwp->outedges) != 0);
               e = (EdgetreeSearch(h, node3, nwp->outedges) != 0);
               edgecount = (a + b + c + d + e);

               switch(edgecount)
               {  
               case 0:   /* 012 */
                 ++t012;

               case 1:   /* 021C, 021U, 021D, 102 */
                 {
                     if ((b == 1) || (d == 1))
                       ++t021C;
                     if (c == 1)
                       ++t021U;
                     if (e == 1)
                       ++t021D;
                     if (a == 1)
                       ++t102;
                 }
                 break;
          
               case 2:   /* 030C, 030T, 111U, 111D */
                 {
                     if ((b + d) == 2)       
                       ++t030C;
                     if (((b + e) == 2) || ((c + d) == 2) || ((c + e) == 2))
                       ++t030T;
                     if (((a + b) == 2) || ((a + e) == 2) || ((d + e) == 2))
                       ++t111U;
                     if (((a + c) == 2) || ((a + d) == 2) || ((b + c) == 2))
                       ++t111D;
                 }
                 break;
            
               case 3:   /* 120C, 120U, 120D, 201 */
                 {
                     if (a == 1)
                         {
                            if (((b + d) == 2) || ((c + e) == 2))
                              ++t120C;
                            if ((b + e) == 2)
                              ++t120U;
                            if ((c + d) == 2)
                              ++t120D;
                            if (((b + c) == 2) || ((d + e) == 2))
                              ++t201;
                         }
                     else 
                         {
                            if (b == 1)
                                {
                                    if (((c + d) == 2) || ((d + e) == 2))
                                      ++t120C;
                                    if ((c + e) == 2)
                                      ++t120D;
                                } 
                            else 
                                {
                                    ++t120U;
                                }
                         } 
                 }
                 break;
             
               case 4:   /* 210 */
                 ++t210;
                 break;
             
               case 5:   /* 300 */            
                 ++t300;
                 break;
               }

               switch(edgecount)
               {  
               case 1:   /* 102, 021D, 021U, 021C */
                 --t012;
                 break;
           
               case 2:   /* 030C, 030T, 111U, 111D */ 
                 {
                     if (((a + c) == 2) || ((a + e) == 2) || ((b + d) == 2) || 
                        ((c + e) == 2)) 
                       --t021C;
                     if (((a + d) == 2) || ((b + e) == 2))
                       --t021U;
                     if (((a + b) == 2) || ((c + d) == 2))
                       --t021D;
                     if (((b + c) == 2) || ((d + e) == 2))
                       --t102;
                 } 
                 break;

               case 3:   /* 201, 120D, 120U, 120C */
                 {
                     if (a == 1)
	    	           {
	    	           if ((c + e) == 2)       
                             --t030C;
                           if (((c + d) == 2) || ((b + e) == 2) || ((b + d) == 2))
                             --t030T;
                           if ((b + c) == 2)
                             --t111U;
                           if ((d + e) == 2)
	                      --t111D;
	 	               }
 		             else
 		               {
 		                   if (b == 1)
 		 	                 {
 		 	                   if ((c + d) == 2)
 			                     --t111U;
 			                   if (((c + e) == 2) || ((d + e) == 2))
			                     --t111D;
 		 	                 }
    	    	                   else
		    	                 --t111U;
    		               }
                   }
                   break;
          
               case 4:   /* 210 */
                 {
                     if (a == 1)
                         {
                            if (((b + c + e) == 3) || ((c + d + e) == 3))
                              --t120C;
                            if ((b + c + d) == 3)
                              --t120U;
                            if ((b + d + e) == 3)
                              --t120D;
                         }
                     else 
                         {
                            if ((b + c + d + e) == 4)
                              --t201;
                         } 
                 }
                 break;
              
               case 5:   /* 300 */            
                 --t210;
                 break;
               }
               }
             }    /* ******  move to next node3 ******** */
        }
        else 
          t012 = t012 + (N_NODES - 2);  

      for(j = 0; j < N_CHANGE_STATS; j++)
        { 
	    triadtype = (Vertex)mtp->inputparams[j]; 

            switch(triadtype)
            {
                case 1:  t003 = (t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
			 t003 = t003+(t111U+t111D+t021C+t021U+t021D+t102+t012);
 	                 mtp->dstats[j] += edgeflag ? -(double)t003 : (double)t003;
		      break;
 	        case 2:   mtp->dstats[j] += edgeflag ? -(double)t012 : (double)t012;
	              break;
	        case 3:   mtp->dstats[j] += edgeflag ? -(double)t102 : (double)t102;
                      break;
	        case 4:   mtp->dstats[j] += edgeflag ? -(double)t021D : (double)t021D;
                      break;
	        case 5:   mtp->dstats[j] += edgeflag ? -(double)t021U : (double)t021U;
                      break;
 	        case 6:   mtp->dstats[j] += edgeflag ? -(double)t021C : (double)t021C;
                      break;
	        case 7:	  mtp->dstats[j] += edgeflag ? -(double)t111D : (double)t111D;
                      break;
	        case 8:	  mtp->dstats[j] += edgeflag ? -(double)t111U : (double)t111U;
                      break;
	        case 9:	  mtp->dstats[j] += edgeflag ? -(double)t030T : (double)t030T;
                      break;
	        case 10:   mtp->dstats[j] += edgeflag ? -(double)t030C : (double)t030C;
                      break;
	        case 11:  mtp->dstats[j] += edgeflag ? -(double)t201 : (double)t201;
                      break;
	        case 12:  mtp->dstats[j] += edgeflag ? -(double)t120D : (double)t120D;
                      break;
	        case 13:  mtp->dstats[j] += edgeflag ? -(double)t120U : (double)t120U;
                      break;
	        case 14:  mtp->dstats[j] += edgeflag ? -(double)t120C : (double)t120C;
                      break;
	        case 15:  mtp->dstats[j] += edgeflag ? -(double)t210 : (double)t210;
                      break;
  	        case 16:  mtp->dstats[j] += edgeflag ? -(double)t300 : (double)t300;
		      break;
            }
        }
      if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */ 
    }
    }else{
/*  undirected */
    for (i=0; i<ntoggles; i++) 
     {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= N_NODES; node3++)
             { 
             if (node3 != h && node3 != t)
               {   
               a = (EdgetreeSearch(MIN(node3,t), MAX(node3,t), nwp->outedges) != 0);
               b = (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0);
               edgecount = (a + b);

               switch(edgecount)
               {  
               case 0:   /* 012 */
		{
                 ++t102;
                 --t012;
		}
                break;

               case 1:   /* 021C, 021U, 021D, 102 */
		{
                 ++t201;
                 --t102;
		}
                break;
          
               case 2:   /* 030C, 030T, 111U, 111D */
		{
                 ++t300;
                 --t201;
		}
                break;
            
               }
	       }

             }    /* ******  move to next node3 ******** */
        }
        else 
          t102 = t102 + (N_NODES - 2);  

      for(j = 0; j < N_CHANGE_STATS; j++)
        { 
	    triadtype = (Vertex)mtp->inputparams[j]; 

            switch(triadtype)
            {
                case 1:  t003 = (t102+t201+t300);
 	                 mtp->dstats[j] += edgeflag ? -(double)t003 : (double)t003;
		      break;
 	        case 2:  mtp->dstats[j] += edgeflag ? -(double)t102 : (double)t102;
		      break;
	        case 3:  mtp->dstats[j] += edgeflag ? -(double)t201 : (double)t201;
                      break;
	        case 4:  mtp->dstats[j] += edgeflag ? -(double)t300 : (double)t300;
                      break;
            }
        }
     if (i+1 < ntoggles)
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */ 
     } /* i loop */
    }
    
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_triangle
*****************/
CHANGESTAT_FN(d_triangle) { 
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr, echange;
  
  ninputs = mtp->ninputparams;
  nstats = N_CHANGE_STATS;
  
  if(ninputs>0){
    /* match on attributes */
    if(nstats>1){
      for (j=0; j<nstats; j++){
	mtp->dstats[j] = 0.0;
      }
    }else{
      *(mtp->dstats) = 0.0;
    }
    for (i=0; i<ntoggles; i++) 
      {
	hattr = mtp->attrib[heads[i]-1];
	if(hattr == mtp->attrib[tails[i]-1]){
	  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	  change = 0;
	  for(e = EdgetreeMinimum(nwp->outedges, t);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	    {
	      if(hattr == mtp->attrib[node3-1]){
		if (nwp->directed_flag)
		  {
		    if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
		      ++change;
		    if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
		      ++change;
		  }
		else
		  {
		    if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		      ++change;
		    }
		  }
	      }
	    }
	  
	  for(e = EdgetreeMinimum(nwp->inedges, t); 
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	    {
	      if(hattr == mtp->attrib[node3-1]){
		if (nwp->directed_flag)
		  {
		    if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
		      ++change;
		    if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
		      ++change;
		  }
		else
		  {
		    if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		      ++change;
		    }
		  }
	      }
	    }
	  
	  echange =  edgeflag ? -(double)change : change;
	  if(nstats>1){
	    for (j=0; j<nstats; j++){
	      mtp->dstats[j] += ((hattr==mtp->inputparams[j]) ? echange : 0.0); 
	    }
	  }else{
	    *(mtp->dstats) += echange;
	  }
	  
	}
	
	if (i+1 < ntoggles) 
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }else{
    /* no attribute matching */
    *(mtp->dstats) = 0.0;
    for (i=0; i<ntoggles; i++) 
      {
	edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	
	change = 0;
	
	for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	  {
	    if (nwp->directed_flag)
	      {
		if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
		  ++change;
		if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
		  ++change;
	      }
	    else
	      {
		if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		  ++change;
		}
	      }
	  }
	
	for(e = EdgetreeMinimum(nwp->inedges, t); 
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	  {
	    if (nwp->directed_flag)
	      {
		if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
		  ++change;
		if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
		  ++change;
	      }
	    else
	      {
		if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0){
		  ++change;
		}
	      }
	  }
	
	*(mtp->dstats) += edgeflag ? -(double)change : change;
	
	if (i+1 < ntoggles) 
	  TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_tripercent
*****************/
CHANGESTAT_FN(d_tripercent) { 
  Edge e;
  Vertex h, t, change, node3;
  double hd, td=0.0, numchange;
  double tripercent=0.0, newtripercent=0.0;
  int edgeflag, i, j;
  int nedges, nnodes, nstats;
  Vertex *head, *tail;
  double *num2star, *numtri;
  double *newnumtri, *newnum2star;
  int ninputs;
  double hattr, eps=0.00000001;
  head = (Vertex *) malloc(sizeof(Vertex) * nwp->nedges);
  tail = (Vertex *) malloc(sizeof(Vertex) * nwp->nedges);
  
  nstats = N_CHANGE_STATS;
  ninputs = mtp->ninputparams;
  nnodes = N_NODES;
  
  num2star = (double *) malloc(sizeof(double) * nstats);
  numtri = (double *) malloc(sizeof(double) * nstats);
  newnum2star = (double *) malloc(sizeof(double) * nstats);
  newnumtri = (double *) malloc(sizeof(double) * nstats);
  for (i=0; i<nstats; i++){
    num2star[i] = 0.0;
    numtri[i] = 0.0;
  }
  
  nedges=0;
  if(ninputs>0){
    /* match on attributes */
    for (h=1; h<=nnodes; h++) 
    {
      for(e = EdgetreeMinimum(nwp->outedges, h);
	    (t = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
      {
        if(fabs(mtp->attrib[h-1] - mtp->attrib[t-1])<eps){
          head[nedges] = h;
          tail[nedges] = t;
          nedges++;
        }
      }
    }
  }else{
    /* no attribute matching */
    for (h=1; h<=nnodes; h++) 
    {
      for(e = EdgetreeMinimum(nwp->outedges, h);
	    (t = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
      {
        head[nedges] = h;
        tail[nedges] = t;
        nedges++;
      }
    }
  }
  if(ninputs>0){
    /* match on attributes */
    for (i=0; i<nedges; i++) {
      h=head[i];
      t=tail[i];
      hattr = mtp->attrib[h-1];
      hd = -1;
      for(e = EdgetreeMinimum(nwp->outedges, h);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
      {
        if(fabs(hattr - mtp->attrib[node3-1])<eps){++hd;}
      }
      if (!nwp->directed_flag){
        for(e = EdgetreeMinimum(nwp->inedges, h);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
        {
          if(hattr == mtp->attrib[node3-1]){++hd;}
        }
        td = -1;
        for(e = EdgetreeMinimum(nwp->outedges, t);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
        {
          if(hattr == mtp->attrib[node3-1]){++td;}
        }
        for(e = EdgetreeMinimum(nwp->inedges, t);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
        {
          if(hattr == mtp->attrib[node3-1]){++td;}
        }
      }
      /* calculate the change in the number of 2-stars */
      if (nwp->directed_flag){
        numchange = hd;
      }else{
        numchange = hd + td;
      }
      if(nstats>1){
        for (j=0; j<nstats; j++){
          num2star[j] += (hattr==mtp->inputparams[j]) ? numchange : 0.0;
        }
      }else{
        num2star[0] += numchange;
      }
      
      change = 0;
      
      for(e = EdgetreeMinimum(nwp->outedges, t);
      (node3 = nwp->outedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        if(fabs(hattr - mtp->attrib[node3-1])<eps){
          if (nwp->directed_flag)
          {
            if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
              ++change;
            if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
              ++change;
          }
          else
          {
            if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
              ++change;
          }
        }
      }
      
      for(e = EdgetreeMinimum(nwp->inedges, t); 
      (node3 = nwp->inedges[e].value) != 0;
      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        if(fabs(hattr - mtp->attrib[node3-1])<eps){
          if (nwp->directed_flag)
          {
            if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
              ++change;
            if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
              ++change;
          }
          else
          {
            if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
              ++change;
          }
        }
      }
      
      /* diff=T (and more than one category?)  */
      if(nstats>1){
        for (j=0; j<nstats; j++){
          numtri[j] += (hattr==mtp->inputparams[j]) ? (double)change : 0.0;
        }
      }else{
        numtri[0] += (double)change;
      }
      
      if (i+1 < nedges) 
        ToggleEdge(head[i], tail[i], nwp);  /* Toggle this edge if more to come */
    }
  }else{
    /* no attribute matching */
    for (i=0; i<nedges; i++) 
    {
      h=head[i];
      t=tail[i];
      /* calculate the change in the number of 2-stars */
      if (nwp->directed_flag){
        hd = nwp->outdegree[h] - 1; 
        numchange = hd;
      }else{
        hd = nwp->outdegree[h] + nwp->indegree[h] - 1;
        td = nwp->outdegree[t] + nwp->indegree[t] - 1;
        numchange = hd + td;
      }
      num2star[0] += numchange;
      change = 0;
      for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        if (nwp->directed_flag)
	      {
          if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
            ++change;
          if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
            ++change;
	      }
        else
	      {
          if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
            ++change;
	      }
      }
      
      for(e = EdgetreeMinimum(nwp->inedges, t); 
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        if (nwp->directed_flag)
	      {
          if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
            ++change;
          if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
            ++change;
	      }
        else
	      {
          if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
            ++change;
	      }
      }
      numtri[0] += (double)change;
      if (i+1 < nedges) 
        ToggleEdge(head[i], tail[i], nwp);  /* Toggle this edge if more to come */
    }
  }
  
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(head[i], tail[i], nwp); 
  
  free(head);
  free(tail);
  
  for (j=0; j<nstats; j++){
    newnumtri[j]=numtri[j];
    newnum2star[j]=num2star[j];
  }
  
  if(ninputs>0){
    /* match on attributes */
    for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hattr = mtp->attrib[h-1];
      if(fabs(hattr - mtp->attrib[t-1])<eps){
        hd = -edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, h);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
        {
          if(fabs(hattr - mtp->attrib[node3-1])<eps){++hd;}
        }
        if (!nwp->directed_flag){
          for(e = EdgetreeMinimum(nwp->inedges, h);
          (node3 = nwp->inedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
          {
            if(fabs(hattr - mtp->attrib[node3-1])<eps){++hd;}
          }
          td = - edgeflag;
          for(e = EdgetreeMinimum(nwp->outedges, t);
          (node3 = nwp->outedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
          {
            if(fabs(hattr - mtp->attrib[node3-1])<eps){++td;}
          }
          for(e = EdgetreeMinimum(nwp->inedges, t);
          (node3 = nwp->inedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
          {
            if(fabs(hattr - mtp->attrib[node3-1])<eps){++td;}
          }
        }
        /* calculate the change in the number of 2-stars */
        if (nwp->directed_flag){
          numchange = hd;
        }else{
          numchange = hd + td;
        }
        /* diff=T (and more than one category?)  */
        if(nstats>1){
          for (j=0; j<nstats; j++){
            newnum2star[j] += (hattr==mtp->inputparams[j]) ? 
            (edgeflag ? - numchange : numchange) : 0.0;
          }
        }else{
          newnum2star[0] += (edgeflag ? - numchange : numchange);
        }
        
        change = 0;
        
        for(e = EdgetreeMinimum(nwp->outedges, t);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
        {
          if(fabs(hattr - mtp->attrib[node3-1])<eps){
            if (nwp->directed_flag)
            {
              if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
                ++change;
              if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
                ++change;
            }
            else
            {
              if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
                ++change;
            }
          }
        }
        
        for(e = EdgetreeMinimum(nwp->inedges, t); 
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
        {
          if(fabs(hattr - mtp->attrib[node3-1])<eps){
            if (nwp->directed_flag)
            {
              if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
                ++change;
              if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
                ++change;
            }
            else
            {
              if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
                ++change;
            }
          }
        }
        /* diff=T (and more than one category?)  */
        if(nstats>1){
          for (j=0; j<nstats; j++){
            newnumtri[j] += (hattr==mtp->inputparams[j]) ? 
            (edgeflag ? -(double)change : change) : 0.0;
          }
        }else{
          newnumtri[0] += (edgeflag ? -(double)change : change);
        }
      }
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    /* no attribute matching */
    for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /* calculate the change in the number of 2-stars */
      if (nwp->directed_flag){
        hd = nwp->outdegree[h] - edgeflag; 
        numchange = hd;
      }else{
        hd = nwp->outdegree[h] + nwp->indegree[h] - edgeflag;
        td = nwp->outdegree[t] + nwp->indegree[t] - edgeflag;
        numchange = hd + td;
      }
      
      newnum2star[0] += (edgeflag ? - numchange : numchange);
      
      change = 0;
      
      for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        if (nwp->directed_flag)
	      {
          if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
            ++change;
          if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
            ++change;
	      }
        else
	      {
          if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
            ++change;
	      }
      }
      
      for(e = EdgetreeMinimum(nwp->inedges, t); 
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        if (nwp->directed_flag)
	      {
          if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
            ++change;
          if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
            ++change;
	      }
        else
	      {
          if (EdgetreeSearch(MIN(node3,h), MAX(node3,h), nwp->outedges) != 0)
            ++change;
	      }
      }
      
      newnumtri[0] += edgeflag ? -(double)change : change;
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  for (j=0; j<nstats; j++){
    tripercent = 0.0;
    newtripercent = 0.0;
    if(num2star[j]>0.0){      tripercent =    numtri[j]/num2star[j];}
    if(newnum2star[j]>0.0){newtripercent = newnumtri[j]/newnum2star[j];}
    if(nstats>1){
      mtp->dstats[j] = (newtripercent - tripercent)*300.0;
    }else{
      *(mtp->dstats) = (newtripercent - tripercent)*300.0;
    }
  }
  free(num2star);
  free(numtri);
  free(newnum2star);
  free(newnumtri);
                                        
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}

/*****************
 changestat: d_ttriple
*****************/
CHANGESTAT_FN(d_ttriple) { 
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr;
  
  ninputs = mtp->ninputparams;
  nstats = N_CHANGE_STATS;
  if(ninputs>0){
    /* match on attributes */
    if(nstats>1){
      for (j=0; j<nstats; j++){
        mtp->dstats[j] = 0.0;
      }
    }else{
      *(mtp->dstats) = 0.0;
    }
    for (i=0; i<ntoggles; i++) 
    {
      hattr = mtp->attrib[heads[i]-1];
      if(hattr == mtp->attrib[tails[i]-1]){
        edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
        change = 0;
        
        for(e = EdgetreeMinimum(nwp->outedges, t);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
        {
          if(hattr == mtp->attrib[node3-1]){
            if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
              ++change;
          }
        }
        
        for(e = EdgetreeMinimum(nwp->inedges, t); 
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
        {
          if(hattr == mtp->attrib[node3-1]){
            if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
              ++change;
            if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
              ++change;
          }
        }
        
        if(nstats>1){
          for (j=0; j<nstats; j++){
            mtp->dstats[j] += (hattr==mtp->inputparams[j]) ? 
            (edgeflag ? -(double)change : (double)change) : 0.0;
          }
        }else{
          *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
        }
        
      }
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }else{
    /* no attribute matching */
    *(mtp->dstats) = 0.0;
    for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      
      change = 0;
      
      for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
      {
        if (EdgetreeSearch(node3, h, nwp->inedges) != 0)
          ++change;
      }
      
      for(e = EdgetreeMinimum(nwp->inedges, t); 
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
      {
        if (EdgetreeSearch(node3, h, nwp->outedges) != 0)
          ++change;
        if (EdgetreeSearch(node3, h,  nwp->inedges) != 0)
          ++change;
      }
      
      *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
      
      if (i+1 < ntoggles) 
        TOGGLE(heads[i], tails[i]);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    TOGGLE(heads[i], tails[i]); 
}



