#include "basechangeStats.h"


/*****************
 double my_choose

 Simple routine to return simple binomial coefficients quickly, 
 avoiding costly call to choose() function.  Note:  my_choose is
 usually not called directly; use CHOOSE macro instead.
*****************/

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
 void d_absdiff

*****************/

void d_absdiff (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp) {
  double change;
  Vertex h, t;
  int i, edgeflag;

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      edgeflag = (EdgetreeSearch(h=heads[i],t=tails[i],nwp->outedges) != 0);
      change = abs(mtp->attrib[h-1] - mtp->attrib[t-1]);
      *(mtp->dstats) += edgeflag ? -change : change;
      if (i+1 < ntoggles) 
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}
  

/*****************
 void d_absdiffcat
*****************/
void d_absdiffcat (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  double change, NAsubstitute, hval, tval;
  Vertex h, t, ninputs;
  int i, j, edgeflag=0;
  /* double checksum;*/ 
  
  ninputs = mtp->ninputparams - nwp->nnodes;
  NAsubstitute = mtp->inputparams[ninputs-1];
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    edgeflag = (EdgetreeSearch(h=heads[i],t=tails[i],nwp->outedges) != 0);
    hval = mtp->attrib[h-1];
    tval = mtp->attrib[t-1];
    if (hval == NAsubstitute ||  tval == NAsubstitute) change = NAsubstitute;
    else change = abs(hval - tval);
	  if (change>0) {
      for (/*checksum=0.0,*/ j=0; j<ninputs; j++) {
        mtp->dstats[j] += (change==mtp->inputparams[j]) ? 
             (edgeflag ? -1.0 : 1.0) : 0.0;
        /*checksum += (change==mtp->inputparams[j]) ? 
             (edgeflag ? -1.0 : 1.0) : 0.0; */
      }
      /*if (abs(checksum) != 1.0) Rprintf("Whoops! checksum=%f\n",checksum);*/
    }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedkstar

*****************/
void d_boundedkstar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  double change, hod, tod;
  double newhod, newtod;
  int edgeflag, i, j, k, bound;
  int p = mtp->nstats;
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
	ToggleEdge(heads[i], tails[i], nwp) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedistar
*****************/
void d_boundedistar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  double change, tod;
  double newtod;
  int edgeflag, i, j, k, bound;
  int p = mtp->nstats;
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
	ToggleEdge(heads[i], tails[i], nwp) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedostar

*****************/
void d_boundedostar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  double change, hod;
  double newhod;
  int edgeflag, i, j, k, bound;
  int p = mtp->nstats;
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
	ToggleEdge(heads[i], tails[i], nwp) ;    /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedtriangle

*****************/
void d_boundedtriangle (int ntoggles, Vertex *heads, Vertex *tails, 
			ModelTerm *mtp, Network *nwp) {
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp);
}

/*****************
 void d_triangle

*****************/
void d_triangle (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr, echange;
  
  ninputs = mtp->ninputparams;
  nstats = mtp->nstats;
  
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
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_degree

*****************/
void d_degree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  int i, j, echange;
  Vertex h, t, hd, td=0, deg, *id, *od;
  TreeNode *oe;  

  oe=nwp->outedges;
  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
      hd = od[h] + id[h];
      td = od[t] + id[t];
      for(j = 0; j < mtp->nstats; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
          mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
          mtp->dstats[j] += (td + echange == deg) - (td == deg);
	}
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_edges
*****************/
void d_edges(int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int edgeflag, i;
  Vertex h, t;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      *(mtp->dstats) += edgeflag ? - 1 : 1;
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_density
*****************/
void d_density(int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int edgeflag, i;
  Vertex h, t, ndyads;
  
  ndyads = (nwp->nnodes)*(nwp->nnodes-1);
  if(!nwp->directed_flag){
    ndyads = ndyads / 2;
  }

  *(mtp->dstats) = 0.0;
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      *(mtp->dstats) += edgeflag ? - 1 : 1;
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = *(mtp->dstats) / ndyads;
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_isolates

*****************/
void d_isolates (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  int i, echange;
  Vertex h, t, hd, td=0, *id, *od;
  TreeNode *oe;  

  oe=nwp->outedges;
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}



/*****************
 void d_geodegree
*****************/
void d_geodegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int i, echange;
  double alpha;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = 
	(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      *(mtp->dstats) += echange*(exp(-alpha*hd)+exp(-alpha*td));
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  *(mtp->dstats) = (*(mtp->dstats))*(exp(-alpha)-1.0);
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwd
*****************/
void d_gwd (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp) {
  int i, echange;
  double alpha;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = 
	(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      *(mtp->dstats) += echange*(exp(-alpha*hd)+exp(-alpha*td)); 
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_geospartner
*****************/
void d_geospartner (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp) {
  Edge e, f;
  int i, echange;
  int L2ht, L2hu, L2ut;
  Vertex h, t, u, v;
  double alpha, cumchange;
  
  *(mtp->dstats) = 0.0;
  alpha = mtp->inputparams[0];
  
  for (i=0; i<ntoggles; i++){      
    cumchange=0.0;
    L2ht=0;
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=0;
	L2ut=0;
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
	cumchange += exp(-alpha*L2hu)+exp(-alpha*L2ut);
      }
    }
    /* step through inedges of t */
    for(e = EdgetreeMinimum(nwp->inedges, t);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=0;
	L2ut=0;
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
	cumchange += exp(-alpha*L2hu)+exp(-alpha*L2ut);
      }
    }
    cumchange  = cumchange*(exp(-alpha*echange)-1.0);
    cumchange += echange*(exp(-alpha*L2ht)-1.0);
    *(mtp->dstats) -= cumchange;
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwdegreealpha
*****************/
void d_gwdegreealpha (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  change = 0.0;
  alpha = mtp->inputparams[0];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      if(hd!=0){
        change += echange*(1.0-pow(oneexpa,(double)hd));
      }
      if(td!=0){
        change += echange*(1.0-pow(oneexpa,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*exp(alpha);
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwdegreelambda
*****************/
void d_gwdegreelambda (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  lambda = mtp->inputparams[0];
  oneexpl = 1.0-1.0/lambda;
  
  change = 0.0;
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      change += 4.0*echange;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      if(hd!=0){
        change -= echange*(1.0-pow(oneexpl,(double)hd));
      }
      if(td!=0){
        change -= echange*(1.0-pow(oneexpl,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwdegree
*****************/
void d_gwdegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      change += 4.0*echange;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      if(hd!=0){
        change -= echange*(1.0-pow(oneexpd,(double)hd));
      }
      if(td!=0){
        change -= echange*(1.0-pow(oneexpd,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwodegree
*****************/
void d_gwodegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex h, t, hd=0, *od;
  
  od=nwp->outdegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      change += 2.0*echange;
      hd = od[h] + (echange - 1)/2;
      if(hd!=0){
        change -= echange*(1.0-pow(oneexpd,(double)hd));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwidegree
*****************/
void d_gwidegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double decay, oneexpd, change;
  Vertex h, t, td=0, *id;
  
  id=nwp->indegree;
  decay = mtp->inputparams[0];
  oneexpd = 1.0-exp(-decay);
  
  change = 0.0;
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      change += 2.0*echange;
      td = id[t] + (echange - 1)/2;
      if(td!=0){
        change -= echange*(1.0-pow(oneexpd,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_icvar
*****************/
void d_icvar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, edgeflag, ichange, change;
  Vertex nnodes, h, t, *id;
  
  id=nwp->indegree;
  nnodes = nwp->nnodes;
  
  change = 0;
  for (i=0; i<ntoggles; i++)
    {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      if(edgeflag){
        ichange = -(2*(nnodes*(id[t]-1) - nwp->nedges+1) + nnodes - 1);
      }else{
        ichange =   2*(nnodes* id[t]    - nwp->nedges  ) + nnodes - 1;
      }
//   Rprintf("h %d t %d nnodes %d  nwp->nedges %d id[t] %d  ic %d\n",h,t, nnodes,  nwp->nedges, id[t], ichange);
//      change += edgeflag ? (-ichange) : ichange;
      change += ichange;
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*1.0/(nnodes*(nnodes-1));
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_idc
*****************/
void d_idc (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp)  {
  int i, edgeflag, ichange, change;
  Vertex k, nnodes, maxidegree0, maxidegree1;
  Vertex h, t, *id;
  
  id=nwp->indegree;
  nnodes = nwp->nnodes;
  
  change = 0;
  for (i=0; i<ntoggles; i++)
    {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      if(edgeflag){
	maxidegree0 = id[t];
	maxidegree1 = id[t]-1;
        for (k=1; k<=nnodes; k++){
	 if(id[k] > maxidegree0) {maxidegree0=id[k];}
	 if(k != t && id[k] > maxidegree1) {maxidegree1=id[k];}
	}
        ichange = nnodes*(maxidegree1-maxidegree0) + 1;
      }else{
	maxidegree0 = 0;
	maxidegree1 = id[t]+1;
        for (k=1; k<=nnodes; k++){
	 if(id[k] > maxidegree0) {maxidegree0=id[k];}
	 if(id[k] > maxidegree1) {maxidegree1=id[k];}
	}
        ichange = nnodes*(maxidegree1-maxidegree0) - 1;
      }
//   Rprintf("h %d t %d nnodes %d  nwp->nedges %d id[t] %d  ic %d\n",h,t, nnodes,  nwp->nedges, id[t], ichange);
      change += ichange;
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*1.0/((nnodes-1)*(nnodes-1));
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_altkstar
*****************/
void d_altkstar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex h, t, hd, td=0, *id, *od;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  change = 0.0;
  lambda = mtp->inputparams[0];
  oneexpl = 1.0-1.0/lambda;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + id[h] + (echange - 1)/2;
      td = od[t] + id[t] + (echange - 1)/2;
      if(hd!=0){
        change += echange*(1.0-pow(oneexpl,(double)hd));
      }
      if(td!=0){
        change += echange*(1.0-pow(oneexpl,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*lambda;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_altostar
*****************/
void d_altostar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex h, t, hd=0, *od;
  
  od=nwp->outdegree;
  change = 0.0;
  lambda = mtp->inputparams[0];
  oneexpl = 1.0-1.0/lambda;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + (echange - 1)/2;
      if(hd!=0){
        change += echange*(1.0-pow(oneexpl,(double)hd));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*lambda;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_altistar
*****************/
void d_altistar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, echange=0;
  double lambda, oneexpl, change;
  Vertex h, t, td=0, *id;
  
  id=nwp->indegree;
  change = 0.0;
  lambda = mtp->inputparams[0];
  oneexpl = 1.0-1.0/lambda;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      td = id[t] + (echange - 1)/2;
      if(td!=0){
        change += echange*(1.0-pow(oneexpl,(double)td));
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  *(mtp->dstats) = change*lambda;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwesp
*****************/
void d_gwesp (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp)  {
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
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_gwdsp
****************/
void d_gwdsp (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
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
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_dyadcov (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  double val;
  Vertex h, t;
  int i, edgeflag, refedgeflag;
  long int nactors, nrows;
  
  nrows = (long int)(mtp->inputparams[0]);
  nactors = nwp->bipartite;
  
  if(!nwp->directed_flag){

  for(i=0;i<3;i++)
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial state of the edge and its reflection*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      refedgeflag = (EdgetreeSearch(t, h, nwp->outedges) != 0);
      
      /*Get the dyadic covariate*/
      val = mtp->attrib[(t-1-nactors)+(h-1)*nrows];
      
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
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }else{

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Get the covariate value*/
      val = mtp->attrib[(t-1-nactors)+(h-1)*nrows];
      /*Update the change statistic, based on the toggle type*/
      *(mtp->dstats) += edgeflag ? -val : val;
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_simmeliandynamic (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  Edge e;
  long int nnodes;
  Vertex h, t, node3, change;
  int i, edgeflag, edgeflagth;
  
  nnodes = ((long int)(mtp->inputparams[0]));
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      edgeflagth = (EdgetreeSearch(t, h, nwp->outedges) == 0);

      if(!edgeflagth){
       /*Check to see if this will form a Simmelian */
       change = 0;
   
       for(e = EdgetreeMinimum(nwp->outedges, t);
           (node3 = nwp->outedges[e].value) != 0;
           e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
       {
         if (EdgetreeSearch(node3, h, nwp->outedges) != 0 
          && EdgetreeSearch(h, node3, nwp->outedges) != 0 
          && EdgetreeSearch(node3, t, nwp->outedges) != 0 
            ){
          /*So this will form a h,t,node3 Simmelian */
          /*Check to see if a prior Simmelian existed */
          if( mtp->inputparams[1+(t-1)+(node3-1)*nnodes]==1 &&
            mtp->inputparams[1+(h-1)+(node3-1)*nnodes]==1 &&
            mtp->inputparams[1+(node3-1)+(h-1)*nnodes]==1 &&
            mtp->inputparams[1+(t-1)+(h-1)*nnodes]==1 &&
            mtp->inputparams[1+(h-1)+(t-1)*nnodes]==1 &&
            mtp->inputparams[1+(node3-1)+(t-1)*nnodes]==1 
	    ){
		++change;
	     }
	 }
       }
      
       change = 6*change;
       *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
      }
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_intransitivedynamic (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  Edge e;
  long int nnodes;
  Vertex h, t, node3;
  double change;
  int i, edgeflag;
  
  nnodes = ((long int)(mtp->inputparams[0]));
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (node3 != h){
       if (EdgetreeSearch(h, node3, nwp->outedges) == 0){
       /* So form a h -> t -> node3 intransitive triad */
       /* Check to see if a prior intransitive triad existed */
         if( mtp->inputparams[1+(h-1)+(t-1)*nnodes]==1 &&
          mtp->inputparams[1+(t-1)+(node3-1)*nnodes]==1 &&
          mtp->inputparams[1+(h-1)+(node3-1)*nnodes]==0
	 ){
//           Rprintf("h %d t %d node3 %d nnodes %d\n",h,t, node3, nnodes);
          change = change - 1.0;
	 }
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, t);
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      if (node3 != h){
       if (EdgetreeSearch(h, node3, nwp->outedges) != 0){
       /* So dissolve a h -> node3 -> t intransitive triad */
       /*Check to see if a prior intransitive triad existed */
         if( mtp->inputparams[1+(h-1)+(t-1)*nnodes]==0 &&
          mtp->inputparams[1+(h-1)+(node3-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(t-1)*nnodes]==1
	 ){
          change = change + 1.0;
	 }
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, h);
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of hail */
    {
      if (node3 != t){
       if (EdgetreeSearch(node3, t, nwp->outedges) == 0){
       /* So form a node3 -> h -> t intransitive triad */
       /*Check to see if a prior intransitive triad existed */
         if( mtp->inputparams[1+(h-1)+(t-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(h-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(t-1)*nnodes]==0
	 ){
          change = change - 1.0;
	 }
       }
      }
    }
    
    *(mtp->dstats) += edgeflag ? -change : change;
//  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change);
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_transitivedynamic (int ntoggles, Vertex *heads, Vertex *tails, 
	                  ModelTerm *mtp, Network *nwp)  {
  Edge e;
  long int nnodes;
  Vertex h, t, node3;
  double change;
  int i, edgeflag;
  
  nnodes = ((long int)(mtp->inputparams[0]));
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (node3 != h){
       if (EdgetreeSearch(h, node3, nwp->outedges) == 0){
       /* So form a h -> t -> node3 intransitive triad */
       /* Check to see if a prior intransitive triad existed */
         if( !(mtp->inputparams[1+(h-1)+(t-1)*nnodes]==1 &&
          mtp->inputparams[1+(t-1)+(node3-1)*nnodes]==1 &&
          mtp->inputparams[1+(h-1)+(node3-1)*nnodes]==0)
	 ){
//           Rprintf("h %d t %d node3 %d nnodes %d\n",h,t, node3, nnodes);
          change = change - 1.0;
	 }
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, t);
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
      if (node3 != h){
       if (EdgetreeSearch(h, node3, nwp->outedges) != 0){
       /* So dissolve a h -> node3 -> t intransitive triad */
       /*Check to see if a prior transitive triad existed */
         if( !(mtp->inputparams[1+(h-1)+(t-1)*nnodes]==0 &&
          mtp->inputparams[1+(h-1)+(node3-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(t-1)*nnodes]==1)
	 ){
          change = change + 1.0;
	 }
       }
      }
    }
    for(e = EdgetreeMinimum(nwp->inedges, h);
	    (node3 = nwp->inedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of hail */
    {
      if (node3 != t){
       if (EdgetreeSearch(node3, t, nwp->outedges) == 0){
       /* So form a node3 -> h -> t intransitive triad */
       /*Check to see if a prior transitive triad existed */
         if( !(mtp->inputparams[1+(h-1)+(t-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(h-1)*nnodes]==1 &&
          mtp->inputparams[1+(node3-1)+(t-1)*nnodes]==0)
	 ){
          change = change - 1.0;
	 }
       }
      }
    }
    
    *(mtp->dstats) += edgeflag ? -change : change;
//  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change);
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_heideriandynamic (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  long int nnodes;
  Vertex h, t;
  int i, edgeflag, edgeflagth;
  int edgeflagp, edgeflagthp;
  
  nnodes = ((long int)(mtp->inputparams[0]));
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      edgeflagth = (EdgetreeSearch(t, h, nwp->outedges) != 0);
      /*Get the prior edge state*/
      edgeflagp  =(int)(mtp->inputparams[1+(h-1)+(t-1)*nnodes]);
      edgeflagthp=(int)(mtp->inputparams[1+(t-1)+(h-1)*nnodes]);

//      if(edgeflagp != edgeflagthp){
//        if(edgeflagth!=edgeflag){ *(mtp->dstats) += 1.0; }
//      }else{
//        if(edgeflagth==edgeflag){ *(mtp->dstats) -= 1.0; }

      if(edgeflagp != edgeflagthp){
        if(edgeflag != edgeflagth){ 
	 *(mtp->dstats) += 1.0;
        }else{
         *(mtp->dstats) -= 1.0;
	}
      }

//      Rprintf("h %d t %d edgeflagp %d  edgeflagthp %d edgeflag %d edgeflagth %d dstats %f\n",h,t,edgeflagp,edgeflagthp, edgeflag, edgeflagth,  *(mtp->dstats) );
      if (i+1 < ntoggles)
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 

//  *(mtp->dstats) = 0.0;
//  Rprintf("nedges %d\n", nwp->nedges );
//  Rprintf("dstats %f\n", *(mtp->dstats) );
}

void d_factor (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  double val, s;
  Vertex h, t;
  long int nlevels;
  int i, j;
  
  nlevels = (long int)(mtp->inputparams[0]);

  for (i=0; i < mtp->nstats; i++){
    mtp->dstats[i] = 0.0;
  }

  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(mtp->nstats); j++) 
    {
      /*Get the covariate value*/
      val = mtp->inputparams[1+j+(h-1)*nlevels];
//          + mtp->inputparams[1+j+(t-1)*nlevels];
      mtp->dstats[j] += s*val;
    }
    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_dsp
*****************/
void d_dsp (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp) 
{
  Edge e, f;
  int i, j, echange;
  int L2hu, L2ut;
  Vertex deg;
  Vertex h, t, u, v;
  
  for (i=0; i < mtp->nstats; i++) 
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
	for(j = 0; j < mtp->nstats; j++){
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
	for(j = 0; j < mtp->nstats; j++){
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
	for(j = 0; j < mtp->nstats; j++){
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
	for(j = 0; j < mtp->nstats; j++){
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += ((L2ut + echange == deg)
			     - (L2ut == deg));
	}
      }
    }
    
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
	   
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_edgecov
*****************/
void d_edgecov (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  double val;
  Vertex h, t;
  long int nactors, nrows;
  int i, edgeflag;
  
  nrows = (long int)(mtp->inputparams[0]);
  nactors = nwp->bipartite;

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial edge state*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Get the covariate value*/
      val = mtp->attrib[(t-1-nactors)+(h-1)*nrows];
//  Rprintf("h %d t %d nactors %d val %f\n",h, t, nactors, val);
      /*Update the change statistic, based on the toggle type*/
      *(mtp->dstats) += edgeflag ? -val : val;
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}
/*****************
 void d_esp
*****************/
void d_esp (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp) 
{
  Edge e, f;
  int i, j, echange;
  int L2ht, L2hu, L2ut;
  Vertex deg;
  Vertex h, t, u, v;
  
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++){      
    L2ht=0;
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
    /* step through outedges of t */
    for(e = EdgetreeMinimum(nwp->outedges, t);
	(u = nwp->outedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->outedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=0;
	L2ut=0;
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
	for(j = 0; j < mtp->nstats; j++){
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += ((L2hu + echange == deg)
			     - (L2hu == deg));
	  mtp->dstats[j] += ((L2ut + echange == deg)
			     - (L2ut == deg));
	}
      }
    }
    /* step through inedges of t */
    for(e = EdgetreeMinimum(nwp->inedges, t);
	(u = nwp->inedges[e].value) != 0;
	e = EdgetreeSuccessor(nwp->inedges, e)){
      if (EdgetreeSearch(MIN(u,h), MAX(u,h), nwp->outedges) != 0){
	L2ht++;
	L2hu=0;
	L2ut=0;
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
	for(j = 0; j < mtp->nstats; j++){
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += ((L2hu + echange == deg)
			     - (L2hu == deg));
	  mtp->dstats[j] += ((L2ut + echange == deg)
			     - (L2ut == deg));
	}
      }
    }
    for(j = 0; j < mtp->nstats; j++){
      deg = (Vertex)mtp->inputparams[j];
      mtp->dstats[j] += echange*((L2ht == deg) - (0 == deg));
    }
    
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_idegree
*****************/
void d_idegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int echange, i, j;
  Edge e;
  Vertex h, t, node3, deg, td=0;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
  
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
	  
	  for(j=0; j < mtp->nstats; j++) 
	    {
	      deg = (Vertex)mtp->inputparams[j];
	      mtp->dstats[j] += (td + echange == deg) - (td == deg);
	    }
	}
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	td = nwp->indegree[t];
	
	for(j=0; j < mtp->nstats; j++) 
	  {
	    deg = (Vertex)mtp->inputparams[j];
	    mtp->dstats[j] += (td + echange == deg) - (td == deg);
	  }
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_istar
*****************/
void d_istar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  double change, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++)
      {
	/* edgeflag is 1 if edge exists and will disappear
           edgeflag is 0 if edge DNE and will appear */
	edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	hattr = mtp->attrib[h-1];
	if(hattr == mtp->attrib[t-1]){
	  td = - edgeflag;
	  for(e = EdgetreeMinimum(nwp->inedges, t);
	      (node3 = nwp->inedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
	    {
	      if(hattr == mtp->attrib[node3-1]){++td;}
	    }
	  
	  for(j=0; j < mtp->nstats; j++) 
	    {
	      kmo = ((int)mtp->inputparams[j]) - 1;
	      change = CHOOSE(td, kmo); 
	      mtp->dstats[j] += (edgeflag ? - change : change); 
	    }
	}
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	/* edgeflag is 1 if edge exists and will disappear
	   edgeflag is 0 if edge DNE and will appear */
	edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	td = nwp->indegree[t] - edgeflag;
	
	for(j=0; j < mtp->nstats; j++) 
	  {
	    kmo = ((int)mtp->inputparams[j]) - 1;
	    change = CHOOSE(td, kmo); 
	    mtp->dstats[j] += (edgeflag ? - change : change); 
	  }
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_kstar
*****************/
void d_kstar (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  double change, hd, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
    
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++)
    {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hattr = mtp->attrib[h-1];
      if(hattr == mtp->attrib[t-1]){
        hd = - edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, h);
        (node3 = nwp->outedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
        {
          if(hattr == mtp->attrib[node3-1]){++hd;}
        }
        for(e = EdgetreeMinimum(nwp->inedges, h);
        (node3 = nwp->inedges[e].value) != 0;
        e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
        {
          if(hattr == mtp->attrib[node3-1]){++hd;}
        }
        td = - edgeflag;
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
        
        for(j=0; j < mtp->nstats; j++) 
        {
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }else{
    for (i=0; i < ntoggles; i++)
    {
      /* edgeflag is 1 if edge exists and will disappear
      edgeflag is 0 if edge DNE and will appear */
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      hd = nwp->outdegree[h] + nwp->indegree[h] - edgeflag; 
      td = nwp->outdegree[t] + nwp->indegree[t] - edgeflag;
      for(j=0; j < mtp->nstats; j++) 
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}


/*****************
 void d_nodematch
*****************/
void d_nodematch (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  double matchval, checksum=0.0;
  Vertex h, t, ninputs;
  int i, j, edgeflag=0, matchflag;
  
  ninputs = mtp->ninputparams - nwp->nnodes;
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      h=heads[i];
      t=tails[i];
      matchflag = ((matchval=mtp->attrib[h-1]) == mtp->attrib[t-1]);
      if (matchflag) 
	edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
      if (ninputs==0) /* diff=F in network statistic specification */
	{
	  *(mtp->dstats) += matchflag ? (edgeflag ? -1.0 : 1.0) : 0.0;
	}
      else /* diff=T (and more than one category?)  */
	{
	  for (checksum=0.0, j=0; j<ninputs; j++) 
	    {
	      mtp->dstats[j] += (matchflag && matchval==mtp->inputparams[j]) ? 
		(edgeflag ? -1.0 : 1.0) : 0.0;
	      checksum += mtp->dstats[j];
	    }
	}
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_mutual

 (1,1) -> anything = -1
 anything -> (1,1) = +1
*****************/
void d_mutual (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp) {
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_asymmetric
*****************/
void d_asymmetric (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int i, edgeflag, refedgeflag;
  
  *(mtp->dstats) = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      refedgeflag = (EdgetreeSearch(t, h, nwp->outedges) != 0);
      
      *(mtp->dstats) += ((edgeflag==refedgeflag) ? 1.0 : -1.0); 
      
      if (i+1 < ntoggles) 
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_nodemain
*****************/
void d_nodemain (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp) {
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_sendercov (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp) {
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_receivercov (int ntoggles, Vertex *heads, Vertex *tails, 
		    ModelTerm *mtp, Network *nwp) {
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_hamming (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int i, nhedge, discord;
  
  nhedge = nwp[1].nedges;
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    /*Get the initial state of the edge and its alter in x0*/
//  edgeflag =(EdgetreeSearch(h=heads[i], t=tails[i], nwp[0].outedges) != 0);
    discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
//    if(!nwp[0].directed_flag && h < t){
//      hh = t;
//      ht = h;
//    }else{
//      hh = h;
//      ht = t;
//    }
    // if we will dissolve an edge discord=-1
    // discord = edgeflag ? -1 : 1;
    

//  so moving away one step
//    discord = (edgeflag0!=edgeflag) ? -1 : 1;

//  Rprintf("h %d t %d edgeflag %d edgeflag0 %d discord %d\n",h, t, edgeflag, edgeflag0, discord);
//  if(nhedge>0)
//  Rprintf("h %d t %d discord %d nhedge %d\n",h, t, discord, nhedge);

    /*Update the change statistics, as appropriate*/
//    *(mtp->dstats) += ((edgeflag0!=edgeflag) ? -1.0 : 1.0);

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

void d_hammingmix (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int i, j, nhedge, edgeflag, discord;
  int matchvalh, matchvalt;
  int nstats;
  
  nhedge =  mtp->inputparams[0];
  nstats = mtp->nstats;
//  Rprintf("nstats %d nhedge %d i0 %f i1 %f i2 %f i3 %f\n",nstats, nhedge, mtp->inputparams[0],
//                                 mtp->inputparams[1],
//                                 mtp->inputparams[2],
//                                 mtp->inputparams[3]
//		  );
  for (i=0; i < mtp->nstats; i++)
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
//   Rprintf("h %d t %d matchvalh %d matchvalt %d edgeflag %d discord %d j %d p0 %f p1 %f\n",h,t,matchvalh,matchvalt,edgeflag,discord,j,mtp->inputparams[2*nhedge+  j], mtp->inputparams[2*nhedge+ nstats+j]);
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

void d_hammingfixmix (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int i, nhedge, discord;
  int matchvalh, matchvalt;
  
  nhedge = mtp->inputparams[0];
//  Rprintf("nhedge %d\n", nhedge);
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
//  Rprintf("Warning: hammingfixmix can only be used with ConstantEdges terms.\n");
//  Rprintf("nhedge %d i0 %f i1 %f i2 %f i3 %f\n", nhedge, mtp->inputparams[0],
//                                 mtp->inputparams[1],
//                                 mtp->inputparams[2],
//                                 mtp->inputparams[3]
//		  );
     
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

void d_duration (int ntoggles, Vertex *heads, Vertex *tails, 
                 ModelTerm *mtp, Network *nwp) {
  Vertex h, t, hh, ht;
  int i, k, nhedge, edgeflag, discord, lookmore;
  int ndyads, nnodes;
  double change=0.0;
  
  nhedge = (int)((mtp->inputparams[0]));
  nnodes = nwp->nnodes;
  ndyads = nwp->nnodes * nwp->nnodes;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    /*Get the initial state of the edge and its reflection*/
    edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    if(!nwp->directed_flag && h < t){
      hh = t;
      ht = h;
    }else{
      hh = h;
      ht = t;
    }
    discord = edgeflag ? 1 : -1;
    
    k=0;
    lookmore=1;
    while(lookmore && k < nhedge){
      if(hh == (Vertex)(mtp->attrib[         k]) &&
        ht == (Vertex)(mtp->attrib[nhedge + k])
      ){
        /*Get dissolution rate */
        lookmore=2*nhedge+ndyads+(hh-1)+(ht-1)*nnodes;
        change=discord*(double)(mtp->attrib[lookmore]);
        lookmore=0;
      }else{
        ++k;
      }
    }
    if(lookmore){
      /* Get formation rate */
      lookmore=2*nhedge+(hh-1)+(ht-1)*nnodes;
      change=-discord*(double)(mtp->attrib[lookmore]);
    }
    /*Update the change statistics, as appropriate*/
    *(mtp->dstats) += change;
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_nodefactor
*****************/
void d_nodefactor (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  double s, factorval;
  Vertex h, t, nlevels;
  int i, j;
  
  nlevels = (mtp->ninputparams) - nwp->nnodes;
  for (i=0; i < mtp->nstats; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(mtp->nstats); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[h-1] != factorval) ? 0.0 : s)
	    + ((mtp->attrib[t-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_nodeifactor
*****************/
void d_nodeifactor (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  double s, factorval;
  Vertex h, t, nlevels;
  int i, j;
  
  nlevels = (mtp->ninputparams) - nwp->nnodes;
  for (i=0; i < mtp->nstats; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(mtp->nstats); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[t-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_nodeofactor
*****************/
void d_nodeofactor (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  double s, factorval;
  Vertex h, t, nlevels;
  int i, j;
  
  nlevels = (mtp->ninputparams) - nwp->nnodes;
  for (i=0; i < mtp->nstats; i++){
    mtp->dstats[i] = 0.0;
  }
  for (i=0; i<ntoggles; i++) 
  {
    s = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0) ? -1.0 : 1.0;
    for (j=0; j<(mtp->nstats); j++) 
    {
      factorval = (mtp->inputparams[j]);
      mtp->dstats[j] += 
	    ((mtp->attrib[h-1] != factorval) ? 0.0 : s);
      
    }
    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_odegree
*****************/
void d_odegree (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int echange, i, j;
  Edge e;
  Vertex h, t, node3, deg, td=0;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
  
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
	  
	  for(j=0; j < mtp->nstats; j++) 
	    {
	      deg = (Vertex)mtp->inputparams[j];
	      mtp->dstats[j] += (td + echange == deg) - (td == deg);
	    }
	}
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
	td = nwp->outdegree[h];
	
	for(j=0; j < mtp->nstats; j++) 
	  {
	    deg = (Vertex)mtp->inputparams[j];
	    mtp->dstats[j] += (td + echange == deg) - (td == deg);
	  }
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_ostar
*****************/
void d_ostar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  double change, td=0.0;
  int edgeflag, i, j, kmo;
  Edge e;
  Vertex h, t, node3;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
  
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  if(ninputs>nstats){
    /* match on attributes */
    for (i=0; i < ntoggles; i++)
      {
	/* edgeflag is 1 if edge exists and will disappear
           edgeflag is 0 if edge DNE and will appear */
	edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	hattr = mtp->attrib[t-1];
	if(hattr == mtp->attrib[h-1]){
	  td = - edgeflag;
	  for(e = EdgetreeMinimum(nwp->outedges, h);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
	    {
	      if(hattr == mtp->attrib[node3-1]){++td;}
	    }
	  
	  for(j=0; j < mtp->nstats; j++) 
	    {
	      kmo = ((int)mtp->inputparams[j]) - 1;
	      change = CHOOSE(td, kmo); 
	      mtp->dstats[j] += (edgeflag ? - change : change); 
	    }
	}
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }else{
    for (i=0; i < ntoggles; i++)
      {
	/* edgeflag is 1 if edge exists and will disappear
	   edgeflag is 0 if edge DNE and will appear */
	edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
	td = nwp->outdegree[h] - edgeflag;
	
	for(j=0; j < mtp->nstats; j++) 
	  {
	    kmo = ((int)mtp->inputparams[j]) - 1;
	    change = CHOOSE(td, kmo); 
	    mtp->dstats[j] += (edgeflag ? - change : change); 
	  }
	if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_smalldiff
*****************/
void d_smalldiff (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) {
    h=heads[i];
    t=tails[i];
    *(mtp->dstats) += (abs(mtp->attrib[h-1] - mtp->attrib[t-1])
    > *(mtp->inputparams)) ? 0.0 :
    ((EdgetreeSearch(h, t, nwp->outedges) != 0) ? -1.0 : 1.0); 
    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_ctriad
*****************/
void d_ctriad (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr;
  
  ninputs = mtp->ninputparams;
  nstats = mtp->nstats;
  
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_ttriad
*****************/
void d_ttriad (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i, j;
  int ninputs, nstats;
  double hattr;
  
  ninputs = mtp->ninputparams;
  nstats = mtp->nstats;
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_tripercent
*****************/
void d_tripercent (int ntoggles, Vertex *heads, Vertex *tails, 
ModelTerm *mtp, Network *nwp) {
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
  
  nstats = mtp->nstats;
  ninputs = mtp->ninputparams;
  nnodes = nwp->nnodes;
  
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
        if(abs(mtp->attrib[h-1] - mtp->attrib[t-1])<eps){
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
        if(abs(hattr - mtp->attrib[node3-1])<eps){++hd;}
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
        if(abs(hattr - mtp->attrib[node3-1])<eps){
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
        if(abs(hattr - mtp->attrib[node3-1])<eps){
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
      if(abs(hattr - mtp->attrib[t-1])<eps){
        hd = -edgeflag;
        for(e = EdgetreeMinimum(nwp->outedges, h);
	      (node3 = nwp->outedges[e].value) != 0;
	      e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of head */
        {
          if(abs(hattr - mtp->attrib[node3-1])<eps){++hd;}
        }
        if (!nwp->directed_flag){
          for(e = EdgetreeMinimum(nwp->inedges, h);
          (node3 = nwp->inedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of head */
          {
            if(abs(hattr - mtp->attrib[node3-1])<eps){++hd;}
          }
          td = - edgeflag;
          for(e = EdgetreeMinimum(nwp->outedges, t);
          (node3 = nwp->outedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
          {
            if(abs(hattr - mtp->attrib[node3-1])<eps){++td;}
          }
          for(e = EdgetreeMinimum(nwp->inedges, t);
          (node3 = nwp->inedges[e].value) != 0;
          e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
          {
            if(abs(hattr - mtp->attrib[node3-1])<eps){++td;}
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
          if(abs(hattr - mtp->attrib[node3-1])<eps){
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
          if(abs(hattr - mtp->attrib[node3-1])<eps){
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void CountTriangle
*****************/
Vertex CountTriangles (Vertex h, Vertex t, int outcount, int incount, 
		       Network *nwp) {
  Edge e;
  Vertex change=0;
  Vertex k;
  
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

/*****************
 void d_boundeddegree
*****************/
void d_boundeddegree (int ntoggles, Vertex *heads, Vertex *tails, 
		      ModelTerm *mtp, Network *nwp) 
{
  int i, j, echange;
  Vertex h, t, hd, td=0, deg, *id, *od;
  TreeNode *oe;  
  int nstats = (int)mtp->nstats;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  oe=nwp->outedges;
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedodegree
*****************/
void d_boundedodegree (int ntoggles, Vertex *heads, Vertex *tails, 
	               ModelTerm *mtp, Network *nwp) {
  int i, j, echange;
  Vertex h, t, hd=0, deg;
  int nstats = (int)mtp->nstats;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = nwp->outdegree[h];
      for(j = 0; j < mtp->nstats; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
	}
      mtp->dstats[nstats-1] += (hd + echange >= bound) - (hd >= bound);
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_boundedidegree
*****************/
void d_boundedidegree (int ntoggles, Vertex *heads, Vertex *tails, 
	               ModelTerm *mtp, Network *nwp) {
  int i, j, echange;
  Vertex h, t, hd=0, deg;
  int nstats = (int)mtp->nstats;
  Vertex bound = (Vertex)mtp->inputparams[nstats-1];
  
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = nwp->indegree[h];
      for(j = 0; j < mtp->nstats; j++) 
	{
	  deg = (Vertex)mtp->inputparams[j];
	  mtp->dstats[j] += (hd + echange == deg) - (hd == deg);
	}
      mtp->dstats[nstats-1] += (hd + echange >= bound) - (hd >= bound);
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_sender (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp) {
  int i, j, echange;
  Vertex h, t, deg;
  
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      if(h == 1){
       echange = -echange;
       for (j=0; j < mtp->nstats; j++){
         deg = (Vertex)mtp->inputparams[j];
         if(deg != 1){mtp->dstats[j] += echange;}
       }
      }else{
       j=0;
       deg = (Vertex)mtp->inputparams[j];
       while(deg != h && j < mtp->nstats){
	j++;
	deg = (Vertex)mtp->inputparams[j];
       }
       if(j < mtp->nstats){mtp->dstats[j] += echange;}
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_receiver (int ntoggles, Vertex *heads, Vertex *tails, 
	         ModelTerm *mtp, Network *nwp) {
  int i, j, echange;
  Vertex h, t, deg;
  
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  
  for (i=0; i<ntoggles; i++)
    {      
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      if(t == 1){
       echange = -echange;
       for (j=0; j < mtp->nstats; j++){ 
         deg = (Vertex)mtp->inputparams[j];
         if(deg != 1){mtp->dstats[j] += echange;}
       }
      }else{
       j=0;
       deg = (Vertex)mtp->inputparams[j];
       while(deg != t && j < mtp->nstats){
	j++;
	deg = (Vertex)mtp->inputparams[j];
       }
       if(j < mtp->nstats){mtp->dstats[j] += echange;}
      }
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_sociality (int ntoggles, Vertex *heads, Vertex *tails, 
	           ModelTerm *mtp, Network *nwp) {
  int i, j, echange;
  Vertex h, t, deg;
  int ninputs, nstats;
  double hattr;
  
  ninputs = (int)mtp->ninputparams;
  nstats  = (int)mtp->nstats;
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
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
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
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
      }
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_localtriangle (int ntoggles, Vertex *heads, Vertex *tails, 
		    ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, change, node3;
  int edgeflag, i;
  int ninputs, nstats;
  long int nmat;
  
  ninputs = mtp->ninputparams;
  nstats = mtp->nstats;
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
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*
d_berninhom - Changescores for inhomogeneous Bernoulli graphs
*/
void d_berninhom (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)
{
  int edgeflag, i, nstats;
  Vertex h, t, n;

  nstats=(int)mtp->nstats;
  n=nwp->nnodes;
  for(i=0;i<nstats;i++)
    mtp->dstats[i] = 0.0;
    
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Lower trianglization of adjacency matrix implies ith element corresponds
      to (row-1) + (col-1)*(n-1) - choose(col,2).*/
      /*Rprintf("t=%ld, h=%ld, cell=%ld, nstats=%d, state=%d\n",t,h, (t-1)+(h-1)*(n-1)-h*(h-1)/2-1, nstats, edgeflag);
      Rprintf("\tdstats content=%f\n",mtp->dstats[(t-1)+(h-1)*(n-1)-h*(h-1)/2-1]);*/
      mtp->dstats[(t-1)+(h-1)*(n-1)-h*(h-1)/2-1] += edgeflag ? - 1 : 1;

      if (i+1 < ntoggles)
        ToggleEdge(h, t, nwp);  /* Toggle this edge if more to come */
    }
  i--;
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp);
}


/*
d_spatial - Changescores for inhomogeneous Bernoulli graphs
*/
void d_spatial (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)
{
  int edgeflag, i, nstats;
  Vertex h, t, n;
  double llr,pb,alpha,gamma;

  nstats=(int)mtp->nstats;
  n=nwp->nnodes;
  pb=mtp->inputparams[0];
  alpha=mtp->inputparams[1];
  gamma=mtp->inputparams[2];
  for(i=0;i<nstats;i++)
    mtp->dstats[i] = 0.0;
    
  for (i=0; i < ntoggles; i++)
    {
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      /*Lower trianglization of adjacency matrix implies ith element corresponds
      to (row-1) + (col-1)*(n-1) - choose(col,2).*/
      llr = -log(((1+exp(pb))*pow(1+exp(alpha)*(mtp->inputparams[2+(t-1)+(h-1)*(n-1)-h*(h-1)/2]),exp(gamma)))/exp(pb)-1);
      *(mtp->dstats) += edgeflag ? - llr : llr;

      if (i+1 < ntoggles)
        ToggleEdge(h, t, nwp);  /* Toggle this edge if more to come */
    }
  i--;
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp);
}


/*
Cycle stuff starts here!
*/
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

void d_cycle (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)
{
  int edgeflag,i,j,k,nstats,directed;
  Vertex h, t;
  long int maxlen;
  double *countv,emult;
  
  /*Perform initial setup*/
  directed=(int)(mtp->inputparams[0]);
  maxlen=(long int)(mtp->inputparams[1]);
  nstats=(int)mtp->nstats;
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
    ToggleEdge(heads[i], tails[i], nwp);
}


/*****************
 void d_nodemix

*****************/
void d_nodemix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp) {
  Vertex h, t, ninputs, ninputs2;
  int i, j, edgeflag=0, matchflag;
  double rtype, ctype;

  ninputs = mtp->ninputparams - nwp->nnodes;
  ninputs2 = ninputs/2;
  for (i=0; i < mtp->nstats; i++)
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++)
    {
      h=heads[i];
      t=tails[i];
      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0); /*Get edge state*/
      matchflag=0;
      /*Find the node covariate values (types) for the head and tail*/
      rtype=MIN(mtp->inputparams[h+ninputs-1],mtp->inputparams[t+ninputs-1]);
      ctype=MAX(mtp->inputparams[h+ninputs-1],mtp->inputparams[t+ninputs-1]);
      /*Find the right statistic to update*/
      for(j=0;(j<ninputs2)&&(!matchflag);j++){
        if((mtp->inputparams[j          ]==rtype)&&
           (mtp->inputparams[j+ninputs2]==ctype)){
            mtp->dstats[j] += (edgeflag ? -1.0 : 1.0);
            matchflag++;
        }
      }
      if (i+1 < ntoggles)
        ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--;
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp);
}

/*****************
 void d_mix
*****************/
void d_mix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  int matchvalh, matchvalt;
  int i, j, edgeflag=0, nstats;

  nstats = mtp->nstats;
  for (i=0; i < mtp->nstats; i++)
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++)
    {
      h=heads[i];
      t=tails[i];
      matchvalh = mtp->inputparams[h-1+2*nstats];
      matchvalt = mtp->inputparams[t-1+2*nstats];
      edgeflag=(EdgetreeSearch(h, t, nwp[0].outedges) != 0); /*Get edge state*/
      for (j=0; j<nstats; j++) 
	  {
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
 void d_simmelian
*****************/
void d_simmelian (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
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
     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_simmelianties (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
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
//   if(firstht){++change;}
      
   change = 2*change;
   *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
   }
   
   if (i+1 < ntoggles) 
     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

///*****************
// void d_nearsimmelian
//*****************/
//void d_nearsimmelian (int ntoggles, Vertex *heads, Vertex *tails, 
//                ModelTerm *mtp, Network *nwp) {
//  Edge e;
//  long int n;
//  Vertex h, t, change, node3;
//  int edgeflag, i;
//
//  n=g->nnodes;
//  
// *(mtp->dstats) = 0.0;
// for (i=0; i<ntoggles; i++) 
// {
//  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
//   
//  sc = edgeflag;
//  change = 0;
//   
//  for(i=1;i<=n;i++)
//    if((i!=h)&&(i!=t)){
//
//   for(e = EdgetreeMinimum(nwp->outedges, t);
//       (node3 = nwp->outedges[e].value) != 0;
//       e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
//   {
//     sc = edgeflag + 1;
//     sc += EdgetreeSearch(t, h, nwp->outedges) != 0;
//     sc += EdgetreeSearch(node3, h, nwp->outedges) != 0 ;
//     sc += EdgetreeSearch(h, node3, nwp->outedges) != 0 ;
//     sc += EdgetreeSearch(node3, t, nwp->outedges) != 0 ;
//
//     if (sc == 6 ){--change} /* so edgeflag==1 and we distroy a simmelian */
//     if (sc == 5 && edgeflag == 0 ){--change}
//     if (sc == 5 && edgeflag == 1 ){++change}
//     if (sc == 4 && edgeflag == 0 ){++change}
//   }
//      
//   *(mtp->dstats) += edgeflag ? -(double)change : (double)change;
//   }
//   
//   if (i+1 < ntoggles) 
//     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//  }
//  
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}

/*****************
 void d_nearsimmelian
*****************/
void d_nearsimmelian (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
  Vertex h, t, node3;
  double change;
  int edgeflag, i, edgeflagth, sc;

 *(mtp->dstats) = 0.0;

 for (i=0; i<ntoggles; i++) 
 {
  edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
  edgeflagth = (EdgetreeSearch(t, h, nwp->outedges) == 0);
   
  for(node3=1;node3<=nwp->nnodes;node3++){
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
     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_hiertriad
*****************/
void d_hiertriad (int ntoggles, Vertex *heads, Vertex *tails, 
                  ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, node3;
  double pos3;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    
//           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag);
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (EdgetreeSearch(h, node3, nwp->outedges) != 0){
        pos3 = numposthree (node3, nwp);
        *(mtp->dstats) -= pos3;
        ToggleEdge(h, t, nwp);
        pos3 = numposthree (node3, nwp);
        *(mtp->dstats) += pos3;
        ToggleEdge(h, t, nwp);
      }
    }
    
    pos3 = numposthree (t, nwp);
    *(mtp->dstats) -= pos3;
    ToggleEdge(h, t, nwp);
    pos3 = numposthree (t, nwp);
    *(mtp->dstats) += pos3;
    ToggleEdge(h, t, nwp);
// Rprintf("h %d t %d pos3 %f ideg %d\n",h,t, pos3, ideg);

    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 void d_hiertriaddegree
*****************/
void d_hiertriaddegree (int ntoggles, Vertex *heads, Vertex *tails, 
                  ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, node3, ideg;
  double pos3;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    
//           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag);
    for(e = EdgetreeMinimum(nwp->outedges, t);
	    (node3 = nwp->outedges[e].value) != 0;
	    e = EdgetreeSuccessor(nwp->outedges, e)) /* step through outedges of tail */
    {
      if (EdgetreeSearch(h, node3, nwp->outedges) != 0){
        ideg = nwp->indegree[node3];
        pos3 = numposthree (node3, nwp);
        *(mtp->dstats) -= pos3*ideg;
        ToggleEdge(h, t, nwp);
        ideg = nwp->indegree[node3];
        pos3 = numposthree (node3, nwp);
        *(mtp->dstats) += pos3*ideg;
        ToggleEdge(h, t, nwp);
      }
    }
    
    ideg = nwp->indegree[t];
    pos3 = numposthree (t, nwp);
    *(mtp->dstats) -= pos3*ideg;
    ToggleEdge(h, t, nwp);
    ideg = nwp->indegree[t];
    pos3 = numposthree (t, nwp);
    *(mtp->dstats) += pos3*ideg;
    ToggleEdge(h, t, nwp);
// Rprintf("h %d t %d pos3 %f ideg %d\n",h,t, pos3, ideg);

    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

double numposthree (Vertex t, Network *nwp) {
  Edge e, f;
  Vertex pos, node2, node3;
  double dpos;
  
  pos = 0;
  for(e = EdgetreeMinimum(nwp->inedges, t); 
   (node2 = nwp->inedges[e].value) != 0;
    e = EdgetreeSuccessor(nwp->inedges, e)) /* step through inedges of tail */
    {
     for(f = EdgetreeMinimum(nwp->inedges, node2); 
	    (node3 = nwp->inedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->inedges, f)) /* step through inedges of node2 */
      { if (EdgetreeSearch(node3, t, nwp->outedges) != 0){++pos;} }
     for(f = EdgetreeMinimum(nwp->outedges, node2); 
	    (node3 = nwp->outedges[f].value) != 0;
	    f = EdgetreeSuccessor(nwp->outedges, f)) /* step through inedges of node2 */
      { if (EdgetreeSearch(node3, t, nwp->outedges) != 0){++pos;} }
    }
  dpos = pos / 2.0;
  return dpos;
}

void d_intransitive (int ntoggles, Vertex *heads, Vertex *tails, 
                  ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, node2;
  double change;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
//           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag);
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
//  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change);

    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_transitive (int ntoggles, Vertex *heads, Vertex *tails, 
                  ModelTerm *mtp, Network *nwp) {
  Edge e;
  Vertex h, t, node2;
  double change;
  int edgeflag, i;
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
  {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    change = 0.0;
    
//           Rprintf("h %d t %d edgeflag %d\n",h,t, edgeflag);
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
//  Rprintf("h %d t %d edgeflag %d change %f\n",h,t, edgeflag, change);

    if (i+1 < ntoggles) 
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

void d_hammingdyadcov (int ntoggles, Vertex *heads, Vertex *tails, 
                ModelTerm *mtp, Network *nwp) {
  Vertex h, t;
  double val;
  long int nnodes, nactors, nevents, n0edge;
  int i, j, discord;
  
  n0edge =  mtp->inputparams[0];
  nnodes = nwp[0].nnodes;
  nactors = nwp[0].bipartite;
  nevents = nwp[0].nnodes - nactors;
//  Rprintf("nactors %d i0 %f i1 %f i2 %f i3 %f\n", nactors,
//                                 mtp->inputparams[0],
//                                 mtp->inputparams[1],
//                                 mtp->inputparams[2],
//                                 mtp->inputparams[3]
//		  );
//  for (i=0; i<1000; i++) {
//  Rprintf("i %d inp %f\n", i, mtp->inputparams[i]);
//  }

  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial discord state*/
//    edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[0].outedges) != 0);
      discord=(EdgetreeSearch(h=heads[i], t=tails[i], nwp[1].outedges) != 0);
      /*Get the covariate value*/
      val = mtp->inputparams[1+(t-nactors-1)*nactors+(h-1)+2*n0edge];
      /*Update the change statistic, based on the toggle type*/
      *(mtp->dstats) += discord ? -val : val;
 // Rprintf("nnodes %d n0edge %d h %d t %d discord %d val %f\n",nnodes, n0edge, h, t-nactors, discord, val);
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
