#include "changestats.ihs.h"


/********************  changestats:  A    ***********/
/*****************
 changestat: d_akappa
*****************/
void d_akappa (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp)  {
  int i, j, echange=0;
  double nedges, change, iar0, far0;
  Vertex h, t, hd, iak2, fak2, nnodes, *od;
  Vertex nactors;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  od=nwp->outdegree;
  nnodes = nwp->nnodes;
  nactors = nwp->bipartite;
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    iak2=0;
    for (j=1; j<=nactors; j++) {      
      fak2 = od[j];
      iak2 += fak2*(fak2-1);
    }
    hd = od[h] + (echange-1)/2;
    fak2 = iak2 + echange*2*hd;
    nedges = (double)(nwp->nedges);
    iar0 = (nwp->nedges==0) ? 0.0 : (iak2*1.0/nedges);
    far0 = (((nwp->nedges)+echange)==0) ? 0.0 : (fak2*1.0/(nedges+echange));
    change += far0 - iar0;
//   Rprintf("h %d t %d nnodes %d nedges %f iak2 %d fak2 %d iar0 %f far0 %f change %f\n",h,t, nnodes,  nedges, iak2, fak2, iar0, far0, change);
//   Rprintf("h %d t %d nnodes %d nedges %f iek2 %d fek2 %d ier0 %f fer0 %f change %f\n",h,t, nnodes,  nedges, iek2, fek2, ier0, fer0, change);
      
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }                                         
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_altistar
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
 changestat: d_altostar
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

/********************  changestats:  B    ***********/
/*****************
 changestat: d_berninhom
 Changescores for inhomogeneous Bernoulli graphs
*****************/
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

/*****************
 changestat: d_biduration
*****************/
void d_biduration (int ntoggles, Vertex *heads, Vertex *tails,
		   ModelTerm *mtp, Network *nwp)
{
  Vertex h, t, hh, ht;
  int i, k, nprevedge, edgeflag, discord, lookmore;
  Vertex nevents, nactors;
  double change=0.0;

  nprevedge = (int)((mtp->inputparams[0]));
  nactors = (int)((mtp->inputparams[1]));
//  nactors = nwp->bipartite;
//  nevents = nwp->nnodes - nwp->bipartite;
  nevents = nwp->nnodes - nactors;

//  Rprintf("nprevedge %d\n", nprevedge);
//  Rprintf("nwp->bipartite %d nevents %f\n", nwp->bipartite, nwp->bipartite);
//  Rprintf("nactors %d nevents %d\n", nactors, nevents);
//  for (k=1; k<=nprevedge; k++) {
//      Rprintf("k %d x0.h %d x0.t %d\n", k,
//	      (Vertex)(mtp->attrib[          k]),
//	      (Vertex)(mtp->attrib[nprevedge+k]));
//  }

  *(mtp->dstats) = 0.0;

  for (i=0; i<ntoggles; i++) 
    {
      /*Get the initial state of the edge and its reflection*/
      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      if(h > t){
        hh = t;
        ht = h;
      }else{
        hh = h;
        ht = t;
      }
      // 1 = edge removed -1 = edge added
      discord = edgeflag ? 1 : -1;

// Rprintf("nprevedge %d\n", nprevedge);

      k=1;
      lookmore=1;
      while(lookmore && k <= nprevedge){
//       Rprintf("h %d t %d hh %d ht %d\n",
//	      h,t,hh,ht);
//       Rprintf("h %d t %d hh %d ht %d\n",
//	      ht,hh,(Vertex)(mtp->attrib[            k]),
//	            (Vertex)(mtp->attrib[nprevedge + k]));
       if(ht == (Vertex)(mtp->attrib[            k]) &&
          hh == (Vertex)(mtp->attrib[nprevedge + k])
	 ){
//       Rprintf("h %d t %d hh %d ht %d\n",
//	      h,t,(Vertex)(mtp->attrib[            k]),
//	          (Vertex)(mtp->attrib[nprevedge + k]));

           /*If the proposed edge existed in x0, get dissolution rate */
	   /* lookmore is location of (hh, ht) in the dissolve matrix */
//           lookmore=2*nprevedge+(hh-1)+(ht-1-nevents)*nactors+1;
           lookmore=2*nprevedge+(hh-1)+(ht-1-nactors)*nactors+1;
	   /* change is NEGATIVE the dissolve rate if an edge exists */
//       Rprintf("hh %d ht %d lookmore %d att %f\n",hh, ht, lookmore, (double)(mtp->attrib[lookmore]));
//           change=-discord*(double)(mtp->attrib[lookmore]);
           change=-discord;
//       Rprintf("lookmore %d change %f discord %d \n",lookmore, change, discord);
	   lookmore=0;
           /*Update the change statistics, as appropriate*/
	   /*change is the number of edges dissolved*/
           *(mtp->dstats) += (double)change;
//      if(change!=0){
//       Rprintf("h %d t %d hh %d ht %d edgeflag %d beta %f\n",
//	      h,t,hh,ht, edgeflag, (double)(mtp->attrib[lookmore]));
//       Rprintf("edgeflag %d change %f discord %d\n",
//	      edgeflag, change, discord);
//      }
          }else{
		   ++k;
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
 changestat: d_bimix
*****************/
void d_bimix (int ntoggles, Vertex *heads, Vertex *tails,
		              ModelTerm *mtp, Network *nwp) {

  int matchvalh, matchvalt;
  Vertex h, t, nnodes, nstats;
  int i, j, edgeflag=0;

  nstats = mtp->nstats;
  nnodes = (mtp->ninputparams)-2*nstats;
//  Rprintf("ninput %d nstats %d\n", mtp->ninputparams, mtp->nstats);
//  Rprintf("nodes %d\n", nnodes);
// for (h=0; h<nnodes; h++){
//   Rprintf("match h %d att h %f\n", 
//	    h, mtp->attrib[h+nstats]);
// }
// for (h=0; h<nstats; h++){
//   Rprintf("match h %d in h %f in t %f\n", 
//	    h, mtp->inputparams[h], mtp->inputparams[h+nstats]);
// }
  for (i=0; i < nstats; i++) 
    mtp->dstats[i] = 0.0;

  for (i=0; i<ntoggles; i++) 
    {
      h=heads[i];
      t=tails[i];
      matchvalh = mtp->attrib[h-1+nstats];
      matchvalt = mtp->attrib[t-1+nstats];
//     if(matchvalt < matchvalh){
//       matchswap = matchvalh;
//       matchvalh = matchvalt;
//       matchvalt = matchswap;
//      }       
      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
      for (j=0; j<nstats; j++) 
	  {
           if(matchvalh==mtp->inputparams[nstats+j] &&
	      matchvalt==mtp->inputparams[       j]
	     ){mtp->dstats[j] += edgeflag ? -1.0 : 1.0;}
	  }
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_bkappa
*****************/
void d_bkappa (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp)  {
  int i, j, echange=0;
  double nedges, change, iar0, far0, ier0, fer0;
  Vertex h, t, hd, td=0, iak2, fak2, iek2, fek2, nnodes, *id, *od;
  Vertex nevents, nactors;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  nnodes = nwp->nnodes;
  nactors = nwp->bipartite;
  nevents = nnodes - nactors;
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    iak2=0;
    for (j=1; j<=nactors; j++) {      
      fak2 = od[j];
      iak2 += fak2*(fak2-1);
    }
    iek2=0;
    for (j=nactors+1; j<=nnodes; j++) {      
      fek2 = id[j];
      iek2 += fek2*(fek2-1);
    }
    hd = od[h] + (echange-1)/2;
    td = id[t] + (echange-1)/2;
    fak2 = iak2 + echange*2*hd;
    fek2 = iek2 + echange*2*td;
    nedges = (double)(nwp->nedges);
    iar0 = (nwp->nedges==0) ? 0.0 : (iak2*1.0/nedges);
    far0 = (((nwp->nedges)+echange)==0) ? 0.0 : (fak2*1.0/(nedges+echange));
    ier0 = (nwp->nedges==0) ? 0.0 : (iek2*1.0/nedges);
    fer0 = (((nwp->nedges)+echange)==0) ? 0.0 : (fek2*1.0/(nedges+echange));
    change += sqrt(far0*fer0) - sqrt(iar0*ier0);
//   Rprintf("h %d t %d nnodes %d nedges %f iar0 %f far0 %f change %f\n",h,t, nnodes,  nedges, iar0, far0, change);
//   Rprintf("h %d t %d nnodes %d nedges %f ier0 %f fer0 %f change %f\n",h,t, nnodes,  nedges, ier0, fer0, change);
      
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_degreep
*****************/
void d_degreep (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  int i, j, echange;
  Vertex head, tail, headdeg, taildeg, deg, *id, *od;
  TreeNode *oe=nwp->outedges;

  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;  
  for (i=0; i<ntoggles; i++) {      
    echange=(EdgetreeSearch(head=heads[i], tail=tails[i], oe)==0)? 1:-1;
    headdeg = od[head] + id[head];
    taildeg = od[tail] + id[tail];
    for(j = 0; j < mtp->nstats; j++) {
      deg = (Vertex)mtp->inputparams[j];
      mtp->dstats[j] += ((headdeg + echange == deg) - (headdeg == deg))/(double)nwp->nnodes;
      mtp->dstats[j] += ((taildeg + echange == deg) - (taildeg == deg))/(double)nwp->nnodes;
    }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_degreep_by_attr
*****************/
void d_degreep_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	        ModelTerm *mtp, Network *nwp) 
{
  /* The inputparams are assumed to be set up as follows:
  The first 2*nstats values are in pairs:  (degree, attrvalue)
  The values following the first 2*nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr, testattr;
  Vertex head, tail, headdeg, taildeg, d, *id, *od;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  for (i=0; i < mtp->nstats; i++) 
    mtp->dstats[i] = 0.0;
  for (i=0; i<ntoggles; i++) {
    echange=(EdgetreeSearch(head=heads[i], tail=tails[i], oe)==0)? 1:-1;
    headdeg = od[head] + id[head];
    taildeg = od[tail] + id[tail];
    headattr = mtp->inputparams[2*mtp->nstats + head - 1]; 
    tailattr = mtp->inputparams[2*mtp->nstats + tail - 1]; 
    for(j = 0; j < mtp->nstats; j++) {
      d = (Vertex)mtp->inputparams[2*j];
      testattr = mtp->inputparams[2*j + 1]; 
      if (headattr == testattr)  /* we have head attr match */
        mtp->dstats[j] += ((headdeg + echange == d) - (headdeg == d))/(double)nwp->nnodes;
      if (tailattr == testattr)  /* we have tail attr match */
        mtp->dstats[j] += ((taildeg + echange == d) - (taildeg == d))/(double)nwp->nnodes;
    }
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_degreep_w_homophily
*****************/
void d_degreep_w_homophily (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) 
{
  /*  The inputparams are assumed to be set up as follows:
  The first nstats values are the values of degree
  The values following the first nstats values are the nodal attributes.
  */
  int i, j, echange, headattr, tailattr;
  Vertex head, tail, headdeg, taildeg, deg, tmp;
  TreeNode *ie=nwp->inedges, *oe=nwp->outedges;
  double *nodeattr;
  Edge e;

  nodeattr = mtp->inputparams + mtp->nstats - 1;  
  for (i=0; i < mtp->nstats; i++) 
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
      for(j = 0; j < mtp->nstats; j++) {
        deg = (Vertex)mtp->inputparams[j];
        mtp->dstats[j] += ((headdeg + echange == deg) - (headdeg == deg))/(double)nwp->nnodes;
	mtp->dstats[j] += ((taildeg + echange == deg) - (taildeg == deg))/(double)nwp->nnodes;
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
 changestat: d_dissolve
*****************/
void d_dissolve(int ntoggles, Vertex *heads, Vertex *tails, 
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
 changestat: d_duration
*****************/
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

/********************  changestats:  C    ***********/

/********************  changestats:  D    ***********/

/********************  changestats:  E    ***********/
/*****************
 changestat: d_ekappa
*****************/
void d_ekappa (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp)  {
  int i, j, echange=0;
  double nedges, change, ier0, fer0;
  Vertex h, t, td=0, iek2, fek2, nnodes, *id;
  Vertex nactors;
  TreeNode *oe;  
  
  oe=nwp->outedges;
  id=nwp->indegree;
  nnodes = nwp->nnodes;
  nactors = nwp->bipartite;
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    iek2=0;
    for (j=nactors+1; j<=nnodes; j++) {      
      fek2 = id[j];
      iek2 += fek2*(fek2-1);
    }
    td = id[t] + (echange-1)/2;
    fek2 = iek2 + echange*2*td;
    nedges = (double)(nwp->nedges);
    ier0 = (nwp->nedges==0) ? 0.0 : (iek2*1.0/nedges);
    fer0 = (((nwp->nedges)+echange)==0) ? 0.0 : (fek2*1.0/(nedges+echange));
    change += fer0 - ier0;
//   Rprintf("h %d t %d nnodes %d nedges %f iak2 %d fak2 %d iar0 %f far0 %f change %f\n",h,t, nnodes,  nedges, iak2, fak2, iar0, far0, change);
//   Rprintf("h %d t %d nnodes %d nedges %f iek2 %d fek2 %d ier0 %f fer0 %f change %f\n",h,t, nnodes,  nedges, iek2, fek2, ier0, fer0, change);
      
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/********************  changestats:  F    ***********/
/*****************
 changestat: d_factor
*****************/
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
 changestat: d_formation
*****************/
void d_formation(int ntoggles, Vertex *heads, Vertex *tails, 
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

/********************  changestats:  G    ***********/
/*****************
 changestat: d_geodegree
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
 changestat: d_geospartner
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
 changestat: d_gwactor
*****************/
void d_gwactor (int ntoggles, Vertex *heads, Vertex *tails,
		              ModelTerm *mtp, Network *nwp)
{
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex h, t, hd=0, *od;
  int nactors, nevents;

  nactors = (int)mtp->inputparams[0];
  nevents = (nwp->nnodes) - nactors;

//id=nwp->indegree;
  od=nwp->outdegree;
  change = 0.0;
  alpha = mtp->inputparams[1];
  oneexpa = 1.0-exp(-alpha);
  
  *(mtp->dstats) = 0.0;
  for (i=0; i<ntoggles; i++) 
    {
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
      hd = od[h] + (echange - 1)/2;
      if(hd!=0){
        change += echange*(1.0-pow(oneexpa,(double)hd));
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
 changestat: d_gwd
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
 changestat: d_gwdegree706
*****************/
void d_gwdegree706 (int ntoggles, Vertex *heads, Vertex *tails, 
       ModelTerm *mtp, Network *nwp)  {
        /* Slight modification to the parameterization in d_gwdegree */
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
 changestat: d_gwdegreealpha
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
 changestat: d_gwdegreelambda
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
      /* The above line may be an error -- should it be 2.0*echange? */
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
 changestat: d_gwevent
*****************/
void d_gwevent (int ntoggles, Vertex *heads, Vertex *tails,
		              ModelTerm *mtp, Network *nwp) {
  int i, echange=0;
  double alpha, oneexpa, change;
  Vertex h, t, td=0, *id;
  int nactors, nevents;

  nevents = (int)mtp->inputparams[0];
  nactors = nwp->nnodes - nevents;

  id=nwp->indegree;
//od=nwp->outdegree;
  change = 0.0;
  alpha = mtp->inputparams[1];
  oneexpa = 1.0-exp(-alpha);
  
  for (i=0; i<ntoggles; i++) 
    {
      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
 //   hd = od[h] + id[t] + (echange - 1)/2;
      td = id[t] + (echange - 1)/2;
//      if(hd!=0){
//        change += echange*(1.0-pow(oneexpa,(double)hd));
//      }
      if(td!=0){
        change += echange*(1.0-pow(oneexpa,(double)td));
      }
//Rprintf("h %d t %d td %d echange %d change %f\n", h, t, td, echange,change);
//  Rprintf(" h %d t %d hd %d echange %d change %f\n", h, t, hd, echange,change);
//  Rprintf(" od[h] %d id[h] %d od[t] %d id[t] %d\n", od[h], id[h], od[t], id[t]);
      
      if (i+1 < ntoggles)
	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
    }

  *(mtp->dstats) = change*exp(alpha);

//  Rprintf("alpha  %f hd %d td %d change %f\n", alpha, hd, td, *(mtp->dstats));
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/********************  changestats:  H    ***********/
/*****************
 changestat: d_heideriandynamic
*****************/
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

/*****************
 changestat: d_hiertriad
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
 changestat: d_hiertriaddegree
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

/*****************
  numposthree:  called by d_hiertriad*
*****************/
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

/********************  changestats:  I    ***********/
/*****************
 changestat: d_icvar
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
 changestat: d_idc
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
 changestat: d_intransitivedynamic
*****************/
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

/*****************
 changestat: d_intransitivity
*****************/
void d_intransitivity (int ntoggles, Vertex *heads, Vertex *tails, 
	             ModelTerm *mtp, Network *nwp) 
{
  int i, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex node3, h, t;

  *(mtp->dstats) = 0.0;

  if (nwp->directed_flag) {
// directed version
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
          for (node3=1; node3 <= nwp->nnodes; node3++)
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
          t012 = t012 + (nwp->nnodes - 2);  

//        t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
//        t003 = t003-(t111U+t111D+t021C+t021U+t021D+t102+t012);
	b = t021C+t030C+t111D+t111U+t120C+t201+t210;
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;

      if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */ 
    }
    }else{
//  undirected
    for (i=0; i<ntoggles; i++) 
     {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= nwp->nnodes; node3++)
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
          t102 = t102 + (nwp->nnodes - 2);  

        t003 = -(t102+t201+t300);
	b = t300; 
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;
  
     if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */ 
     } // i loop
    }
    
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/********************  changestats:  K    ***********/
/*****************
 changestat: d_kappa
*****************/
void d_kappa (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp)  {
  int i, j, echange=0;
  double nedges, change, ir0, fr0;
  Vertex h, t, hd, td=0, ik2, fk2, nnodes, *id, *od;
  TreeNode *oe=nwp->outedges;
  
  id=nwp->indegree;
  od=nwp->outdegree;
  nnodes = nwp->nnodes;
  
  change = 0.0;
  for (i=0; i<ntoggles; i++) {      
    echange = (EdgetreeSearch(h=heads[i], t=tails[i], oe) == 0) ? 1 : -1;
    ik2=0;
    for (j=1; j<=nnodes; j++) {      
      fk2 = od[j] + id[j];
      ik2 += fk2*(fk2-1);
    }
    hd = od[h] + id[h] + (echange-1)/2;
    td = od[t] + id[t] + (echange-1)/2;
    fk2 = ik2 + echange*2*(hd+td);
    nedges = (double)(nwp->nedges);
    ir0 = (nwp->nedges==0) ? 0.0 : (ik2*0.5/nedges);
    fr0 = (((nwp->nedges)+echange)==0) ? 0.0 : (fk2*0.5/(nedges+echange));
    change += fr0 - ir0;
//   Rprintf("h %d t %d nnodes %d nedges %f ik2 %d fk2 %d ir0 %f fr0 %f change %f\n",h,t, nnodes,  nedges, ik2, fk2, ir0, fr0, change);
      
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  *(mtp->dstats) = change;
  
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/*****************
 changestat: d_monopolymixmat
 ...has nothing to do with the board game; should be read 
 "mono, poly mixing matrix"
*****************/
void d_monopolymixmat(int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp) {
  int edgeflag, i;
  Edge e;
  Vertex h, t, Fdeg, Mdeg, otherF, otherM;
  /* m/p means monogamous/polygamous; F/M means female/male */
  int mFmM, mFpM, pFmM;   
  /* NB: pFpM would be redundant since the total of all 4 is #edges */
  Vertex *od=nwp->outdegree, *id=nwp->indegree;
  TreeNode *oe = nwp->outedges, *ie = nwp->inedges;

  mtp->dstats[0] = mtp->dstats[1] = mtp->dstats[2] = 0.0;
  for (i=0; i < ntoggles; i++) {
    edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
    Fdeg = od[h];
    Mdeg = id[t];
    /* Calculate contribution from change of (F,M) edge only */
    mFmM = (Fdeg==0 && Mdeg==0) - (Fdeg==1 && Mdeg==1 && edgeflag);
    mFpM = (Fdeg==0 && Mdeg>0) - (Fdeg==1 && Mdeg>1 && edgeflag); 
    pFmM = (Mdeg==0 && Fdeg>0) - (Mdeg==1 && Fdeg>1 && edgeflag);
    /* Now calculate contribution from other partners of F or M */
    if(Fdeg - edgeflag == 1) {/* Only case that concerns us */
      for(e = EdgetreeMinimum(oe, h);
      (otherM = oe[e].value) != 0 && otherM == t;
      e = EdgetreeSuccessor(oe, e)); /* This finds otherM */
      if (id[otherM] > 1) {
        mFpM += (Fdeg==1 ? -1 : 1);
      } else {
        mFmM += (Fdeg==1 ? -1 : 1);
        pFmM += (Fdeg==1 ? 1 : -1);
      }
    }
    if(Mdeg - edgeflag == 1) {/* Similarly for Mdeg */
      for(e = EdgetreeMinimum(ie, t);
      (otherF = ie[e].value) != 0 && otherF == h;
      e = EdgetreeSuccessor(ie, e)); /* This finds otherF */
      if (od[otherF] > 1) { /* otherF is poly */
        pFmM += (Mdeg==1 ? -1 : 1);
      } else { /*otherF is mono */
        mFmM += (Mdeg==1 ? -1 : 1);
        mFpM += (Mdeg==1 ? 1 : -1);
      }
    }
    mtp->dstats[0] += (double) mFmM;
    mtp->dstats[1] += (double) mFpM;
    mtp->dstats[2] += (double) pFmM;    
    if (i+1 < ntoggles)
      ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
  }
  i--; 
  while (--i >= 0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

/********************  changestats:  L    ***********/

/********************  changestats:  M    ***********/

/********************  changestats:  N    ***********/

/********************  changestats:  O    ***********/

/********************  changestats:  R    ***********/

/********************  changestats:  S    ***********/
/*****************
 changestat: d_simmeliandynamic
*****************/
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

/*****************
 changestat: d_spatial
  Changescores for inhomogeneous Bernoulli graphs
*****************/                            
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

/********************  changestats:  T    ***********/
/*****************
 changestat: d_transitivedynamic
*****************/
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

/*****************
 changestat: d_transitivity
*****************/
void d_transitivity (int ntoggles, Vertex *heads, Vertex *tails, 
	             ModelTerm *mtp, Network *nwp) 
{
  int i, edgeflag, a, b, c, d, e, edgecount, t300, 
  t210, t120C, t120U, t120D, t201, t030C, t030T, t111U, 
  t111D, t021C, t021U, t021D, t102, t012, t003;
  Vertex node3, h, t;

  *(mtp->dstats) = 0.0;

  if (nwp->directed_flag) {
// directed version
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
          for (node3=1; node3 <= nwp->nnodes; node3++)
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
          t012 = t012 + (nwp->nnodes - 2);  

        t003 = -(t300+t210+t120C+t120U+t120D+t201+t030C+t030T);
        t003 = t003-(t111U+t111D+t021C+t021U+t021D+t102+t012);
	b = t003+t012+t021U+t021D+t102+t030T+t120U+t120D+t300;
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;

      if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */ 
    }
    }else{
//  undirected
    for (i=0; i<ntoggles; i++) 
     {      
      edgeflag = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
      t300 = 0; t201 = 0; t102 = 0; t012 = 0;

      if ((EdgetreeMinimum(nwp->outedges, t) != 0) || 
          (EdgetreeMinimum(nwp->inedges, t) != 0) || 
          (EdgetreeMinimum(nwp->outedges, h) != 0) ||
          (EdgetreeMinimum(nwp->inedges, h) != 0)) {      

          /* ****** loop through node3 ****** */
          for (node3=1; node3 <= nwp->nnodes; node3++)
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
          t102 = t102 + (nwp->nnodes - 2);  

        t003 = (t102+t201+t300);
	b = -t300; 
	*(mtp->dstats) += edgeflag ? -(double)b : (double)b;
  
     if (i+1 < ntoggles)
	  ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */ 
     } // i loop
    }
    
  i--; 
  while (--i>=0)  /*  Undo all previous toggles. */
    ToggleEdge(heads[i], tails[i], nwp); 
}

///*****************
// changestat: d_nearsimmelian
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


// *****************
// changestat: d_actor
// *****************
//void d_actor (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  int i, j, echange;
//  Vertex h, t, deg, nevents, nactors;
//  
//  nactors = (Vertex)inp->inputparams[0];
//  nevents = (nwp->nnodes) - nactors;
////  Rprintf("nevents %d\n", nevents);
////  Rprintf("nactors %d\n", nactors);
////  Rprintf("inp->nstats %d\n", inp->nstats);
////  Rprintf("ntggles %d\n", ntoggles);
//
//  for (i=0; i < inp->nstats; i++) 
//    inp->dstats[i] = 0.0;
//  
//  for (i=0; i<ntoggles; i++)
//    {      
//      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
//      j=0;
//      deg = (Vertex)inp->inputparams[j+1];
//      while(deg != h && j < inp->nstats){
//       j++;
//       deg = (Vertex)inp->inputparams[j+1];
//      }
////  Rprintf("deg %d j %d\n", deg, j);
//      if(j < inp->nstats){inp->dstats[j] += echange;}
////  Rprintf("heads[i] %d tails[i] %d echange %d j %d deg %d\n", heads[i], tails[i], echange,j,inp->dstats[j]);
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
// *****************
// changestat: d_actor
// *****************
//void d_event (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  int i, j, echange;
//  Vertex h, t, deg, nevents, nactors;
//  
//  nactors = (Vertex)inp->inputparams[0];
//  nevents = (nwp->nnodes) - nactors;
////  Rprintf("nevents %d\n", nevents);
////  Rprintf("nactors %d\n", nactors);
//  
//  for (i=0; i < inp->nstats; i++) 
//    inp->dstats[i] = 0.0;
//  
//  for (i=0; i<ntoggles; i++)
//    {      
//      echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
//      j=0;
//      deg = (Vertex)inp->inputparams[j+1];
//      while(deg != t && j < inp->nstats){
//       j++;
//       deg = (Vertex)inp->inputparams[j+1];
//      }
//      if(j < inp->nstats){inp->dstats[j] += echange;}
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//



///*****************
// changestat: d_esa
//*****************/
//void d_esa (int ntoggles, Vertex *heads, Vertex *tails, 
//	      struct OptionInput *inp, Gptr g) 
//{
//  Edge e, f;
//  int i, j, echange;
//  int L2hu, L2ut;
//  Vertex deg;
//  Vertex h, t, u, v;
//  int nevents, nactors;
//
//  nactors = (int)inp->inputparams[0];
//  nevents = (nwp->nnodes) - nactors;
//
//  for (i=0; i < inp->nstats; i++) 
//    inp->dstats[i] = 0.0;
//  
//    for (i=0; i<ntoggles; i++){      
//     echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
//     //
//     // Next for actor shared event counts
//     //
////   // step through inedges of t
////   for(e = EdgetreeMinimum(nwp->inedges, t);
////       (u = nwp->inedges[e].value) != 0;
////       e = EdgetreeSuccessor(nwp->inedges, e)){
////	 if (u != h){
////          L2hu=0;
////          // step through outedges of u
////          for(f = EdgetreeMinimum(nwp->outedges, u);
////	      (v = nwp->outedges[f].value) != 0;
////           f = EdgetreeSuccessor(nwp->outedges, f)){
////               if(EdgetreeSearch(h,v,nwp->outedges)!= 0) L2hu++;
////	  }
////          for(j = 0; j < inp->nstats; j++){
////	    deg = (Vertex)inp->inputparams[j+1];
////            inp->dstats[j] += ((L2hu + echange == deg)
////			     - (L2hu == deg));
////	  }
////         }
////   }
////
////   Next for event shared actor counts
////
//   // step through outedges of h 
//   for(e = EdgetreeMinimum(nwp->outedges, h);
//       (u = nwp->outedges[e].value) != 0;
//       e = EdgetreeSuccessor(nwp->outedges, e)){
//	 if (u != t){
//          L2ut=0;
//          // step through inedges of u
//          for(f = EdgetreeMinimum(nwp->inedges, u); 
//	   (v = nwp->inedges[f].value) != 0;
//           f = EdgetreeSuccessor(nwp->inedges, f)){
//            if(EdgetreeSearch(v,t,nwp->outedges)!= 0) L2ut++;
//	  }
//          for(j = 0; j < inp->nstats; j++){
//	    deg = (Vertex)inp->inputparams[j+1];
//            inp->dstats[j] += ((L2ut + echange == deg)
//			     - (L2ut == deg));
//	  }
//         }
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
//
///*****************
// changestat: d_ase
//*****************/
//void d_ase (int ntoggles, Vertex *heads, Vertex *tails, 
//	      struct OptionInput *inp, Gptr g) 
//{
//  Edge e, f;
//  int i, j, echange;
//  int L2hu, L2ut;
//  Vertex deg;
//  Vertex h, t, u, v;
//  int nevents, nactors;
//
//  nactors = (int)inp->inputparams[0];
//  nevents = (nwp->nnodes) - nactors;
//
//  for (i=0; i < inp->nstats; i++) 
//    inp->dstats[i] = 0.0;
//  
//    for (i=0; i<ntoggles; i++){      
//     echange = (EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) == 0) ? 1 : -1;
//     //
//     // Next for actor shared event counts
//     //
//     // step through inedges of t
//     for(e = EdgetreeMinimum(nwp->inedges, t);
//         (u = nwp->inedges[e].value) != 0;
//         e = EdgetreeSuccessor(nwp->inedges, e)){
//  	 if (u != h){
//            L2hu=0;
//            // step through outedges of u
//            for(f = EdgetreeMinimum(nwp->outedges, u);
//  	      (v = nwp->outedges[f].value) != 0;
//             f = EdgetreeSuccessor(nwp->outedges, f)){
//                 if(EdgetreeSearch(h,v,nwp->outedges)!= 0) L2hu++;
//  	    }
//            for(j = 0; j < inp->nstats; j++){
//  	      deg = (Vertex)inp->inputparams[j+1];
//              inp->dstats[j] += ((L2hu + echange == deg)
//  			     - (L2hu == deg));
//  	    }
//         }
//    }
//     
//   if (i+1 < ntoggles)
//     ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//  }
//
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_bimix
//*****************/
//void d_bimix (int ntoggles, Vertex *heads, Vertex *tails, 
//	        struct OptionInput *inp, Gptr g) {
//
//  int matchvalh, matchvalt, matchswap;
//  Vertex h, t, nnodes, nstats;
//  int i, j, edgeflag=0;
//
//  nstats = inp->nstats;
//  nnodes = (inp->ninputparams)-2*nstats;
////  Rprintf("ninput %d nstats %d\n", inp->ninputparams, inp->nstats);
////  Rprintf("nodes %d\n", nnodes);
//// for (h=0; h<nnodes; h++){
////   Rprintf("match h %d att h %f\n", 
////	    h, inp->attrib[h+nstats]);
//// }
//// for (h=0; h<nstats; h++){
////   Rprintf("match h %d in h %f in t %f\n", 
////	    h, inp->inputparams[h], inp->inputparams[h+nstats]);
//// }
//  for (i=0; i < nstats; i++) 
//    inp->dstats[i] = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      h=heads[i];
//      t=tails[i];
//      matchvalh = inp->attrib[h-1+nstats];
//      matchvalt = inp->attrib[t-1+nstats];
////     if(matchvalt < matchvalh){
////       matchswap = matchvalh;
////       matchvalh = matchvalt;
////       matchvalt = matchswap;
////      }       
//      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
//      for (j=0; j<nstats; j++) 
//	  {
//           if(matchvalh==inp->inputparams[nstats+j] &&
//	      matchvalt==inp->inputparams[       j]
//	     ){inp->dstats[j] += edgeflag ? -1.0 : 1.0;}
//	  }
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_dynamic
//*****************/
//void d_dynamic (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  Vertex h, t, hh, ht;
//  int i, k, nprevedge, edgeflag, discord, lookmore;
//  int nevents, nactors;
//  int nmat, nnodes;
//  double change=0.0;
//
//  nprevedge = (int)((inp->inputparams[0]));
//  nnodes = (nwp->nnodes);
//  nactors = (int)(inp->inputparams[1]);
//  nevents = nnodes - nactors;
//
//  nmat = nevents*nactors;
//
////  Rprintf("nprevedge %d\n", nprevedge);
////  Rprintf("nactors %d nevents %d\n", nactors, nevents);
////  for (k=1; k<=nprevedge; k++) {
////      Rprintf("k %d x0.h %d x0.t %d\n", k,
////	      (Vertex)(inp->attrib[          k]),
////	      (Vertex)(inp->attrib[nprevedge+k]));
////  }
//
//  for (i=0; i < inp->nstats; i++)
//   inp->dstats[i] = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      /*Get the initial state of the edge and its reflection*/
//      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
//      if(h > t){
//	hh = t;
//	ht = h;
//      }else{
//	hh = h;
//	ht = t;
//      }
//      discord = edgeflag ? 1 : -1;
//
//      k=1;
//      lookmore=1;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,hh,ht);
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      hh,ht,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//       if(ht == (Vertex)(inp->attrib[            k]) &&
//          hh == (Vertex)(inp->attrib[nprevedge + k])
//	 ){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//
//           /*If the proposed edge existed in x0, get dissolution rate */
//	   /* lookmore is location of (hh, ht) in the dissolve matrix */
//           lookmore=2*nprevedge+nmat+(hh-1)+(ht-1-nactors)*nevents;
//	   /* change is NEGATIVE the dissolve rate if an edge exists */
////       Rprintf("hh %d ht %d lookmore %d att %f\n",hh, ht, lookmore, (double)(inp->attrib[lookmore]));
//           change=-discord*(double)(inp->attrib[lookmore]);
////       Rprintf("lookmore %d change %f\n",lookmore, change);
//	   lookmore=0;
//           /*Update the change statistics, as appropriate*/
//	   /*change is the number of edges dissolved*/
//           inp->dstats[0] += change;
////      if(change!=0){
////       Rprintf("h %d t %d hh %d ht %d edgeflag %d beta %f\n",
////	      h,t,hh,ht, edgeflag, (double)(inp->attrib[lookmore]));
////       Rprintf("edgeflag %d change %f discord %d\n",
////	      edgeflag, change, discord);
////      }
//          }else{
//		   ++k;
//	  }
//      }
//
//      if(lookmore && nprevedge > 0){
//       /*If the proposed edge existed in x0, get dissolution rate */
//       /* lookmore, so no edge exisited in x0, Get formation rate */
//       /* lookmore is location of (hh, ht) in the form matrix */
//       lookmore=2*nprevedge+(hh-1-nactors)+(ht-1)*nevents;
//       /* change is the form number of formed edges */
//       change=-discord*(double)(inp->attrib[lookmore]);
//       /*Update the change statistics, as appropriate*/
//       inp->dstats[1] += change;
//      }
//
////      if(lookmore && change>0){
////       Rprintf("h %d t %d hh %d ht %d edgeflag %d alpha\n",
////	      h,t,hh,ht,edgeflag);
////       Rprintf("edgeflag %d change %f discord %d\n",
////	      edgeflag, change, discord);
////      }
////      /*Update the change statistics, as appropriate*/
////      *(inp->dstats) += change;
//
////  Rprintf("heads[i] %d tails[i] %d change %f \n", 
////		      heads[i], tails[i], change);
//      
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_biduration
//*****************/
//void d_biduration (int ntoggles, Vertex *heads, Vertex *tails, 
//	        struct OptionInput *inp, Gptr g) {
//
//  int matchvalh, matchvalt;
//  int k, nprevedge, discord, lookmore;
//  Vertex h, t, nnodes, nstats;
//  Vertex hh, ht;
//  int i, j, edgeflag=0;
//
//  nstats = inp->nstats;
////  nnodes = (inp->ninputparams)-2*nstats;
//  nnodes = (nwp->nnodes);
//  nprevedge = (int)(inp->inputparams[0]);
//  Rprintf("ninput %d nstats %d\n", inp->ninputparams, inp->nstats);
//  Rprintf("nodes %d nprevedge %d\n", nnodes, nprevedge);
//// for (h=1; h<=nnodes; h++){
////   Rprintf("match h %d att h %f\n", 
////	    h, inp->attrib[h+nstats]);
//// }
//// for (h=0; h<nstats; h++){
////   Rprintf("match h %d in h %f in t %f\n", 
////	    h, inp->inputparams[h], inp->inputparams[h+nstats]);
//// }
//  for (i=0; i < nstats; i++) 
//    inp->dstats[i] = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      h=heads[i];
//      t=tails[i];
//      matchvalh = inp->attrib[h+nstats];
//      matchvalt = inp->attrib[t+nstats];
////     if(matchvalt < matchvalh){
////       matchswap = matchvalh;
////       matchvalh = matchvalt;
////       matchvalt = matchswap;
////      }       
//      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
//      // 1 = edge removed -1 = edge added
//      discord = edgeflag ? 1 : -1;
//      k=1;
//      lookmore=1;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d h %d t %d\n",
////	        h,t,h,t);
////       Rprintf("h %d t %d h %d t %d\n",
////	      t,h,(Vertex)(inp->attrib[            k]),
////	            (Vertex)(inp->attrib[nprevedge + k]));
////     Check to see if it was a previous edge
//       if(t == (Vertex)(inp->attrib[nstats+nnodes           + k]) &&
//          h == (Vertex)(inp->attrib[nstats+nnodes+nprevedge + k])
//	 ){
////        Rprintf("h %d t %d h attr %d t attr %d\n",
////	           h,t,(Vertex)(inp->attrib[nstats+nnodes+nprevedge +   k]),
////	          (Vertex)(inp->attrib[nstats+nnodes+ k]));
////        Check to see which mixing group it is in
//          j=0;
//          while(lookmore && j < nstats){
//            if(matchvalh==inp->inputparams[nstats+j+1] &&
//	       matchvalt==inp->inputparams[      +j+1]
//	      ){
//
//              /*If the proposed edge existed in x0, get dissolution rate */
//	      /* lookmore is location of (h, t) in the dissolve matrix */
////            lookmore=2*nprevedge+(h-1)+(t-1-nactors)*nevents;
//	      /* change is NEGATIVE the dissolve rate if an edge exists */
////          Rprintf("h %d t %d lookmore %d att %f\n",h, t, lookmore, (double)(inp->attrib[lookmore]));
////              change=-discord*(double)(inp->attrib[lookmore]);
//           /*Update the change statistics, as appropriate*/
//	   /*change is the number of edges dissolved*/
//	        inp->dstats[j] -= discord;
//	        lookmore=0;
//	       }
//	       j++;
//	  }
//	  lookmore=0;
//       }
//       k++;
//      }
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_endure
//*****************/
//void d_endure (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  Vertex h, t, hh, ht;
//  int i, k, nprevedge, edgeflag, discord, lookmore;
//  int nevents, nactors;
//  int nnodes;
//  double change=0.0;
//
//  nprevedge = (int)((inp->inputparams[0]));
//  nnodes = (nwp->nnodes);
//  nactors = (int)(inp->inputparams[1]);
//  nevents = nnodes - nactors;
//
////  Rprintf("nprevedge %d\n", nprevedge);
////  Rprintf("nactors %d nevents %d\n", nactors, nevents);
////  for (k=1; k<=nprevedge; k++) {
////      Rprintf("k %d x0.h %d x0.t %d\n", k,
////	      (Vertex)(inp->attrib[          k]),
////	      (Vertex)(inp->attrib[nprevedge+k]));
////  }
//
//  *(inp->dstats) = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      /*Get the initial state of the edge and its reflection*/
//      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
//      if(h > t){
//	hh = t;
//	ht = h;
//      }else{
//	hh = h;
//	ht = t;
//      }
//      // 1 = edge removed -1 = edge added
//      discord = edgeflag ? 1 : -1;
//
//// Rprintf("nprevedge %d\n", nprevedge);
//
//      k=1;
//      lookmore=1;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,hh,ht);
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      ht,hh,(Vertex)(inp->attrib[            k]),
////	            (Vertex)(inp->attrib[nprevedge + k]));
//       if(ht == (Vertex)(inp->attrib[            k]) &&
//          hh == (Vertex)(inp->attrib[nprevedge + k])
//	 ){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//
//           /*If the proposed edge existed in x0, get dissolution rate */
//	   /* lookmore is location of (hh, ht) in the dissolve matrix */
//	   /* change is NEGATIVE the dissolve rate if an edge exists */
////       Rprintf("hh %d ht %d lookmore %d att %f\n",hh, ht, lookmore, (double)(inp->attrib[lookmore]));
////           change=-discord*(double)(inp->attrib[lookmore]);
//           change=-discord;
////       Rprintf("lookmore %d change %f discord %d \n",lookmore, change, discord);
//	   lookmore=0;
//           /*Update the change statistics, as appropriate*/
//	   /*change is the number of edges dissolved*/
//           *(inp->dstats) += change;
////      if(change!=0){
////       Rprintf("h %d t %d hh %d ht %d edgeflag %d beta %f\n",
////	      h,t,hh,ht, edgeflag, (double)(inp->attrib[lookmore]));
////       Rprintf("edgeflag %d change %f discord %d\n",
////	      edgeflag, change, discord);
////      }
//          }else{
//		   ++k;
//	  }
//      }
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_bichange
//*****************/
//void d_bichange (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  Vertex h, t, hh, ht;
//  int i, k, nprevedge, edgeflag, discord, lookmore;
//  int nevents, nactors;
//  int nnodes;
//  double change=0.0;
//
//  nprevedge = (int)((inp->inputparams[0]));
//  nnodes = (nwp->nnodes);
//  nactors = (int)(inp->inputparams[1]);
//  nevents = nnodes - nactors;
//
////  Rprintf("nprevedge %d\n", nprevedge);
////  Rprintf("nactors %d nevents %d\n", nactors, nevents);
////  for (k=1; k<=nprevedge; k++) {
////      Rprintf("k %d x0.h %d x0.t %d\n", k,
////	      (Vertex)(inp->attrib[          k]),
////	      (Vertex)(inp->attrib[nprevedge+k]));
////  }
//
//  *(inp->dstats) = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      /*Get the initial state of the edge and its reflection*/
//      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
//      if(h > t){
//	hh = t;
//	ht = h;
//      }else{
//	hh = h;
//	ht = t;
//      }
//
//      k=1;
//      lookmore=1;
//      change=0.0;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,hh,ht);
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      hh,ht,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//       if(ht == (Vertex)(inp->attrib[            k]) &&
//          hh == (Vertex)(inp->attrib[nprevedge + k])
//	 ){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//
//           /*If the proposed edge existed in x0, get dissolution rate */
//	   /* lookmore is location of (hh, ht) in the dissolve matrix */
//           lookmore=2*nprevedge+(hh-1)+(ht-1-nactors)*nevents;
//	   /* change is NEGATIVE the dissolve rate if an edge exists */
////       Rprintf("hh %d ht %d lookmore %d att %f\n",hh, ht, lookmore, (double)(inp->attrib[lookmore]));
//           change+=edgeflag*(double)(inp->attrib[lookmore]);
////       Rprintf("lookmore %d change %f\n",lookmore, change);
//	   lookmore=0;
////      if(change!=0){
////       Rprintf("h %d t %d hh %d ht %d edgeflag %d beta %f\n",
////	      h,t,hh,ht, edgeflag, (double)(inp->attrib[lookmore]));
////       Rprintf("edgeflag %d change %f discord %d\n",
////	      edgeflag, change, discord);
////      }
//          }else{
//		   ++k;
//	  }
//
//      if(lookmore && edgeflag){
//       /* lookmore, so no edge exisited in x0, Get formation rate */
//       /* lookmore is location of (hh, ht) in the form matrix */
//       lookmore=2*nprevedge+(hh-1)+(ht-1-nactors)*nevents;
//       /* change is the form number of formed edges */
//       change+=(1-edgeflag)*(double)(inp->attrib[lookmore]);
//      }
//      /*Update the change statistics, as appropriate*/
//      /*change is the number of edges dissolved*/
//      *(inp->dstats) -= change;
//      }
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_step
//*****************/
//void d_step (int ntoggles, Vertex *heads, Vertex *tails, 
//	     struct OptionInput *inp, Gptr g) {
//
//  int matchvalh, matchvalt;
//  int k, nprevedge, discord, lookmore;
//  Vertex h, t, nnodes, nstats;
//  Vertex hh, ht;
//  Vertex ch, ct, cv, ph, pt, pv;
//  int i, j, edgeflag=0;
//
//  nstats = inp->nstats;
////  nnodes = (inp->ninputparams)-2*nstats;
//  ph1 = *(nwp->ph1);
//  pt1 = *(nwp->pt1);
//  pv1 = *(nwp->pv1);
//  ph2 = *(nwp->ph2);
//  pt2 = *(nwp->pt2);
//  pv2 = *(nwp->pv2);
//  if(pv1 != (EdgetreeSearch(ph1, pt1, nwp->outedges) != 0)
//  && pv1 != (EdgetreeSearch(ph1, pt1, nwp->outedges) != 0)){
//	  // Proposal made so update current
//   *(nwp->ch1) = ph1;
//   *(nwp->ct1) = pt1;
//   *(nwp->cv1) = pv1;
//   *(nwp->ch2) = ph2;
//   *(nwp->ct2) = pt2;
//   *(nwp->cv2) = pv2;
//  }
//  ch1 = *(nwp->ch1);
//  ct1 = *(nwp->ct1);
//  cv1 = *(nwp->cv1);
//  ch2 = *(nwp->ch2);
//  ct2 = *(nwp->ct2);
//  cv2 = *(nwp->cv2);
//
//  nnodes = (nwp->nnodes);
//  nprevedge = (int)(inp->inputparams[0]);
//  Rprintf("ninput %d nstats %d\n", inp->ninputparams, inp->nstats);
//  Rprintf("nodes %d nprevedge %d\n", nnodes, nprevedge);
//// for (h=1; h<=nnodes; h++){
////   Rprintf("match h %d att h %f\n", 
////	    h, inp->attrib[h+nstats]);
//// }
//// for (h=0; h<nstats; h++){
////   Rprintf("match h %d in h %f in t %f\n", 
////	    h, inp->inputparams[h], inp->inputparams[h+nstats]);
//// }
//  for (i=0; i < nstats; i++) 
//    inp->dstats[i] = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      h=heads[i];
//      t=tails[i];
//      edgeflag=(EdgetreeSearch(h, t, nwp->outedges) != 0);
//      // 1 = edge removed -1 = edge added
//      discord = edgeflag ? 1 : -1;
//      k=1;
//      lookmore=1;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d h %d t %d\n",
////	        h,t,h,t);
////       Rprintf("h %d t %d h %d t %d\n",
////	      t,h,(Vertex)(inp->attrib[            k]),
////	            (Vertex)(inp->attrib[nprevedge + k]));
////     Check to see if it was a previous edge
//       if(t == (Vertex)(inp->attrib[nstats+nnodes           + k]) &&
//          h == (Vertex)(inp->attrib[nstats+nnodes+nprevedge + k])
//	 ){
////        Rprintf("h %d t %d h attr %d t attr %d\n",
////	           h,t,(Vertex)(inp->attrib[nstats+nnodes+nprevedge +   k]),
////	          (Vertex)(inp->attrib[nstats+nnodes+ k]));
////        Check to see which mixing group it is in
//          j=0;
//          while(lookmore && j < nstats){
//            if(matchvalh==inp->inputparams[nstats+j+1] &&
//	       matchvalt==inp->inputparams[      +j+1]
//	      ){
//
//              /*If the proposed edge existed in x0, get dissolution rate */
//	      /* lookmore is location of (h, t) in the dissolve matrix */
////            lookmore=2*nprevedge+(h-1)+(t-1-nactors)*nevents;
//	      /* change is NEGATIVE the dissolve rate if an edge exists */
////          Rprintf("h %d t %d lookmore %d att %f\n",h, t, lookmore, (double)(inp->attrib[lookmore]));
////              change=-discord*(double)(inp->attrib[lookmore]);
//           /*Update the change statistics, as appropriate*/
//	   /*change is the number of edges dissolved*/
//	        inp->dstats[j] -= discord;
//	        lookmore=0;
//	       }
//	       j++;
//	  }
//	  lookmore=0;
//       }
//       k++;
//      }
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
///*****************
// changestat: d_biduration_old
//*****************/
//void d_biduration_old (int ntoggles, Vertex *heads, Vertex *tails, 
//	         struct OptionInput *inp, Gptr g) 
//{
//  Vertex h, t, hh, ht;
//  int i, k, nprevedge, edgeflag, discord, lookmore;
//  int nevents, nactors;
//  int nnodes;
//  double change=0.0;
//
//  nprevedge = (int)((inp->inputparams[0]));
//  nnodes = (nwp->nnodes);
//  nactors = (int)(inp->inputparams[1]);
//  nevents = nnodes - nactors;
//
////  Rprintf("nprevedge %d\n", nprevedge);
////  Rprintf("nactors %d nevents %d\n", nactors, nevents);
////  for (k=1; k<=nprevedge; k++) {
////      Rprintf("k %d x0.h %d x0.t %d\n", k,
////	      (Vertex)(inp->attrib[          k]),
////	      (Vertex)(inp->attrib[nprevedge+k]));
////  }
//
//  *(inp->dstats) = 0.0;
//
//  for (i=0; i<ntoggles; i++) 
//    {
//      /*Get the initial state of the edge and its reflection*/
//      edgeflag=(EdgetreeSearch(h=heads[i], t=tails[i], nwp->outedges) != 0);
//      if(h > t){
//	hh = t;
//	ht = h;
//      }else{
//	hh = h;
//	ht = t;
//      }
//      // 1 = edge removed -1 = edge added
//      discord = edgeflag ? 1 : -1;
//
//// Rprintf("nprevedge %d\n", nprevedge);
//
//      k=1;
//      lookmore=1;
//      while(lookmore && k <= nprevedge){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,hh,ht);
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      ht,hh,(Vertex)(inp->attrib[            k]),
////	            (Vertex)(inp->attrib[nprevedge + k]));
//       if(ht == (Vertex)(inp->attrib[            k]) &&
//          hh == (Vertex)(inp->attrib[nprevedge + k])
//	 ){
////       Rprintf("h %d t %d hh %d ht %d\n",
////	      h,t,(Vertex)(inp->attrib[            k]),
////	          (Vertex)(inp->attrib[nprevedge + k]));
//
//           /*If the proposed edge existed in x0, get dissolution rate */
//	   /* lookmore is location of (hh, ht) in the dissolve matrix */
//           lookmore=2*nprevedge+(hh-1)+(ht-1-nactors)*nevents;
//	   /* change is NEGATIVE the dissolve rate if an edge exists */
////       Rprintf("hh %d ht %d lookmore %d att %f\n",hh, ht, lookmore, (double)(inp->attrib[lookmore]));
////           change=-discord*(double)(inp->attrib[lookmore]);
//           change=-discord;
////       Rprintf("lookmore %d change %f discord %d \n",lookmore, change, discord);
//	   lookmore=0;
//           /*Update the change statistics, as appropriate*/
//	   /*change is the number of edges dissolved*/
//           *(inp->dstats) += change;
////      if(change!=0){
////       Rprintf("h %d t %d hh %d ht %d edgeflag %d beta %f\n",
////	      h,t,hh,ht, edgeflag, (double)(inp->attrib[lookmore]));
////       Rprintf("edgeflag %d change %f discord %d\n",
////	      edgeflag, change, discord);
////      }
//          }else{
//		   ++k;
//	  }
//      }
//
//      if (i+1 < ntoggles)
//	ToggleEdge(heads[i], tails[i], nwp);  /* Toggle this edge if more to come */
//    }
//  i--; 
//  while (--i>=0)  /*  Undo all previous toggles. */
//    ToggleEdge(heads[i], tails[i], nwp); 
//}
//
 

