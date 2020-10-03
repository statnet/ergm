/*  File src/changestats_dgw_sp.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */

#include "changestats_dgw_sp.h"

/**************************
 dsp Calculation functions
**************************/

/*
Changescore for ESPs based on two-paths in undirected graphs i.e. configurations for edge i<->j such that i<->k<->j (where <-> here denotes an undirected edge).

UTP:
L2th - count t<->k<->h
L2tk - for each t<->k neq h: k<->h, count u such that k<->u<->h
L2hk - for each h<->k neq t: k<->t, count u such that k<->u<->t

This function will only work properly with undirected graphs, and should only be called in that case.
*/
void dspUTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  //Rprintf("dvec:\n");
  //for(i=0;i<nd;i++)
  //  Rprintf("%d ",(int)dvec[i]);
  //Rprintf("\n");

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head,e,u){
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tu + echange == deg) - (L2tu == deg));
          cs[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head,e,u){
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tu + echange == deg) - (L2tu == deg));
          cs[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);

  //Rprintf("cs:\n");
  //for(i=0;i<nd;i++)
  //  Rprintf("%d ",(int)cs[i]);
  //Rprintf("\n");

}



/*
Changescore for dsps based on outgoing two-paths, i.e. configurations for non-edge i->j such that i->k->j.

This function should only be used in the directed case
*/
void dspOTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2tk, L2kh; /*Two-path counts for t->h, t->k, and k->h edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //ZERO_ALL_CHANGESTATS
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
  cs[i]=0.0;
  
  FOR_EACH_TOGGLE(i) {
    //Rprintf("Working on toggle %d (%d,%d)\n",i,TAIL(i),HEAD(i));
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    
    /* step through outedges of head (i.e., k: t->k)*/
    STEP_THROUGH_OUTEDGES(head, e, k) {
      
      if(k!=tail){ /*Only use contingent cases*/
        L2tk=0;
        /* step through inedges of k, incl. (head,k) itself */
        STEP_THROUGH_INEDGES(k, f, u) {
          L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/
        }
        /*Update the changestat for the t->k edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tk + echange == deg) - (L2tk == deg));
        }
      }
    }
    /* step through inedges of tail (i.e., k: k->t)*/
    STEP_THROUGH_INEDGES(tail, e, k) {
      if (k!=head){ /*Only use contingent cases*/
    
        L2kh=0;
        /* step through outedges of k , incl. (k,tail) itself */
        STEP_THROUGH_OUTEDGES(k, f, u) {
          L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/
        }
        /*Update the changestat for the k->t edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2kh + echange == deg) - (L2kh == deg));
        }
      }
    }
    
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for DSPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.
IE cyclical shared partners
 
ITP:
L2th - count j->k->i
L2hk - for each j->k neq i: k->i, count u such that k->u->j 
L2kt - for each k->i neq j: j->k, count u such that i->u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void dspITP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2hk, L2kt; /*Two-path counts for t->h, h->k, and k->t edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //ZERO_ALL_CHANGESTATS
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through outedges of head (i.e., k: h->k)*/
    STEP_THROUGH_OUTEDGES(head, e, k) {
      if((k!=tail)){ /*Only use contingent cases*/
        /*We have a h->k->t two-path, so add it to our count.*/
        
        L2kt=0;
        /*Now, count # u such that k->u->h (so that we know k's ESP value)*/
        STEP_THROUGH_INEDGES(k, f, u) {
          L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/
        }
        /*Update the changestat for the h->k edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2kt + echange == deg) - (L2kt == deg));
        }
      }
    }
    /* step through inedges of tail (i.e., k: k->t)*/
    STEP_THROUGH_INEDGES(tail, e, k) {
      if((k!=head)){ /*Only use contingent cases*/
        L2hk=0;
        /*Now, count # u such that t->u->k (so that we know k's ESP value)*/
        STEP_THROUGH_OUTEDGES(k, f, u) {
          L2hk+=IS_OUTEDGE(u,head); /*Increment if there is a cyclic edge*/
        }
        /*Update the changestat for the k->t edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2hk + echange == deg) - (L2hk == deg));
        }
      }
    }
    
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for DSPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

OSP:
L2th - count t->k, h->k
L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void dspOSP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2tk; /*Two-path counts for t->h, t->k, and k->t edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //ZERO_ALL_CHANGESTATS
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/
    STEP_THROUGH_INEDGES(head, e, k) {
      if(k!=tail){
        /*Do we have a t->k,h->k SP?  If so, add it to our count.*/
        L2tk=0;
        /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/
        STEP_THROUGH_OUTEDGES(k, f, u) {
          if (u != tail)
            /*Increment if there is an OSP  */
            L2tk+=(IS_OUTEDGE(tail,u));  
        }
        /*Update the changestat for the t->k edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          // twice (one for each direction of t<->k)
          cs[j] += (2*((L2tk + echange == deg) - (L2tk == deg)));
        }
        
      }
    }
    
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

ISP:
L2th - count k->t, k->h
L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void dspISP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2kh; /*Two-path counts for t->h, h->k, and k->h edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //ZERO_ALL_CHANGESTATS
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/
    STEP_THROUGH_OUTEDGES(tail, e, k) {
      if(k!=head){
        L2kh=0;
        /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/
        STEP_THROUGH_INEDGES(k, f, u) {
          if(u!=head)
            L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/
        }
        /*Update the changestat for the k->h edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          // twice (one for each direction of h<->k)
          cs[j] += (2 * ((L2kh + echange == deg) - (L2kh == deg)));
        }
        
      }
    }
    
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on recursive two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j).

RTP:
L2th - count t<->k<->h
L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t
L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t
L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h
L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void dspRTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange, htedge;
  int L2th,L2tk,L2kt,L2hk,L2kh; /*Two-path counts for various edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/
    /* step through inedges of tail (k->t: k!=h,h->t,k<->h)*/
    STEP_THROUGH_INEDGES(tail,e,k){
      if(k!=head){
        /*Do we have a t<->k<->h TP?  If so, add it to our count.*/
        L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head));
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/
          L2kt=0;
          /*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=tail)&&(IS_OUTEDGE(u,k)))
              L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/
          }
          /*Update the changestat for the k->t edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2kt + echange == deg) - (L2kt == deg));
          }
        }
      }
    }
    /* step through outedges of tail (t->k: k!=h,h->t,k<->h)*/
    STEP_THROUGH_OUTEDGES(tail,e,k) {
      if(k!=head){
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/
          L2tk=0;
          /*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=tail)&&(IS_OUTEDGE(u,k)))
              L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/
          }
          /*Update the changestat for the t->k edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2tk + echange == deg) - (L2tk == deg));
          }
        }
      }
    }
    /* step through inedges of head (k->h: k!=t,h->t,k<->t)*/
    STEP_THROUGH_INEDGES(head,e,k){
      if(k!=tail){
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/
          L2kh=0;
          /*Now, count # u such that k<->u<->h (to get k->h's ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=head)&&(IS_OUTEDGE(u,k)))
              L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/
          }
          /*Update the changestat for the k->h edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2kh + echange == deg) - (L2kh == deg));
          }
        }
      }
    }
    /* step through outedges of head (h->k: k!=t,h->t,k<->t)*/
    STEP_THROUGH_OUTEDGES(head,e,k){
      if(k!=tail){
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/
          L2hk=0;
          /*Now, count # u such that k<->u<->h (to get h->k's ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=head)&&(IS_OUTEDGE(u,k)))
              L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/
          }
          /*Update the changestat for the h->k edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2hk + echange == deg) - (L2hk == deg));
          }
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*****************
 changestat: d_dsp
*****************/
/*
Note that d_esp is a meta-function, dispatching actual changescore 
calculation to one of the esp*_calc routines, based on the selected shared 
partner type code.

Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Recursive two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
D_CHANGESTAT_FN(d_ddsp) { 
  int type;
  double *dvec,*cs;
  
  /*Set things up*/
  ZERO_ALL_CHANGESTATS(i);
  type=(int)INPUT_PARAM[0];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/
  cs=CHANGE_STAT;               /*Grab the pointer to the CS vector*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: dspUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPOTP: dspOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPITP: dspITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPRTP: dspRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPOSP: dspOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPISP: dspISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


/*****************
 changestat: d_gwdsp
*****************/

/*
Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Recursive two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
D_CHANGESTAT_FN(d_dgwdsp) { 
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[0];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[1];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[2];   /*Get the max ESP cutoff to use*/
  cs=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/
  dvec=Calloc(maxesp,double);   /*Allocate memory for the ESP vals*/
  for(i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1.0;

  /*Obtain the DSP changescores (by type)*/
  switch(type){
    case ESPUTP: dspUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPOTP: dspOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPITP: dspITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPRTP: dspRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPOSP: dspOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPISP: dspISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
  }
  
  /*Compute the gwdsp changescore*/
  for(i=0;i<maxesp;i++){
    if(cs[i]!=0.0) {                 /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*cs[i];
      //Rprintf("count %f: %f\n", dvec[i], cs[i]);
    }
  }
  CHANGE_STAT[0]*=exp(alpha);
  Free(cs);
  Free(dvec);
}



/**************************
 ESP Calculation functions
**************************/

/*
Changescore for ESPs based on two-paths in undirected graphs i.e. configurations for edge i<->j such that i<->k<->j (where <-> here denotes an undirected edge).

UTP:
L2th - count t<->k<->h
L2tk - for each t<->k neq h: k<->h, count u such that k<->u<->h
L2hk - for each h<->k neq t: k<->t, count u such that k<->u<->t

This function will only work properly with undirected graphs, and should only be called in that case.
*/
void espUTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tu, L2uh;
  Vertex deg;
  Vertex tail, head, u, v;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  //Rprintf("dvec:\n");
  //for(i=0;i<nd;i++)
  //  Rprintf("%d ",(int)dvec[i]);
  //Rprintf("\n");

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through outedges of head */
    STEP_THROUGH_OUTEDGES(head,e,u){
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tu + echange == deg) - (L2tu == deg));
          cs[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    /* step through inedges of head */
    STEP_THROUGH_INEDGES(head,e,u){
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){
        L2th++;
        L2tu=0;
        L2uh=0;
        /* step through outedges of u */
        STEP_THROUGH_OUTEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        /* step through inedges of u */
        STEP_THROUGH_INEDGES(u,f,v){
          if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;
          if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;
        }
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tu + echange == deg) - (L2tu == deg));
          cs[j] += ((L2uh + echange == deg) - (L2uh == deg));
        }
      }
    }
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);

  //Rprintf("cs:\n");
  //for(i=0;i<nd;i++)
  //  Rprintf("%d ",(int)cs[i]);
  //Rprintf("\n");

}


/*
Changescore for ESPs based on outgoing two-paths, i.e. configurations for edge i->j such that i->k->j.

OTP:
L2th - count i->k->j
L2tk - for each i->k neq j: j->k, count u such that i->u->k
L2kh - for each k->j neq i: k->i, count u such that k->u->j

This function should only be used in the directed case, with espUTP being used in the undirected case.
*/
void espOTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tk, L2kh; /*Two-path counts for t->h, t->k, and k->h edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //Rprintf("Clearing changestats\n");
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    //Rprintf("Working on toggle %d (%d,%d)\n",i,TAIL(i),HEAD(i));
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    //Rprintf("\tEdge change is %d\n",echange);
    /* step through outedges of tail (i.e., k: t->k)*/
    //Rprintf("Walking through outedges of tail\n",echange);
    STEP_THROUGH_OUTEDGES(tail,e,k) {
      if((k!=head)&&(IS_OUTEDGE(k,head))){
        /*We have a t->k->h two-path, so add it to our count.*/
        L2th++;
      }
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/
        //Rprintf("\tk==%d, passed criteria\n",k);
        L2tk=0;
        /*Now, count # u such that t->u->k (to find t->k's ESP value)*/
        STEP_THROUGH_INEDGES(k,f,u){
          //Rprintf("\t\tTrying u==%d\n",u);
          if(u!=tail) 
            L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/
        }
        /*Update the changestat for the t->k edge*/
        //Rprintf("\t2-path count was %d\n",L2tk);
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2tk + echange == deg) - (L2tk == deg));
        }
      }
    }
    /* step through inedges of head (i.e., k: k->h)*/
    //Rprintf("Walking through inedges of head\n",k);
    STEP_THROUGH_INEDGES(head,e,k){
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/
        //Rprintf("\tk==%d, passed criteria\n",k);
        L2kh=0;
        /*Now, count # u such that k->u->j (to find k->h's ESP value)*/
        STEP_THROUGH_OUTEDGES(k,f,u){
          if(u!=head) 
            L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/
        }
        /*Update the changestat for the k->t edge*/
        //Rprintf("\t2-path count was %d\n",L2kh);
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2kh + echange == deg) - (L2kh == deg));
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    //Rprintf("2-path count for (t,h) was %d\n",L2th);
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.

ITP:
L2th - count j->k->i
L2hk - for each j->k neq i: k->i, count u such that k->u->j 
L2kt - for each k->i neq j: j->k, count u such that i->u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void espITP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2hk, L2kt; /*Two-path counts for t->h, h->k, and k->t edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  //Rprintf("Clearing changestats\n");
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    //Rprintf("Working on toggle %d (%d,%d)\n",i,TAIL(i),HEAD(i));
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    //Rprintf("\tEdge change is %d\n",echange);
    /* step through outedges of head (i.e., k: h->k)*/
    //Rprintf("Walking through outedges of head\n",echange);
    STEP_THROUGH_OUTEDGES(head,e,k){
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/
        //Rprintf("\tk==%d, passed criteria\n",k);
        /*We have a h->k->t two-path, so add it to our count.*/
        L2th++;
        L2hk=0;
        /*Now, count # u such that k->u->h (so that we know k's ESP value)*/
        STEP_THROUGH_OUTEDGES(k,f,u){
          //Rprintf("\t\tTrying u==%d\n",u);
          if(u!=head) 
            L2hk+=IS_OUTEDGE(u,head); /*Increment if there is a cyclic edge*/
        }
        /*Update the changestat for the h->k edge*/
        //Rprintf("\t2-path count was %d\n",L2hk);
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2hk + echange == deg) - (L2hk == deg));
        }
      }
    }
    /* step through inedges of tail (i.e., k: k->t)*/
    //Rprintf("Walking through inedges of tail\n",k);
    STEP_THROUGH_INEDGES(tail,e,k){
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/
        //Rprintf("\tk==%d, passed criteria\n",k);
        L2kt=0;
        /*Now, count # u such that t->u->k (so that we know k's ESP value)*/
        STEP_THROUGH_INEDGES(k,f,u){
          if(u!=tail) 
            L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/
        }
        /*Update the changestat for the k->t edge*/
        //Rprintf("\t2-path count was %d\n",L2kt);
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2kt + echange == deg) - (L2kt == deg));
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    //Rprintf("2-path count for (t,h) was %d\n",L2th);
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

OSP:
L2th - count t->k, h->k
L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void espOSP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2tk, L2kt; /*Two-path counts for t->h, t->k, and k->t edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/
    STEP_THROUGH_OUTEDGES(tail,e,k) {
      if(k!=head){
        /*Do we have a t->k,h->k SP?  If so, add it to our count.*/
        L2th+=IS_OUTEDGE(head,k);
        if(IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/
          L2tk=0;
          /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if(u!=tail)
              L2tk+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/
          }
          /*Update the changestat for the t->k edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2tk + echange == deg) - (L2tk == deg));
          }
        }
      }
    }
    /* step through inedges of tail (i.e., k: k->t, k->h, k!=h)*/
    STEP_THROUGH_INEDGES(tail,e,k){
      if((k!=head)&&(IS_OUTEDGE(k,head))){ /*Only stats that could change*/
        L2kt=0;
        /*Now, count # u such that t->u,k->u (to get k->t's ESP value)*/
        STEP_THROUGH_OUTEDGES(k,f,u){ 
          if(u!=tail)
            L2kt+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/
        }
        /*Update the changestat for the k->t edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2kt + echange == deg) - (L2kt == deg));
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

ISP:
L2th - count k->t, k->h
L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void espISP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange;
  int L2th, L2hk, L2kh; /*Two-path counts for t->h, h->k, and k->h edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/
    STEP_THROUGH_INEDGES(head,e,k){
      if(k!=tail){
        /*Do we have a k->t,k->h SP?  If so, add it to our count.*/
        L2th+=IS_OUTEDGE(k,tail);
        if(IS_OUTEDGE(tail,k)){ /*Only consider stats that could change*/
          L2kh=0;
          /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/
          STEP_THROUGH_INEDGES(k,f,u){ 
            if(u!=head)
              L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/
          }
          /*Update the changestat for the k->h edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2kh + echange == deg) - (L2kh == deg));
          }
        }
      }
    }
    /* step through outedges of head (i.e., k: h->k, t->k, k!=t)*/
    STEP_THROUGH_OUTEDGES(head,e,k){
      if((k!=tail)&&(IS_OUTEDGE(tail,k))){ /*Only stats that could change*/
        L2hk=0;
        /*Now, count # u such that u->h,u->k (to get k->h's ESP value)*/
        STEP_THROUGH_INEDGES(k,f,u){ 
          if(u!=head)
            L2hk+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/
        }
        /*Update the changestat for the h->k edge*/
        for(j = 0; j < nd; j++){
          deg = (Vertex)dvec[j];
          cs[j] += ((L2hk + echange == deg) - (L2hk == deg));
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*
Changescore for ESPs based on recursive two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j).

RTP:
L2th - count t<->k<->h
L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t
L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t
L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h
L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
void espRTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs) { 
  Edge e, f;
  int i, j, echange, htedge;
  int L2th,L2tk,L2kt,L2hk,L2kh; /*Two-path counts for various edges*/
  Vertex deg;
  Vertex tail, head, k, u;
  
  for(i=0;i<nd;i++)  /*Initialize the changestat vector*/
    cs[i]=0.0;

  FOR_EACH_TOGGLE(i) {
    L2th=0;
    echange = (IS_OUTEDGE(tail=TAIL(i), head=HEAD(i)) == 0) ? 1 : -1;
    htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/
    /* step through inedges of tail (k->t: k!=h,h->t,k<->h)*/
    STEP_THROUGH_INEDGES(tail,e,k){
      if(k!=head){
        /*Do we have a t<->k<->h TP?  If so, add it to our count.*/
        L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head));
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/
          L2kt=0;
          /*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=tail)&&(IS_OUTEDGE(u,k)))
              L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/
          }
          /*Update the changestat for the k->t edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2kt + echange == deg) - (L2kt == deg));
          }
        }
      }
    }
    /* step through outedges of tail (t->k: k!=h,h->t,k<->h)*/
    STEP_THROUGH_OUTEDGES(tail,e,k) {
      if(k!=head){
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/
          L2tk=0;
          /*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=tail)&&(IS_OUTEDGE(u,k)))
              L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/
          }
          /*Update the changestat for the t->k edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2tk + echange == deg) - (L2tk == deg));
          }
        }
      }
    }
    /* step through inedges of head (k->h: k!=t,h->t,k<->t)*/
    STEP_THROUGH_INEDGES(head,e,k){
      if(k!=tail){
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/
          L2kh=0;
          /*Now, count # u such that k<->u<->h (to get k->h's ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=head)&&(IS_OUTEDGE(u,k)))
              L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/
          }
          /*Update the changestat for the k->h edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2kh + echange == deg) - (L2kh == deg));
          }
        }
      }
    }
    /* step through outedges of head (h->k: k!=t,h->t,k<->t)*/
    STEP_THROUGH_OUTEDGES(head,e,k){
      if(k!=tail){
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/
          L2hk=0;
          /*Now, count # u such that k<->u<->h (to get h->k's ESP value)*/
          STEP_THROUGH_OUTEDGES(k,f,u){ 
            if((u!=head)&&(IS_OUTEDGE(u,k)))
              L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/
          }
          /*Update the changestat for the h->k edge*/
          for(j = 0; j < nd; j++){
            deg = (Vertex)dvec[j];
            cs[j] += ((L2hk + echange == deg) - (L2hk == deg));
          }
        }
      }
    }
    /*Finally, update the changestat for the t->h edge*/
    for(j = 0; j < nd; j++){
      deg = (Vertex)dvec[j];
      cs[j] += echange*(L2th == deg);
    }
    TOGGLE_IF_MORE_TO_COME(i);
  }  
  UNDO_PREVIOUS_TOGGLES(i);
}


/*****************
 changestat: d_esp
*****************/
/*
Note that d_esp is a meta-function, dispatching actual changescore 
calculation to one of the esp*_calc routines, based on the selected shared 
partner type code.

Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Recursive two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
D_CHANGESTAT_FN(d_desp) { 
  int type;
  double *dvec,*cs;
  
  /*Set things up*/
  ZERO_ALL_CHANGESTATS(i);
  type=(int)INPUT_PARAM[0];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/
  cs=CHANGE_STAT;               /*Grab the pointer to the CS vector*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: espUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPOTP: espOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPITP: espITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPRTP: espRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPOSP: espOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
    case ESPISP: espISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


/*****************
 changestat: d_gwesp
*****************/

/*
Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Recursive two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
D_CHANGESTAT_FN(d_dgwesp) { 
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[0];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[1];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[2];   /*Get the max ESP cutoff to use*/
  cs=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/
  dvec=Calloc(maxesp,double);   /*Allocate memory for the ESP vals*/
  for(i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1.0;

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: espUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPOTP: espOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPITP: espITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPRTP: espRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPOSP: espOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
    case ESPISP: espISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs); break;
  }
  
  /*Compute the gwesp changescore*/
  for(i=0;i<maxesp;i++){
    if(cs[i]!=0.0)  {                /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*cs[i];
      //Rprintf("count %f: %f, ChangeStat %f\n", dvec[i], cs[i], CHANGE_STAT[0]);
    }
  }
  CHANGE_STAT[0]*=exp(alpha);

  Free(cs);
  Free(dvec);
}



/*****************
 changestat: d_nsp
*****************/
/*
Note that d_esp is a meta-function, dispatching actual changescore 
calculation to one of the esp*_calc routines, based on the selected shared 
partner type code.

Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Recursive two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
D_CHANGESTAT_FN(d_dnsp) { 
  int i,type;
  double *dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  ZERO_ALL_CHANGESTATS(i);
  type=(int)INPUT_PARAM[0];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/
  cs_esp=Calloc(N_CHANGE_STATS,double);     /*Allocate memory for the DSP changescores*/
  cs_dsp=Calloc(N_CHANGE_STATS,double);     /*Allocate memory for the DSP changescores*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: 
    espUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspUTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPOTP: 
    espOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspOTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPITP: 
    espITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspITP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPRTP: 
    espRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspRTP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPOSP: 
    espOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspOSP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPISP: 
    espISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_esp);
    dspISP_calc(ntoggles,tails,heads,mtp,nwp,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
  
  for(i=0;i<N_CHANGE_STATS;i++)
    CHANGE_STAT[i]=(cs_dsp[i]-cs_esp[i]);

  Free(cs_esp);
  Free(cs_dsp);
}


/*****************
 changestat: d_gwnsp
*****************/

/*
Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Recursive two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
D_CHANGESTAT_FN(d_dgwnsp) { 
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[0];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[1];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[2];   /*Get the max ESP cutoff to use*/
  cs_esp=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/
  dvec=Calloc(maxesp,double);   /*Allocate memory for the ESP vals*/
  for(i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1.0;
  
  cs_dsp=Calloc(maxesp,double);     /*Allocate memory for the ESP changescores*/

  /*Obtain the changescores (by type)*/
  switch(type){
    case ESPUTP: 
      espUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspUTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPOTP: 
      espOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspOTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPITP: 
      espITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspITP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPRTP: 
      espRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspRTP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPOSP: 
      espOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspOSP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
    case ESPISP: 
      espISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_esp);
      dspISP_calc(ntoggles,tails,heads,mtp,nwp,maxesp,dvec,cs_dsp); 
      break;
  }
  
  /*Compute the gwnsp changescore*/
  for(i=0;i<maxesp;i++)
    if((cs_dsp[i]-cs_esp[i])!=0.0)                  /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*(cs_dsp[i]-cs_esp[i]);
  CHANGE_STAT[0]*=exp(alpha);

  Free(cs_esp);
  Free(cs_dsp);
  Free(dvec);  
}



/*****************
 changestat: d_ddspbwrap
*****************/

D_CHANGESTAT_FN(d_ddspbwrap) {
  d_ddsp(ntoggles, tails, heads, mtp, nwp);
  
  // correct for double counting of directed vs. undirected dyads
  for(int ind = 0; ind < N_CHANGE_STATS; ind++) CHANGE_STAT[ind] /= 2.0;
}


/*****************
 changestat: d_dgwdspbwrap
*****************/

D_CHANGESTAT_FN(d_dgwdspbwrap) {
  d_dgwdsp(ntoggles, tails, heads, mtp, nwp);
  
  // correct for double counting of directed vs. undirected dyads
  CHANGE_STAT[0] /= 2.0;
}

