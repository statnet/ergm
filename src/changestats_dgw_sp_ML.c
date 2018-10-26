#include "changestats_dgw_sp_ML.h"
#include "ergm_dyad_hashmap.h"

#define SETUP_calc_dsp							\
  memset(cs, 0, nd*sizeof(double));					\
  SETUP_update_spcache;

#define SETUP_calc_esp							\
  SETUP_calc_dsp;							\
  l3fc = ergm_LayerLogic2(t0, h0, tail, head, ll3, TRUE);		\
  l3rc = DIRECTED ? ergm_LayerLogic2(h0, t0, tail, head, ll3, TRUE) : 0;

#define INC_IF_TWOPATH(ij, t1, h1, t2, h2) if(ergm_LayerLogic2Path(t1,h1,t2,h2, ll1, ll2, any_order)) L2 ## ij ++;

#define UPDATE_CS_1(ij, t1, h1, t2, h2, EXTRA)				\
  {									\
  int c2path = ergm_c_LayerLogic2Path(t1,h1,t2,h2,			\
				      ll1,ll2, any_order,		\
				      l1c,l2c,0,0);			\
  for(unsigned int j = 0; j < nd; j++){					\
    Vertex deg = dvec[j];						\
    cs[j] += ((L2 ## ij	+ c2path == deg)				\
	      - (L2 ## ij == deg)) EXTRA;				\
  }									\
  }

#define UPDATE_CS_2(ij, t1, h1, t2, h2, EXTRA)				\
  {									\
  int c2path = ergm_c_LayerLogic2Path(t1,h1,t2,h2,			\
				      ll1,ll2, any_order,		\
				      0,0,l1c,l2c);			\
  for(unsigned int j = 0; j < nd; j++){					\
    Vertex deg = dvec[j];						\
    cs[j] += ((L2 ## ij + c2path == deg)				\
	      - (L2 ## ij == deg)) EXTRA;				\
  }									\
  }


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
static inline void dspUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_dsp;
  any_order = TRUE;
  
  CALC_with_dirs({
  /* step through edges of head */
  ML_EXEC_THROUGH_EDGES(ll0, h,e,u, {
      if (u!=t){
	unsigned int L2tu = 0;
	if(spcache) L2tu = GETDMUI(t,u,spcache);
	else{
	  /* step through edges of u */
	  ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {
	      /* Confirm that 2-path u - v - t satisfies layer conditions */
	      INC_IF_TWOPATH(tu,t,v,v,u);
	    });
	}
	UPDATE_CS_1(tu,t,h,h,u,);
      }
    });
  ML_EXEC_THROUGH_EDGES(ll0, t,e,u, {
      if (u!=h){
        unsigned int L2uh = 0;
	if(spcache) L2uh = GETDMUI(u,h,spcache);
	else{
	  /* step through edges of u */
	  ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {
	      /* Confirm that 2-path u - v - h satisfies layer conditions */
	      INC_IF_TWOPATH(uh,u,v,v,h);
	    });
	}
	UPDATE_CS_2(uh,u,t,t,h,);
      }
    });
    });
}



/*
Changescore for dsps based on outgoing two-paths, i.e. configurations for non-edge i->j such that i->k->j.

This function should only be used in the directed case
*/
static inline void dspOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_dsp;

  CALC_with_dirs({
  /* step through outedges of head (i.e., k: t->k)*/
  ML_EXEC_THROUGH_OUTEDGES(ll0,h, e, k, {
      if(k!=t){ /*Only use contingent cases*/
        unsigned int L2tk = 0;
	if(spcache) L2tk = GETDMUI(t,k,spcache);
	else{
	  /* step through inedges of k, incl. (h,k) itself */
	  ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {
	      INC_IF_TWOPATH(tk,t,u,u,k);
	    });
	}
        /*Update the changestat for the t->k edge*/
	UPDATE_CS_1(tk,t,h,h,k,);
      }
    });
    /* step through inedges of tail (i.e., k: k->t)*/
    ML_EXEC_THROUGH_INEDGES(ll0,t, e, k, {
      if (k!=h){ /*Only use contingent cases*/
        unsigned int L2kh = 0;
	if(spcache) L2kh = GETDMUI(k,h,spcache);
	else{
	  /* step through outedges of k , incl. (k,tail) itself */
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {
	      INC_IF_TWOPATH(kh,k,u,u,h);
	    });
	}
        /*Update the changestat for the k->t edge*/
	UPDATE_CS_2(kh,k,t,t,h,);
      }
    });
    });
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
static inline void dspITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_dsp;

  CALC_with_dirs({
  /* step through outedges of head (i.e., k: h->k)*/
  ML_EXEC_THROUGH_OUTEDGES(ll0, h, e, k, {
      if((k!=t)){ /*Only use contingent cases*/
        unsigned int L2kt = 0;
	if(spcache) L2kt = GETDMUI(t,k,spcache); // spcache is an OTP cache.
	else{
	  ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {
	      INC_IF_TWOPATH(kt,t,u,u,k);
	    });
	}
        /*Update the changestat for the h->k edge*/
	UPDATE_CS_1(kt,t,h,h,k,);
      }
    });
    /* step through inedges of tail (i.e., k: k->t)*/
    ML_EXEC_THROUGH_INEDGES(ll0,t, e, k, {
      if((k!=h)){ /*Only use contingent cases*/
        unsigned int L2hk = 0;
	if(spcache) L2hk = GETDMUI(k,h,spcache); // spcache is an OTP cache.
	else{
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {
	      INC_IF_TWOPATH(hk,k,u,u,h);
	    });
	}
        /*Update the changestat for the k->t edge*/
	UPDATE_CS_2(hk,k,t,t,h,);
      }
    });
    });
}


/*
Changescore for DSPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

OSP:
L2th - count t->k, h->k
L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
static inline void dspOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_dsp;
  any_order = TRUE;

  CALC_with_dirs({
  ML_EXEC_THROUGH_INEDGES(ll0,h, e, k, {
      if(k!=t){
        unsigned int L2tk = 0;
	if(spcache) L2tk = GETDMUI(t,k,spcache);
	else{
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k, f, u, {
	      if (u!=t)
		/*Increment if there is an OSP  */
		INC_IF_TWOPATH(tk,t,u,k,u);
	    });
	}
        /*Update the changestat for the t->k edge*/
	// twice (one for each direction of t<->k)
	UPDATE_CS_1(tk,t,h,k,h,*2);
      }
    });
    });
}


/*
Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

ISP:
L2th - count k->t, k->h
L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
static inline void dspISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_dsp;
  any_order = TRUE;

  CALC_with_dirs({
  ML_EXEC_THROUGH_OUTEDGES(ll0,t, e, k, {
      if(k!=h){
        unsigned int L2kh = 0;
	if(spcache) L2kh = GETDMUI(k,h,spcache);
	else{
	  ML_EXEC_THROUGH_INEDGES(ll0,k, f, u, {
	      if(u!=h)
		/*Increment if there is an ISP*/
		INC_IF_TWOPATH(kh,u,k,u,h);
	    });
	}
        /*Update the changestat for the k->h edge*/
	// twice (one for each direction of h<->k)
	UPDATE_CS_1(kh,t,h,t,k,*2);
      }
    });
    });
}


/* /\* */
/* Changescore for ESPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j). */

/* RTP: */
/* L2th - count t<->k<->h */
/* L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h */
/* L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h */

/* We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function. */
/* *\/ */
/* static inline void dspRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs) {  */
/*   Edge e, f; */
/*   int j, echange, htedge; */
/*   int L2th,L2tk,L2kt,L2hk,L2kh; /\*Two-path counts for various edges*\/ */
/*   Vertex deg; */
/*   Vertex k, u; */
  
/*   memset(cs, 0, nd*sizeof(double)); */
/*     L2th=0; */
/*     echange = (IS_OUTEDGE(tail,head) == 0) ? 1 : -1; */
/*     htedge=IS_OUTEDGE(head,tail);  /\*Is there an h->t (reciprocating) edge?*\/ */
/*     /\* step through inedges of tail (k->t: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         /\*Do we have a t<->k<->h TP?  If so, add it to our count.*\/ */
/*         L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)); */
/*         if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /\*Only consider stats that could change*\/ */
/*           L2kt=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->t edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kt + echange == deg) - (L2kt == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of tail (t->k: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /\*Only consider stats that could change*\/ */
/*           L2tk=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the t->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2tk + echange == deg) - (L2tk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through inedges of head (k->h: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /\*Only consider stats that could change*\/ */
/*           L2kh=0; */
/*           /\*Now, count # u such that k<->u<->h (to get k->h's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->h edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kh + echange == deg) - (L2kh == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of head (h->k: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /\*Only consider stats that could change*\/ */
/*           L2hk=0; */
/*           /\*Now, count # u such that k<->u<->h (to get h->k's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the h->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2hk + echange == deg) - (L2hk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\*Finally, update the changestat for the t->h edge*\/ */
/*     for(unsigned int j = 0; j < nd; j++){ */
/*       Vertex deg = dvec[j]; */
/*       cs[j] += echange*(L2th == deg); */
/*     } */
/* } */


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
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
C_CHANGESTAT_FN(c_ddsp_ML) { 
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[3]>=0) ? AUX_STORAGE_NUM(3) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[4];
  
  /*Set things up*/
  unsigned int type=(int)INPUT_PARAM[5];     /*Get the ESP type code to be used*/
  double *dvec=INPUT_PARAM+6;           /*Get the pointer to the ESP stats list*/
  double *cs=CHANGE_STAT;               /*Grab the pointer to the CS vector*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP: dspUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break;
  case ESPOTP: dspOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break;
  case ESPITP: dspITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break;
  /* case ESPRTP: dspRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break; */
  case ESPOSP: dspOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break;
  case ESPISP: dspISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs); break;
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
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/

I_CHANGESTAT_FN(i_dgwdsp_ML) {
  Vertex maxesp = (unsigned int)INPUT_PARAM[7]; // Index must be same as for maxesp below.
  ALLOC_STORAGE(maxesp*2, double, storage);
  double *dvec=storage+maxesp;          /*Grab memory for the ESP vals*/
  for(unsigned int i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1;
}

C_CHANGESTAT_FN(c_dgwdsp_ML) {
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[3]>=0) ? AUX_STORAGE_NUM(3) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[4];
  GET_STORAGE(double, storage);
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[5];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[6];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[7];   /*Get the max ESP cutoff to use*/
  cs=storage;                   /*Grab memory for the ESP changescores*/
  dvec=storage+maxesp;          /*Grab memory for the ESP vals*/

  /*Obtain the DSP changescores (by type)*/
  switch(type){
    case ESPUTP: dspUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break;
    case ESPOTP: dspOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break;
    case ESPITP: dspITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break;
    /* case ESPRTP: dspRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break; */
    case ESPOSP: dspOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break;
    case ESPISP: dspISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs); break;
  }
  
  /*Compute the gwdsp changescore*/
  for(i=0;i<maxesp;i++){
    if(cs[i]!=0.0) {                 /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*cs[i];
      //Rprintf("count %f: %f\n", dvec[i], cs[i]);
    }
  }
  CHANGE_STAT[0]*=exp(alpha);
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
static inline void espUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_esp;
  any_order = TRUE;
  
  CALC_with_dirs({
  unsigned int L2th = 0;
  if(spcache) L2th = GETDMUI(t,h,spcache);
  ML_EXEC_THROUGH_EDGES(ll0, h,e,u, {
      if (ML_IS_UNDIRECTED_EDGE(ll0,u,t) != 0){
        if(!spcache && l3c) INC_IF_TWOPATH(th,t,u,u,h);
	unsigned int Buh = ML_GETWT(ll3,u,h);
	unsigned int Btu = ML_GETWT(ll3,t,u);
        unsigned int L2tu = 0;
        unsigned int L2uh = 0;
	if(spcache){
	  L2tu = GETDMUI(t,u,spcache);
	  L2uh = GETDMUI(u,h,spcache);
	}else{
	  ML_EXEC_THROUGH_EDGES(ll0, u,f,v, {
	      if(Buh) INC_IF_TWOPATH(uh,u,v,v,h);
	      if(Btu) INC_IF_TWOPATH(tu,t,v,v,u);
	    });
	}
	if(Buh) UPDATE_CS_2(uh,u,t,t,h,);
	if(Btu) UPDATE_CS_1(tu,t,h,h,u,);
      }
    });
  if(l3c)
    for(unsigned int j = 0; j < nd; j++){
      Vertex deg = dvec[j];
      cs[j] += l3c*(L2th == deg);
    }
    });
}


/*
Changescore for ESPs based on outgoing two-paths, i.e. configurations for edge i->j such that i->k->j.

OTP:
L2th - count i->k->j
L2tk - for each i->k neq j: j->k, count u such that i->u->k
L2kh - for each k->j neq i: k->i, count u such that k->u->j

This function should only be used in the directed case, with espUTP being used in the undirected case.
*/
static inline void espOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_esp;

  CALC_with_dirs({
  if(l3c){
    unsigned int L2th = 0;
    if(spcache) L2th = GETDMUI(t,h,spcache);
    else{
      ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, {
	  if(k!=h)
	    INC_IF_TWOPATH(th,t,k,k,h);
	});
    }
    for(unsigned int j = 0; j < nd; j++){
      Vertex deg = dvec[j];
      cs[j] += l3c*(L2th == deg);
    }
  }
  ML_EXEC_THROUGH_OUTEDGES(ll3,t,e,k, {
      if((k!=h)&&(ML_IS_OUTEDGE(ll0,h,k))){ /*Only use contingent cases*/
        unsigned int L2tk = 0;
	if(spcache) L2tk = GETDMUI(t,k,spcache);
	else{
	  ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {
	      //Rprintf("\t\tTrying u==%d\n",u);
	      if(u!=t) 
		INC_IF_TWOPATH(tk,t,u,u,k);
	    });
	}
	UPDATE_CS_1(tk,t,h,h,k,);
      }
    });
  /* step through inedges of h (i.e., k: k->h)*/
  ML_EXEC_THROUGH_INEDGES(ll3,h,e,k, {
      if((k!=t)&&(ML_IS_OUTEDGE(ll0,k,t))){ /*Only use contingent cases*/
        unsigned int L2kh = 0;
	if(spcache) L2kh = GETDMUI(k,h,spcache);
	else{
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {
	      if(u!=h) 
		INC_IF_TWOPATH(kh,k,u,u,h);
	    });
	}
	UPDATE_CS_2(kh,k,t,t,h,);
      }
    });
    });
}


/*
Changescore for ESPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.

ITP:
L2th - count j->k->i
L2hk - for each j->k neq i: k->i, count u such that k->u->j 
L2kt - for each k->i neq j: j->k, count u such that i->u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
static inline void espITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) { 
   SETUP_calc_esp;
   
  CALC_with_dirs({
  if(l3c){
    unsigned int L2th = 0;
    if(spcache) L2th = GETDMUI(h,t,spcache); // spcache is an OTP cache.
    else{
      ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, {
	  if(k!=h)
	    INC_IF_TWOPATH(th,h,k,k,t);
	});
    }
    for(unsigned int j = 0; j < nd; j++){
      Vertex deg = dvec[j];
      cs[j] += l3c*(L2th == deg);
    }
  }

  ML_EXEC_THROUGH_OUTEDGES(ll3,h,e,k, {
      if((k!=t)&&(ML_IS_OUTEDGE(ll0,k,t))){ /*Only use contingent cases*/
        //Rprintf("\tk==%d, passed criteria\n",k);
        /*We have a h->k->t two-path, so add it to our count.*/
        unsigned int L2hk = 0;
	if(spcache) L2hk = GETDMUI(k,h,spcache); // spcache is an OTP cache.
	else{
        /*Now, count # u such that k->u->h (so that we know k's ESP value)*/
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {
	      //Rprintf("\t\tTrying u==%d\n",u);
	      if(u!=h) 
		INC_IF_TWOPATH(hk,k,u,u,h);
	    });
	}
        /*Update the changestat for the h->k edge*/
        //Rprintf("\t2-path count was %d\n",L2hk);
	UPDATE_CS_2(hk,k,t,t,h,);
      }
    });
    /* step through inedges of t (i.e., k: k->t)*/
    ML_EXEC_THROUGH_INEDGES(ll3,t,e,k, {
	if((k!=h)&&(ML_IS_OUTEDGE(ll0,h,k))){ /*Only use contingent cases*/
	  unsigned int L2kt = 0;
	  if(spcache) L2kt = GETDMUI(t,k,spcache); // spcache is an OTP cache.
	  else{
	    /*Now, count # u such that t->u->k (so that we know k's ESP value)*/
	    ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, {
		if(u!=t) 
		  INC_IF_TWOPATH(kt,t,u,u,k);
	      });
	  }
	  /*Update the changestat for the k->t edge*/
	  //Rprintf("\t2-path count was %d\n",L2kt);
          UPDATE_CS_1(kt,t,h,h,k,);
	}
    });
    });
}


/*
Changescore for ESPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

OSP:
L2th - count t->k, h->k
L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
static inline void espOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) { 
  SETUP_calc_esp;
  
  CALC_with_dirs({
  if(l3c){
    unsigned int L2th = 0;
    if(spcache) L2th = GETDMUI(t,h,spcache);
    else{
      ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, {
	  if(k!=h)
	    INC_IF_TWOPATH(th,t,k,h,k);
	});
    }
    for(unsigned int j = 0; j < nd; j++){
      Vertex deg = dvec[j];
      cs[j] += l3c*(L2th == deg);
    }
  }

  /* step through outedges of t (i.e., k: t->k, k->h, k!=h)*/
  ML_EXEC_THROUGH_OUTEDGES(ll3,t,e,k, {
      if(k!=h && ML_IS_OUTEDGE(ll0,k,h)){ /*Only consider stats that could change*/
	unsigned int L2tk = 0;
	if(spcache) L2tk = GETDMUI(t,k,spcache);
	else{
	  /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, { 
	      if(u!=t)
		/*Increment if there is an OSP*/
		INC_IF_TWOPATH(tk,t,u,k,u);
	    });
	}
	/*Update the changestat for the t->k edge*/
	UPDATE_CS_1(tk,t,h,k,h,);
      }
    });
  /* step through inedges of t (i.e., k: k->t, k->h, k!=h)*/
  ML_EXEC_THROUGH_INEDGES(ll3,t,e,k, {
      if(k!=h && ML_IS_OUTEDGE(ll0,k,h)){ /*Only stats that could change*/
	unsigned int L2kt=0;
	if(spcache) L2kt = GETDMUI(k,t,spcache);
	else{
	  /*Now, count # u such that t->u,k->u (to get k->t's ESP value)*/
	  ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, { 
	      if(u!=t)
		/*Increment if there is an OSP*/
		INC_IF_TWOPATH(kt,t,u,k,u);
	    });
	}
        /*Update the changestat for the k->t edge*/
	UPDATE_CS_1(kt,t,h,k,h,);
      }
    });
    });
}


/*
Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

ISP:
L2th - count k->t, k->h
L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
static inline void espISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) { 
   SETUP_calc_esp;
  
  CALC_with_dirs({
  if(l3c){
    unsigned int L2th = 0;
    if(spcache) L2th = GETDMUI(t,h,spcache);
    else{
      ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, {
	  if(k!=h)
	    INC_IF_TWOPATH(th,k,t,k,h);
	});
    }
    for(unsigned int j = 0; j < nd; j++){
      Vertex deg = dvec[j];
      cs[j] += l3c*(L2th == deg);
    }
  }
  
  /* step through inedges of h (i.e., k: k->h, t->k, k!=t)*/
  ML_EXEC_THROUGH_INEDGES(ll3,h,e,k, {
      if(k!=t && ML_IS_OUTEDGE(ll0,t,k)){
	unsigned int L2kh = 0;
	if(spcache) L2kh = GETDMUI(k,h,spcache);
	else{
	  /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/
	  ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, { 
	      if(u!=h)
		/*Increment if there is an ISP*/
		INC_IF_TWOPATH(kh,u,k,u,h);
	    });
	}
	/*Update the changestat for the k->h edge*/
	UPDATE_CS_1(kh,t,h,t,k,);
      }
    });
    /* step through outedges of h (i.e., k: h->k, t->k, k!=t)*/
  ML_EXEC_THROUGH_OUTEDGES(ll3,h,e,k, {
      if(k!=t && ML_IS_OUTEDGE(ll0,t,k)){ /*Only stats that could change*/
      unsigned int L2hk = 0;
      if(spcache) L2hk = GETDMUI(h,k,spcache);
      else{
	/*Now, count # u such that u->h,u->k (to get k->h's ESP value)*/
	ML_EXEC_THROUGH_INEDGES(ll0,k,f,u, { 
	    if(u!=h)
	      /*Increment if there is an ISP*/
	      INC_IF_TWOPATH(hk,u,k,u,h);
	  });
      }
      /*Update the changestat for the h->k edge*/
      UPDATE_CS_1(hk,t,h,t,k,);
      }
    });
    });
}


/* /\* */
/* Changescore for ESPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j). */

/* RTP: */
/* L2th - count t<->k<->h */
/* L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t */
/* L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h */
/* L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h */

/* We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function. */
/* *\/ */
/* static inline void espRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreDyadMapUInt *spcache, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs) {  */
/*   Edge e, f; */
/*   int j, echange, htedge; */
/*   int L2th,L2tk,L2kt,L2hk,L2kh; /\*Two-path counts for various edges*\/ */
/*   Vertex deg; */
/*   Vertex k, u; */
  
/*   memset(cs, 0, nd*sizeof(double)); */
 
/*     L2th=0; */
/*     echange = (IS_OUTEDGE(t,h) == 0) ? 1 : -1; */
/*     htedge=IS_OUTEDGE(h,t);  /\*Is there an h->t (reciprocating) edge?*\/ */
/*     /\* step through inedges of t (k->t: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         /\*Do we have a t<->k<->h TP?  If so, add it to our count.*\/ */
/*         L2th+=(IS_OUTEDGE(t,k)&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)); */
/*         if(htedge&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)){ /\*Only consider stats that could change*\/ */
/*           L2kt=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2kt+=(IS_OUTEDGE(u,t)&&IS_OUTEDGE(t,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->t edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kt + echange == deg) - (L2kt == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of t (t->k: k!=h,h->t,k<->h)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,t,e,k, { */
/*       if(k!=h){ */
/*         if(htedge&&IS_OUTEDGE(h,k)&&IS_OUTEDGE(k,h)){ /\*Only consider stats that could change*\/ */
/*           L2tk=0; */
/*           /\*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=t)&&(IS_OUTEDGE(u,k))) */
/*               L2tk+=(IS_OUTEDGE(u,t)&&IS_OUTEDGE(t,u));  /\*k<->u<->t?*\/ */
/*           }); */
/*           /\*Update the changestat for the t->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2tk + echange == deg) - (L2tk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through inedges of h (k->h: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_INEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(t,k)&&IS_OUTEDGE(k,t)){ /\*Only consider stats that could change*\/ */
/*           L2kh=0; */
/*           /\*Now, count # u such that k<->u<->h (to get k->h's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2kh+=(IS_OUTEDGE(u,h)&&IS_OUTEDGE(h,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the k->h edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2kh + echange == deg) - (L2kh == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\* step through outedges of h (h->k: k!=t,h->t,k<->t)*\/ */
/*     ML_EXEC_THROUGH_OUTEDGES(ll0,h,e,k, { */
/*       if(k!=t){ */
/*         if(htedge&&IS_OUTEDGE(t,k)&&IS_OUTEDGE(k,t)){ /\*Only consider stats that could change*\/ */
/*           L2hk=0; */
/*           /\*Now, count # u such that k<->u<->h (to get h->k's ESP value)*\/ */
/*           ML_EXEC_THROUGH_OUTEDGES(ll0,k,f,u, {  */
/*             if((u!=h)&&(IS_OUTEDGE(u,k))) */
/*               L2hk+=(IS_OUTEDGE(u,h)&&IS_OUTEDGE(h,u));  /\*k<->u<->h?*\/ */
/*           }); */
/*           /\*Update the changestat for the h->k edge*\/ */
/*           for(unsigned int j = 0; j < nd; j++){ */
/*             Vertex deg = dvec[j]; */
/*             cs[j] += ((L2hk + echange == deg) - (L2hk == deg)); */
/*           } */
/*         } */
/*       } */
/*     }); */
/*     /\*Finally, update the changestat for the t->h edge*\/ */
/*     for(unsigned int j = 0; j < nd; j++){ */
/*       Vertex deg = dvec[j]; */
/*       cs[j] += echange*(L2th == deg); */
/*     } */
/* } */


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
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
C_CHANGESTAT_FN(c_desp_ML) {
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[4]>=0) ? AUX_STORAGE_NUM(4) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[5];

  int type;
  double *dvec,*cs;
  
  /*Set things up*/
  ZERO_ALL_CHANGESTATS(i);
  type=(int)INPUT_PARAM[6];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+7;           /*Get the pointer to the ESP stats list*/
  cs=CHANGE_STAT;               /*Grab the pointer to the CS vector*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: espUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break;
    case ESPOTP: espOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break;
    case ESPITP: espITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break;
    /* case ESPRTP: espRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break; */
    case ESPOSP: espOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break;
    case ESPISP: espISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs); break;
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
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
I_CHANGESTAT_FN(i_dgwesp_ML) {
  Vertex maxesp = (unsigned int)INPUT_PARAM[8]; // Index must be same as for maxesp below.
  ALLOC_STORAGE(maxesp*2, double, storage);
  double *dvec=storage+maxesp;          /*Grab memory for the ESP vals*/
  for(unsigned int i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1;
}

C_CHANGESTAT_FN(c_dgwesp_ML) { 
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[4]>=0) ? AUX_STORAGE_NUM(4) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[5];
  GET_STORAGE(double, storage);
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[6];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[7];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[8];   /*Get the max ESP cutoff to use*/
  cs=storage;                   /*Grab memory for the ESP changescores*/
  dvec=storage+maxesp;          /*Grab memory for the ESP vals*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
    case ESPUTP: espUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break;
    case ESPOTP: espOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break;
    case ESPITP: espITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break;
    /* case ESPRTP: espRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break; */
    case ESPOSP: espOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break;
    case ESPISP: espISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs); break;
  }
  
  /*Compute the gwesp changescore*/
  for(i=0;i<maxesp;i++){
    if(cs[i]!=0.0)  {                /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*cs[i];
      //Rprintf("count %f: %f, ChangeStat %f\n", dvec[i], cs[i], CHANGE_STAT[0]);
    }
  }
  CHANGE_STAT[0]*=exp(alpha);
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
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
I_CHANGESTAT_FN(i_dnsp_ML) {
  ALLOC_STORAGE((int)N_CHANGE_STATS*2, double, storage);
  (void)storage; // Get rid of an unused warning.
}

C_CHANGESTAT_FN(c_dnsp_ML) {
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[4]>=0) ? AUX_STORAGE_NUM(4) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[5];
  GET_STORAGE(double, storage);
  int i,type;
  double *dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  type=(int)INPUT_PARAM[6];     /*Get the ESP type code to be used*/
  dvec=INPUT_PARAM+7;           /*Get the pointer to the ESP stats list*/
  cs_esp=storage;               /*Grab memory for the DSP changescores*/
  cs_dsp=storage+N_CHANGE_STATS;/*Grab memory for the DSP changescores*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: 
    espUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp);
    dspUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPOTP: 
    espOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp);
    dspOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPITP: 
    espITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp);
    dspITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  /* case ESPRTP:  */
  /*   espRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp); */
  /*   dspRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp);  */
  /*   break; */
  case ESPOSP: 
    espOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp);
    dspOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  case ESPISP: 
    espISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,N_CHANGE_STATS,dvec,cs_esp);
    dspISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,N_CHANGE_STATS,dvec,cs_dsp); 
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
  
  for(i=0;i<N_CHANGE_STATS;i++)
    CHANGE_STAT[i]=(cs_dsp[i]-cs_esp[i]);
}


/*****************
 changestat: d_gwnsp
*****************/

/*
Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
I_CHANGESTAT_FN(i_dgwnsp_ML) {
  Vertex maxesp = (unsigned int)INPUT_PARAM[8]; // Index must be same as for maxesp below.
  ALLOC_STORAGE(maxesp*3, double, storage);
  double *dvec=storage+maxesp;          /*Grab memory for the ESP vals*/
  for(unsigned int i=0;i<maxesp;i++)         /*Initialize the ESP vals*/
    dvec[i]=i+1;
}

C_CHANGESTAT_FN(c_dgwnsp_ML) {
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll0, 0);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll1, 1);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll2, 2);
  GET_AUX_STORAGE_NUM(StoreLayerLogic, ll3, 3);
  StoreDyadMapUInt *spcache = (INPUT_PARAM[4]>=0) ? AUX_STORAGE_NUM(4) : NULL;
  unsigned int any_order = (unsigned int) INPUT_PARAM[5];
  GET_STORAGE(double, storage);
  int type;
  Vertex i,maxesp;
  double alpha, oneexpa,*dvec,*cs_esp, *cs_dsp;
  
  /*Set things up*/
  CHANGE_STAT[0] = 0.0;         /*Zero the changestat*/
  alpha = INPUT_PARAM[6];       /*Get alpha*/
  oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  type=(int)INPUT_PARAM[7];     /*Get the ESP type code to be used*/
  maxesp=(int)INPUT_PARAM[8];   /*Get the max ESP cutoff to use*/
  cs_esp=storage;     /*Grab memory for the ESP changescores*/
  dvec=storage+maxesp;   /*Grab memory for the ESP vals*/
  cs_dsp=storage+maxesp+maxesp;     /*Grab memory for the ESP changescores*/

  /*Obtain the changescores (by type)*/
  switch(type){
    case ESPUTP: 
      espUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp);
      dspUTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp); 
      break;
    case ESPOTP: 
      espOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp);
      dspOTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp); 
      break;
    case ESPITP: 
      espITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp);
      dspITP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp); 
      break;
    /* case ESPRTP:  */
    /*   espRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp); */
    /*   dspRTP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp);  */
    /*   break; */
    case ESPOSP: 
      espOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp);
      dspOSP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp); 
      break;
    case ESPISP: 
      espISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,ll3,any_order,maxesp,dvec,cs_esp);
      dspISP_ML_calc(tail,head,mtp,nwp,spcache,ll0,ll1,ll2,any_order,maxesp,dvec,cs_dsp); 
      break;
  }
  
  /*Compute the gwnsp changescore*/
  for(i=0;i<maxesp;i++)
    if((cs_dsp[i]-cs_esp[i])!=0.0)                  /*Filtering to save a few pow() calls*/
      CHANGE_STAT[0]+=(1.0-pow(oneexpa,dvec[i]))*(cs_dsp[i]-cs_esp[i]);
  CHANGE_STAT[0]*=exp(alpha);
}



