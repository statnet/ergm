/*  File src/changestats_dgw_sp.c in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#include "changestats_dgw_sp.h"
#include "ergm_dyad_hashmap.h"
#include "changestats.h"

#define all_calcs(term)                         \
  dvec_calc(term)                               \
       dist_calc(term)                          \
       gw_calc(term)

#define all_calcs2(term)                        \
  dvec_calc2(term)                              \
       dist_calc2(term)                         \
       gw_calc2(term)


#define sp_args tail,head,mtp,nwp,edgestate,spcache,N_CHANGE_STATS,dvec,CHANGE_STAT

#define dvec_calc(term)                                                 \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, int nd, Vertex *dvec, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change(L, {                                                \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += ((L+echange == deg) - (L == deg));                   \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += (echange)*(L == deg);                                \
        }                                                               \
      });                                                               \
  }

#define dvec_calc2(term)                                                \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, int nd, Vertex *dvec, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change(L, {                                                \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += ((L+echange == deg) - (L == deg))*2;                 \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += (echange)*(L == deg);                                \
        }                                                               \
      });                                                               \
  }


#define spd_args tail,head,mtp,nwp,edgestate,spcache,N_CHANGE_STATS,CHANGE_STAT

#define dist_calc(term)                                                 \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, int nd, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change(L, {                                                \
        int nL = L + echange;                                           \
        if(nL > nd) cutoff_error(mtp);                                  \
        if(L) cs[L-1]--;                                                \
        if(nL) cs[nL-1]++;                                              \
      },{                                                               \
        if(L > nd) cutoff_error(mtp);                                   \
        if(L) cs[L-1] += echange;                                       \
      });                                                               \
  }

#define dist_calc2(term)                                                \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, int nd, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change(L, {                                                \
        int nL = L + echange;                                           \
        if(nL > nd) cutoff_error(mtp);                                  \
        if(L) cs[L-1]-=2;                                               \
        if(nL) cs[nL-1]+=2;                                             \
      },{                                                               \
        if(L > nd) cutoff_error(mtp);                                   \
        if(L) cs[L-1] += echange;                                       \
      });                                                               \
  }


#define gwsp_args tail,head,mtp,nwp,edgestate,spcache,alpha,oneexpa

#define gw_calc(term)                                                   \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, double alpha, double oneexpa) { \
    double cumchange = 0;                                               \
    term ## _change(L, {                                                \
        cumchange += pow(oneexpa, L-edgestate);                         \
      },{                                                               \
        if(alpha < 100.0) cumchange += exp(alpha)*(1-pow(oneexpa, L));  \
        else cumchange += L;                                            \
      });                                                               \
    return cumchange;                                                   \
  }


#define gw_calc2(term)                                                  \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreDyadMapUInt *spcache, double alpha, double oneexpa) { \
    double cumchange = 0;                                               \
    term ## _change(L, {                                                \
        cumchange += pow(oneexpa, L-edgestate)*2;                       \
      },{                                                               \
        if(alpha < 100.0) cumchange += exp(alpha)*(1-pow(oneexpa, L));  \
        else cumchange += L;                                            \
      });                                                               \
    return cumchange;                                                   \
  }



#define call_subroutine_path(count, subroutine_path)    \
  {int L = (count);                                     \
    {subroutine_path}}

#define call_subroutine_focus(count, subroutine_focus)  \
  {int L = (count);                                     \
    {subroutine_focus}}


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

#define dspUTP_change(L, subroutine_path, subroutine_focus)     \
  /* step through edges of head */                              \
  EXEC_THROUGH_EDGES(head,e,u, {                                \
      if (u!=tail){                                             \
        int L2tu;                                               \
        if(spcache) L2tu = GETDMUI(tail,u,spcache);             \
        else{                                                   \
          L2tu=0;                                               \
          /* step through edges of u */                         \
          EXEC_THROUGH_EDGES(u,f,v, {                           \
              if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;        \
            });                                                 \
        }                                                       \
        call_subroutine_path(L2tu, subroutine_path);            \
      }                                                         \
    });                                                         \
  EXEC_THROUGH_EDGES(tail,e,u, {                                \
      if (u!=head){                                             \
        int L2uh;                                               \
        if(spcache) L2uh = GETDMUI(u,head,spcache);             \
        else{                                                   \
          L2uh=0;                                               \
          /* step through edges of u */                         \
          EXEC_THROUGH_EDGES(u,f,v, {                           \
              if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;        \
            });                                                 \
        }                                                       \
        call_subroutine_path(L2uh, subroutine_path);            \
      }                                                         \
    });


/*
  Changescore for dsps based on outgoing two-paths, i.e. configurations for non-edge i->j such that i->k->j.

  This function should only be used in the directed case
*/

#define dspOTP_change(L, subroutine_path, subroutine_focus)             \
  /* step through outedges of head (i.e., k: t->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head, e, k, {                                   \
      if(k!=tail){ /*Only use contingent cases*/                        \
        int L2tk;                                                       \
        if(spcache) L2tk = GETDMUI(tail,k,spcache);                     \
	else{                                                           \
	  L2tk=0;                                                       \
	  /* step through inedges of k, incl. (head,k) itself */        \
	  EXEC_THROUGH_INEDGES(k, f, u, {                               \
	      L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2tk, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if (k!=head){ /*Only use contingent cases*/                       \
        int L2kh;                                                       \
	if(spcache) L2kh = GETDMUI(k,head,spcache);                     \
	else{                                                           \
	  L2kh=0;                                                       \
	  /* step through outedges of k , incl. (k,tail) itself */      \
	  EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
	      L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2kh, subroutine_path);                    \
      }                                                                 \
    });


/*
  Changescore for DSPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.
  IE cyclical shared partners
 
  ITP:
  L2th - count j->k->i
  L2hk - for each j->k neq i: k->i, count u such that k->u->j
  L2kt - for each k->i neq j: j->k, count u such that i->u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspITP_change(L, subroutine_path, subroutine_focus)             \
  /* step through outedges of head (i.e., k: h->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head, e, k, {                                   \
      if((k!=tail)){ /*Only use contingent cases*/                      \
        int L2kt;                                                       \
        /*We have a h->k->t two-path, so add it to our count.*/         \
        if(spcache) L2kt = GETDMUI(tail,k,spcache); /* spcache is an OTP cache. */ \
        else{                                                           \
          L2kt=0;                                                       \
          /*Now, count # u such that k->u->h (so that we know k's ESP value)*/ \
          EXEC_THROUGH_INEDGES(k, f, u, {                               \
              L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(L2kt, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if((k!=head)){ /*Only use contingent cases*/                      \
        int L2hk;                                                       \
        if(spcache) L2hk = GETDMUI(k,head,spcache);                     \
        else{                                                           \
          L2hk=0;                                                       \
          /*Now, count # u such that t->u->k (so that we know k's ESP value)*/ \
          EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
              L2hk+=IS_OUTEDGE(u,head); /*Increment if there is a cyclic edge*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(L2hk, subroutine_path);                    \
      }                                                                 \
    });


/*
  Changescore for DSPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  OSP:
  L2th - count t->k, h->k
  L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
  L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspOSP_change(L, subroutine_path, subroutine_focus)             \
  /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/        \
  EXEC_THROUGH_INEDGES(head, e, k, {                                    \
      if(k!=tail){                                                      \
        int L2tk;                                                       \
        /*Do we have a t->k,h->k SP?  If so, add it to our count.*/     \
        if(spcache) L2tk = GETDMUI(tail,k,spcache);                     \
        else{                                                           \
          L2tk=0;                                                       \
          /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/ \
          EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
              if (u != tail)                                            \
                /*Increment if there is an OSP  */                      \
                L2tk+=(IS_OUTEDGE(tail,u));                             \
            });                                                         \
        }                                                               \
        call_subroutine_path(L2tk, subroutine_path);                    \
      }                                                                 \
    });


/*
  Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  ISP:
  L2th - count k->t, k->h
  L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
  L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspISP_change(L, subroutine_path, subroutine_focus)             \
  /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/         \
  EXEC_THROUGH_OUTEDGES(tail, e, k, {                                   \
      int L2kh;                                                         \
      if(k!=head){                                                      \
        if(spcache) L2kh = GETDMUI(k,head,spcache);                     \
        else{                                                           \
          L2kh=0;                                                       \
          /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/ \
          EXEC_THROUGH_INEDGES(k, f, u, {                               \
              if(u!=head)                                               \
                L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(L2kh, subroutine_path);                    \
      }                                                                 \
    });


/*
  Changescore for DSPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j).

  RTP:
  L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h
  L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t

  Thanks to the symmetries involved, this covers all cases.

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define dspRTP_change(L, subroutine_path, subroutine_focus)             \
  int htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/ \
  if(htedge){ /* Otherwise, t->h doesn't make a difference. */          \
    /* step through reciprocated outedges of tail (t->k: k!=h,k<-t)*/   \
    EXEC_THROUGH_OUTEDGES(tail,e,k,{                                    \
        if(k!=head&&IS_OUTEDGE(k,tail)){                                \
          int L2kh;                                                     \
          if(spcache) L2kh = GETDMUI(k,head,spcache);                   \
          else{                                                         \
            L2kh=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get k->h's SP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u,{                               \
                if(u!=tail&&u!=head&&(IS_OUTEDGE(u,k)))                 \
                  L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2kh, subroutine_path);                  \
        }                                                               \
      });                                                               \
    /* step through reciprocated outedges of tail (t->k: k!=h,k<-t)*/   \
    EXEC_THROUGH_OUTEDGES(head,e,k,{                                    \
        if(k!=tail&&IS_OUTEDGE(k,head)){                                \
          int L2kt;                                                     \
          if(spcache) L2kt = GETDMUI(k,tail,spcache);                   \
          else{                                                         \
            L2kt=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get k->t's SP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u,{                               \
                if(u!=head&&u!=tail&&(IS_OUTEDGE(u,k)))                 \
                  L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2kt, subroutine_path);                  \
        }                                                               \
      });                                                               \
  }

all_calcs(dspUTP)
all_calcs(dspOTP)
all_calcs(dspITP)
all_calcs2(dspOSP)
all_calcs2(dspISP)
all_calcs2(dspRTP)

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
C_CHANGESTAT_FN(c_ddsp) { 
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: dspUTP_calc(sp_args); break;
  case ESPOTP: dspOTP_calc(sp_args); break;
  case ESPITP: dspITP_calc(sp_args); break;
  case ESPRTP: dspRTP_calc(sp_args); break;
  case ESPOSP: dspOSP_calc(sp_args); break;
  case ESPISP: dspISP_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_ddspdist) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: dspUTP_dist_calc(spd_args); break;
  case ESPOTP: dspOTP_dist_calc(spd_args); break;
  case ESPITP: dspITP_dist_calc(spd_args); break;
  case ESPRTP: dspRTP_dist_calc(spd_args); break;
  case ESPOSP: dspOSP_dist_calc(spd_args); break;
  case ESPISP: dspISP_dist_calc(spd_args); break;
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
C_CHANGESTAT_FN(c_dgwdsp) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case ESPUTP: cumchange = dspUTP_gw_calc(gwsp_args); break;
  case ESPOTP: cumchange = dspOTP_gw_calc(gwsp_args); break;
  case ESPITP: cumchange = dspITP_gw_calc(gwsp_args); break;
  case ESPRTP: cumchange = dspRTP_gw_calc(gwsp_args); break;
  case ESPOSP: cumchange = dspOSP_gw_calc(gwsp_args); break;
  case ESPISP: cumchange = dspISP_gw_calc(gwsp_args); break;
  }
  
  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
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
#define espUTP_change(L, subroutine_path, subroutine_focus)     \
  int L2th;                                                     \
  if(spcache) L2th = GETDMUI(tail,head,spcache); else L2th=0;   \
  /* step through outedges of head */                           \
  EXEC_THROUGH_EDGES(head,e,u, {                                \
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){                     \
        int L2tu;                                               \
        int L2uh;                                               \
	if(spcache){                                            \
	  L2tu = GETDMUI(tail,u,spcache);                       \
	  L2uh = GETDMUI(u,head,spcache);                       \
	}else{                                                  \
	  L2th++;                                               \
	  L2tu=0;                                               \
	  L2uh=0;                                               \
	  /* step through edges of u */                         \
	  EXEC_THROUGH_EDGES(u,f,v, {                           \
	      if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;        \
	      if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;        \
	    });                                                 \
	}                                                       \
        call_subroutine_path(L2tu, subroutine_path);            \
        call_subroutine_path(L2uh, subroutine_path);            \
      }                                                         \
    });                                                         \
  call_subroutine_focus(L2th, subroutine_focus);


/*
  Changescore for ESPs based on outgoing two-paths, i.e. configurations for edge i->j such that i->k->j.

  OTP:
  L2th - count i->k->j
  L2tk - for each i->k neq j: j->k, count u such that i->u->k
  L2kh - for each k->j neq i: k->i, count u such that k->u->j

  This function should only be used in the directed case, with espUTP being used in the undirected case.
*/
#define espOTP_change(L, subroutine_path, subroutine_focus)             \
  int L2th;                                                             \
  if(spcache) L2th = GETDMUI(tail,head,spcache); else L2th=0;           \
  /* step through outedges of tail (i.e., k: t->k)*/                    \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(!spcache&&(k!=head)&&(IS_OUTEDGE(k,head))){                    \
        /*We have a t->k->h two-path, so add it to our count.*/         \
        L2th++;                                                         \
      }                                                                 \
      int L2tk;                                                         \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        if(spcache) L2tk = GETDMUI(tail,k,spcache);                     \
	else{                                                           \
	  L2tk=0;                                                       \
	  /*Now, count # u such that t->u->k (to find t->k's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=tail)                                               \
		L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2tk, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (i.e., k: k->h)*/                     \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      int L2kh;                                                         \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        if(spcache) L2kh = GETDMUI(k,head,spcache);                     \
	else{                                                           \
	  L2kh=0;                                                       \
	  /*Now, count # u such that k->u->j (to find k->h's ESP value)*/ \
	  EXEC_THROUGH_OUTEDGES(k,f,u, {                                \
	      if(u!=head)                                               \
		L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2kh, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(L2th, subroutine_focus);


/*
  Changescore for ESPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.

  ITP:
  L2th - count j->k->i
  L2hk - for each j->k neq i: k->i, count u such that k->u->j
  L2kt - for each k->i neq j: j->k, count u such that i->u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espITP_change(L, subroutine_path, subroutine_focus)             \
  int L2th;                                                             \
  if(spcache) L2th = GETDMUI(head,tail,spcache); else L2th=0;           \
  /* step through outedges of head (i.e., k: h->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      int L2hk;                                                         \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        if(spcache) L2hk = GETDMUI(k,head,spcache);                     \
	else{                                                           \
	  /*We have a h->k->t two-path, so add it to our count.*/       \
	  L2th++;                                                       \
	  L2hk=0;                                                       \
	  /*Now, count # u such that k->u->h (so that we know k's ESP value)*/ \
	  EXEC_THROUGH_OUTEDGES(k,f,u, {                                \
	      if(u!=head)                                               \
		L2hk+=IS_OUTEDGE(u,head); /*Increment if there is a cyclic edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2hk, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      int L2kt;                                                         \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        if(spcache) L2kt = GETDMUI(tail,k,spcache);                     \
	else{                                                           \
	  L2kt=0;                                                       \
	  /*Now, count # u such that t->u->k (so that we know k's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=tail)                                               \
		L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2kt, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(L2th, subroutine_focus);


/*
  Changescore for ESPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  OSP:
  L2th - count t->k, h->k
  L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
  L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espOSP_change(L, subroutine_path, subroutine_focus)             \
  int L2th;                                                             \
  if(spcache) L2th = GETDMUI(tail,head,spcache); else L2th=0;           \
  /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/        \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
        if(!spcache)                                                    \
	  /*Do we have a t->k,h->k SP?  If so, add it to our count.*/   \
	  L2th+=IS_OUTEDGE(head,k);                                     \
                                                                        \
	if(IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk;                                                     \
          if(spcache) L2tk = GETDMUI(tail,k,spcache);                   \
	  else{                                                         \
	    L2tk=0;                                                     \
	    /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/ \
	    EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
		if(u!=tail)                                             \
		  L2tk+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/ \
	      });                                                       \
	  }                                                             \
          call_subroutine_path(L2tk, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t, k->h, k!=h)*/         \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if((k!=head)&&(IS_OUTEDGE(k,head))){ /*Only stats that could change*/ \
        int L2kt;                                                       \
        if(spcache) L2kt = GETDMUI(k,tail,spcache);                     \
	else{                                                           \
	  L2kt=0;                                                       \
	  /*Now, count # u such that t->u,k->u (to get k->t's ESP value)*/ \
	  EXEC_THROUGH_OUTEDGES(k,f,u, {                                \
	      if(u!=tail)                                               \
		L2kt+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2kt, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(L2th, subroutine_focus);


/*
  Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  ISP:
  L2th - count k->t, k->h
  L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
  L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espISP_change(L, subroutine_path, subroutine_focus)             \
  int L2th;                                                             \
  if(spcache) L2th = GETDMUI(tail,head,spcache); else L2th=0;           \
  /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/         \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
        if(!spcache)                                                    \
	  /*Do we have a k->t,k->h SP?  If so, add it to our count.*/   \
	  L2th+=IS_OUTEDGE(k,tail);                                     \
                                                                        \
	if(IS_OUTEDGE(tail,k)){ /*Only consider stats that could change*/ \
          int L2kh;                                                     \
          if(spcache) L2kh = GETDMUI(k,head,spcache);                   \
	  else{                                                         \
	    L2kh=0;                                                     \
	    /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/ \
	    EXEC_THROUGH_INEDGES(k,f,u, {                               \
		if(u!=head)                                             \
		  L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
	      });                                                       \
	  }                                                             \
          call_subroutine_path(L2kh, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (i.e., k: h->k, t->k, k!=t)*/        \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if((k!=tail)&&(IS_OUTEDGE(tail,k))){ /*Only stats that could change*/ \
        int L2hk;                                                       \
        if(spcache) L2hk = GETDMUI(head,k,spcache);                     \
	else{                                                           \
	  L2hk=0;                                                       \
	  /*Now, count # u such that u->h,u->k (to get k->h's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=head)                                               \
		L2hk+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(L2hk, subroutine_path);                    \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(L2th, subroutine_focus);


/*
  Changescore for ESPs based on reciprocated two-paths, i.e. configurations for edge i->j such that i<->k and j<->k (with k!=j).

  RTP:
  L2th - count t<->k<->h
  L2kt - for each k->t neq h: h->t,k<->h, count u such that k<->u<->t
  L2tk - for each t->k neq h: h->t,k<->h, count u such that k<->u<->t
  L2kh - for each k->h neq t: h->t,k<->t, count u such that k<->u<->h
  L2hk - for each h->k neq t: h->t,k<->t, count u such that k<->u<->h

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espRTP_change(L, subroutine_path, subroutine_focus)             \
  int L2th; /*Two-path counts for various edges*/                       \
  if(spcache) L2th = GETDMUI(tail,head,spcache); else L2th=0;           \
  int htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/ \
  /* step through inedges of tail (k->t: k!=h,h->t,k<->h)*/             \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if(k!=head){                                                      \
        if(!spcache)                                                    \
          /*Do we have a t<->k<->h TP?  If so, add it to our count.*/   \
          L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)); \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2kt;                                                     \
          if(spcache) L2kt = GETDMUI(k,tail,spcache);                   \
          else{                                                         \
            L2kt=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=tail)&&(IS_OUTEDGE(u,k)))                        \
                  L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2kt, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of tail (t->k: k!=h,h->t,k<->h)*/            \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk;                                                     \
          if(spcache) L2tk = GETDMUI(tail,k,spcache);                   \
          else{                                                         \
            L2tk=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=tail)&&(IS_OUTEDGE(u,k)))                        \
                  L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2tk, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (k->h: k!=t,h->t,k<->t)*/             \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2kh;                                                     \
          if(spcache) L2kh = GETDMUI(k,head,spcache);                   \
          else{                                                         \
            L2kh=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get k->h's ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=head)&&(IS_OUTEDGE(u,k)))                        \
                  L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2kh, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (h->k: k!=t,h->t,k<->t)*/            \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2hk;                                                     \
          if(spcache) L2hk = GETDMUI(head,k,spcache);                   \
          else{                                                         \
            L2hk=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get h->k's ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=head)&&(IS_OUTEDGE(u,k)))                        \
                  L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(L2hk, subroutine_path);                  \
        }                                                               \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(L2th, subroutine_focus);

all_calcs(espUTP)
all_calcs(espOTP)
all_calcs(espITP)
all_calcs(espOSP)
all_calcs(espISP)
all_calcs(espRTP)


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
C_CHANGESTAT_FN(c_desp) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: espUTP_calc(sp_args); break;
  case ESPOTP: espOTP_calc(sp_args); break;
  case ESPITP: espITP_calc(sp_args); break;
  case ESPRTP: espRTP_calc(sp_args); break;
  case ESPOSP: espOSP_calc(sp_args); break;
  case ESPISP: espISP_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_despdist) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: espUTP_dist_calc(spd_args); break;
  case ESPOTP: espOTP_dist_calc(spd_args); break;
  case ESPITP: espITP_dist_calc(spd_args); break;
  case ESPRTP: espRTP_dist_calc(spd_args); break;
  case ESPOSP: espOSP_dist_calc(spd_args); break;
  case ESPISP: espISP_dist_calc(spd_args); break;
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
C_CHANGESTAT_FN(c_dgwesp) { 
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: cumchange = espUTP_gw_calc(gwsp_args); break;
  case ESPOTP: cumchange = espOTP_gw_calc(gwsp_args); break;
  case ESPITP: cumchange = espITP_gw_calc(gwsp_args); break;
  case ESPRTP: cumchange = espRTP_gw_calc(gwsp_args); break;
  case ESPOSP: cumchange = espOSP_gw_calc(gwsp_args); break;
  case ESPISP: cumchange = espISP_gw_calc(gwsp_args); break;
  }
  
  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
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
#define NEGATE_CHANGE_STATS for(unsigned int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] *= -1;
C_CHANGESTAT_FN(c_dnsp) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/
  
  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP: 
    espUTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspUTP_calc(sp_args);
    break;
  case ESPOTP: 
    espOTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOTP_calc(sp_args);
    break;
  case ESPITP: 
    espITP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspITP_calc(sp_args);
    break;
  case ESPRTP: 
    espRTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspRTP_calc(sp_args);
    break;
  case ESPOSP: 
    espOSP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOSP_calc(sp_args);
    break;
  case ESPISP: 
    espISP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspISP_calc(sp_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_dnspdist) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/

  /*Obtain the ESP changescores (by type)*/
  switch(type){
  case ESPUTP:
    espUTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspUTP_dist_calc(spd_args);
    break;
  case ESPOTP:
    espOTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOTP_dist_calc(spd_args);
    break;
  case ESPITP:
    espITP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspITP_dist_calc(spd_args);
    break;
  case ESPRTP:
    espRTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspRTP_dist_calc(spd_args);
    break;
  case ESPOSP:
    espOSP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOSP_dist_calc(spd_args);
    break;
  case ESPISP:
    espISP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspISP_dist_calc(spd_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
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
C_CHANGESTAT_FN(c_dgwnsp) {
  /*Set things up*/
  StoreDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double oneexpa = 1.0-exp(-alpha);    /*Precompute (1-exp(-alpha))*/
  int type = IINPUT_PARAM[0];     /*Get the ESP type code to be used*/
  double cumchange = 0;

  /*Obtain the changescores (by type)*/
  switch(type){
  case ESPUTP:
    cumchange = dspUTP_gw_calc(gwsp_args) - espUTP_gw_calc(gwsp_args);
    break;
  case ESPOTP:
    cumchange = dspOTP_gw_calc(gwsp_args) - espOTP_gw_calc(gwsp_args);
    break;
  case ESPITP:
    cumchange = dspITP_gw_calc(gwsp_args) - espITP_gw_calc(gwsp_args);
    break;
  case ESPRTP:
    cumchange = dspRTP_gw_calc(gwsp_args) - espRTP_gw_calc(gwsp_args);
    break;
  case ESPOSP:
    cumchange = dspOSP_gw_calc(gwsp_args) - espOSP_gw_calc(gwsp_args);
    break;
  case ESPISP:
    cumchange = dspISP_gw_calc(gwsp_args) - espISP_gw_calc(gwsp_args);
    break;
  }

  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
}


/*****************
 changestat: c_ddspbwrap
*****************/
C_CHANGESTAT_FN(c_ddspbwrap) {
  c_ddsp(tail, head, mtp, nwp, edgestate);

  // correct for double counting of directed vs. undirected dyads
  for(int ind = 0; ind < N_CHANGE_STATS; ind++) CHANGE_STAT[ind] /= 2.0;
}


C_CHANGESTAT_FN(c_ddspdistbwrap) {
  c_ddspdist(tail, head, mtp, nwp, edgestate);

  // correct for double counting of directed vs. undirected dyads
  for(int ind = 0; ind < N_CHANGE_STATS; ind++) CHANGE_STAT[ind] /= 2.0;
}


/*****************
 changestat: c_dgwdspbwrap
*****************/
C_CHANGESTAT_FN(c_dgwdspbwrap) {
  c_dgwdsp(tail, head, mtp, nwp, edgestate);
  
  // correct for double counting of directed vs. undirected dyads
  CHANGE_STAT[0] /= 2.0;
}

