/*  File src/changestats_dgw_sp.h in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_H_
#define _CHANGESTATS_DGW_SP_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

typedef enum {L2UTP, L2OTP, L2ITP, L2RTP, L2OSP, L2ISP} L2Type;

#define call_subroutine_path(count, subroutine_path)    \
  {int L2 = (L2 ## count);                              \
    {subroutine_path}}

#define call_subroutine_focus(count, subroutine_focus)  \
  {int L2 = (L2 ## count);                              \
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

#define dspUTP_change(subroutine_path, subroutine_focus)        \
  /* step through edges of head */                              \
  EXEC_THROUGH_EDGES(head,e,u, {                                \
      if (u!=tail){                                             \
        int L2tu;                                               \
        if(spcache) L2tu = GETUDMUI(tail,u,spcache);            \
        else{                                                   \
          L2tu=0;                                               \
          /* step through edges of u */                         \
          EXEC_THROUGH_EDGES(u,f,v, {                           \
              if(IS_UNDIRECTED_EDGE(v,tail)!= 0) L2tu++;        \
            });                                                 \
        }                                                       \
        call_subroutine_path(tu, subroutine_path);              \
      }                                                         \
    });                                                         \
  EXEC_THROUGH_EDGES(tail,e,u, {                                \
      if (u!=head){                                             \
        int L2uh;                                               \
        if(spcache) L2uh = GETUDMUI(u,head,spcache);            \
        else{                                                   \
          L2uh=0;                                               \
          /* step through edges of u */                         \
          EXEC_THROUGH_EDGES(u,f,v, {                           \
              if(IS_UNDIRECTED_EDGE(v,head)!= 0) L2uh++;        \
            });                                                 \
        }                                                       \
        call_subroutine_path(uh, subroutine_path);              \
      }                                                         \
    });


/*
  Changescore for dsps based on outgoing two-paths, i.e. configurations for non-edge i->j such that i->k->j.

  This function should only be used in the directed case
*/

#define dspOTP_change(subroutine_path, subroutine_focus)                \
  /* step through outedges of head (i.e., k: t->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head, e, k, {                                   \
      if(k!=tail){ /*Only use contingent cases*/                        \
        int L2tk;                                                       \
        if(spcache) L2tk = GETDDMUI(tail,k,spcache);                    \
	else{                                                           \
	  L2tk=0;                                                       \
	  /* step through inedges of k, incl. (head,k) itself */        \
	  EXEC_THROUGH_INEDGES(k, f, u, {                               \
	      L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(tk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if (k!=head){ /*Only use contingent cases*/                       \
        int L2kh;                                                       \
	if(spcache) L2kh = GETDDMUI(k,head,spcache);                    \
	else{                                                           \
	  L2kh=0;                                                       \
	  /* step through outedges of k , incl. (k,tail) itself */      \
	  EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
	      L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(kh, subroutine_path);                      \
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
#define dspITP_change(subroutine_path, subroutine_focus)                \
  /* step through outedges of head (i.e., k: h->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head, e, k, {                                   \
      if((k!=tail)){ /*Only use contingent cases*/                      \
        int L2kt;                                                       \
        /*We have a h->k->t two-path, so add it to our count.*/         \
        if(spcache) L2kt = GETDDMUI(tail,k,spcache); /* spcache is an OTP cache. */ \
        else{                                                           \
          L2kt=0;                                                       \
          /*Now, count # u such that k->u->h (so that we know k's ESP value)*/ \
          EXEC_THROUGH_INEDGES(k, f, u, {                               \
              L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(kt, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if((k!=head)){ /*Only use contingent cases*/                      \
        int L2hk;                                                       \
        if(spcache) L2hk = GETDDMUI(k,head,spcache);                    \
        else{                                                           \
          L2hk=0;                                                       \
          /*Now, count # u such that t->u->k (so that we know k's ESP value)*/ \
          EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
              L2hk+=IS_OUTEDGE(u,head); /*Increment if there is a cyclic edge*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(hk, subroutine_path);                      \
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
#define dspOSP_change(subroutine_path, subroutine_focus)                \
  /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/        \
  EXEC_THROUGH_INEDGES(head, e, k, {                                    \
      if(k!=tail){                                                      \
        int L2tk;                                                       \
        /*Do we have a t->k,h->k SP?  If so, add it to our count.*/     \
        if(spcache) L2tk = GETUDMUI(tail,k,spcache);                    \
        else{                                                           \
          L2tk=0;                                                       \
          /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/ \
          EXEC_THROUGH_OUTEDGES(k, f, u, {                              \
              if (u != tail)                                            \
                /*Increment if there is an OSP  */                      \
                L2tk+=(IS_OUTEDGE(tail,u));                             \
            });                                                         \
        }                                                               \
        call_subroutine_path(tk, subroutine_path);                      \
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
#define dspISP_change(subroutine_path, subroutine_focus)                \
  /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/         \
  EXEC_THROUGH_OUTEDGES(tail, e, k, {                                   \
      int L2kh;                                                         \
      if(k!=head){                                                      \
        if(spcache) L2kh = GETUDMUI(k,head,spcache);                    \
        else{                                                           \
          L2kh=0;                                                       \
          /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/ \
          EXEC_THROUGH_INEDGES(k, f, u, {                               \
              if(u!=head)                                               \
                L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
            });                                                         \
        }                                                               \
        call_subroutine_path(kh, subroutine_path);                      \
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
#define dspRTP_change(subroutine_path, subroutine_focus)                \
  int htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/ \
  if(htedge){ /* Otherwise, t->h doesn't make a difference. */          \
    /* step through reciprocated outedges of tail (t->k: k!=h,k<-t)*/   \
    EXEC_THROUGH_OUTEDGES(tail,e,k,{                                    \
        if(k!=head&&IS_OUTEDGE(k,tail)){                                \
          int L2kh;                                                     \
          if(spcache) L2kh = GETUDMUI(k,head,spcache);                  \
          else{                                                         \
            L2kh=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get k->h's SP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u,{                               \
                if(u!=tail&&u!=head&&(IS_OUTEDGE(u,k)))                 \
                  L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      });                                                               \
    /* step through reciprocated outedges of tail (t->k: k!=h,k<-t)*/   \
    EXEC_THROUGH_OUTEDGES(head,e,k,{                                    \
        if(k!=tail&&IS_OUTEDGE(k,head)){                                \
          int L2kt;                                                     \
          if(spcache) L2kt = GETUDMUI(k,tail,spcache);                  \
          else{                                                         \
            L2kt=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get k->t's SP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u,{                               \
                if(u!=head&&u!=tail&&(IS_OUTEDGE(u,k)))                 \
                  L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(kt, subroutine_path);                    \
        }                                                               \
      });                                                               \
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
#define espUTP_change(subroutine_path, subroutine_focus)        \
  int L2th;                                                     \
  if(spcache) L2th = GETDDMUI(tail,head,spcache); else L2th=0;  \
  /* step through outedges of head */                           \
  EXEC_THROUGH_EDGES(head,e,u, {                                \
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){                     \
        int L2tu;                                               \
        int L2uh;                                               \
	if(spcache){                                            \
	  L2tu = GETUDMUI(tail,u,spcache);                      \
	  L2uh = GETUDMUI(u,head,spcache);                      \
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
        call_subroutine_path(tu, subroutine_path);              \
        call_subroutine_path(uh, subroutine_path);              \
      }                                                         \
    });                                                         \
  call_subroutine_focus(th, subroutine_focus);


/*
  Changescore for ESPs based on outgoing two-paths, i.e. configurations for edge i->j such that i->k->j.

  OTP:
  L2th - count i->k->j
  L2tk - for each i->k neq j: j->k, count u such that i->u->k
  L2kh - for each k->j neq i: k->i, count u such that k->u->j

  This function should only be used in the directed case, with espUTP being used in the undirected case.
*/
#define espOTP_change(subroutine_path, subroutine_focus)                \
  int L2th;                                                             \
  if(spcache) L2th = GETDDMUI(tail,head,spcache); else L2th=0;          \
  /* step through outedges of tail (i.e., k: t->k)*/                    \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(!spcache&&(k!=head)&&(IS_OUTEDGE(k,head))){                    \
        /*We have a t->k->h two-path, so add it to our count.*/         \
        L2th++;                                                         \
      }                                                                 \
      int L2tk;                                                         \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        if(spcache) L2tk = GETDDMUI(tail,k,spcache);                    \
	else{                                                           \
	  L2tk=0;                                                       \
	  /*Now, count # u such that t->u->k (to find t->k's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=tail)                                               \
		L2tk+=IS_OUTEDGE(tail,u); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(tk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (i.e., k: k->h)*/                     \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      int L2kh;                                                         \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        if(spcache) L2kh = GETDDMUI(k,head,spcache);                    \
	else{                                                           \
	  L2kh=0;                                                       \
	  /*Now, count # u such that k->u->j (to find k->h's ESP value)*/ \
	  EXEC_THROUGH_OUTEDGES(k,f,u, {                                \
	      if(u!=head)                                               \
		L2kh+=IS_OUTEDGE(u,head); /*Increment if there is a trans edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(kh, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);


/*
  Changescore for ESPs based on incoming two-paths, i.e. configurations for edge i->j such that j->k->i.

  ITP:
  L2th - count j->k->i
  L2hk - for each j->k neq i: k->i, count u such that k->u->j
  L2kt - for each k->i neq j: j->k, count u such that i->u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espITP_change(subroutine_path, subroutine_focus)                \
  int L2th;                                                             \
  if(spcache) L2th = GETDDMUI(head,tail,spcache); else L2th=0;          \
  /* step through outedges of head (i.e., k: h->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      int L2hk;                                                         \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        if(spcache) L2hk = GETDDMUI(k,head,spcache);                    \
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
        call_subroutine_path(hk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      int L2kt;                                                         \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        if(spcache) L2kt = GETDDMUI(tail,k,spcache);                    \
	else{                                                           \
	  L2kt=0;                                                       \
	  /*Now, count # u such that t->u->k (so that we know k's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=tail)                                               \
		L2kt+=IS_OUTEDGE(tail,u); /*Increment if there is a cyclic edge*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(kt, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);


/*
  Changescore for ESPs based on outgoing shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  OSP:
  L2th - count t->k, h->k
  L2tk - for each t->k neq h: k->h, count u such that t->u, k->u
  L2kt - for each k->t neq h: k->h, count u such that t->u, k->u

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espOSP_change(subroutine_path, subroutine_focus)                \
  int L2th;                                                             \
  if(spcache) L2th = GETUDMUI(tail,head,spcache); else L2th=0;          \
  /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/        \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
        if(!spcache)                                                    \
	  /*Do we have a t->k,h->k SP?  If so, add it to our count.*/   \
	  L2th+=IS_OUTEDGE(head,k);                                     \
                                                                        \
	if(IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk;                                                     \
          if(spcache) L2tk = GETUDMUI(tail,k,spcache);                  \
	  else{                                                         \
	    L2tk=0;                                                     \
	    /*Now, count # u such that t->u,k->u (to get t->k's ESP value)*/ \
	    EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
		if(u!=tail)                                             \
		  L2tk+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/ \
	      });                                                       \
	  }                                                             \
          call_subroutine_path(tk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t, k->h, k!=h)*/         \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if((k!=head)&&(IS_OUTEDGE(k,head))){ /*Only stats that could change*/ \
        int L2kt;                                                       \
        if(spcache) L2kt = GETUDMUI(k,tail,spcache);                    \
	else{                                                           \
	  L2kt=0;                                                       \
	  /*Now, count # u such that t->u,k->u (to get k->t's ESP value)*/ \
	  EXEC_THROUGH_OUTEDGES(k,f,u, {                                \
	      if(u!=tail)                                               \
		L2kt+=IS_OUTEDGE(tail,u);  /*Increment if there is an OSP*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(kt, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);


/*
  Changescore for ESPs based on incoming shared partners, i.e. configurations for edge i->j such that i->k and j->k (with k!=j).

  ISP:
  L2th - count k->t, k->h
  L2hk - for each h->k neq t: t->k, count u such that u->h, u->k
  L2kh - for each k->h neq t: t->k, count u such that u->h, u->k

  We assume that this is only called for directed graphs - otherwise, use the baseline espUTP function.
*/
#define espISP_change(subroutine_path, subroutine_focus)                \
  int L2th;                                                             \
  if(spcache) L2th = GETUDMUI(tail,head,spcache); else L2th=0;          \
  /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/         \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
        if(!spcache)                                                    \
	  /*Do we have a k->t,k->h SP?  If so, add it to our count.*/   \
	  L2th+=IS_OUTEDGE(k,tail);                                     \
                                                                        \
	if(IS_OUTEDGE(tail,k)){ /*Only consider stats that could change*/ \
          int L2kh;                                                     \
          if(spcache) L2kh = GETUDMUI(k,head,spcache);                  \
	  else{                                                         \
	    L2kh=0;                                                     \
	    /*Now, count # u such that u->h,u->k (to get h>k's ESP value)*/ \
	    EXEC_THROUGH_INEDGES(k,f,u, {                               \
		if(u!=head)                                             \
		  L2kh+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
	      });                                                       \
	  }                                                             \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (i.e., k: h->k, t->k, k!=t)*/        \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if((k!=tail)&&(IS_OUTEDGE(tail,k))){ /*Only stats that could change*/ \
        int L2hk;                                                       \
        if(spcache) L2hk = GETUDMUI(head,k,spcache);                    \
	else{                                                           \
	  L2hk=0;                                                       \
	  /*Now, count # u such that u->h,u->k (to get k->h's ESP value)*/ \
	  EXEC_THROUGH_INEDGES(k,f,u, {                                 \
	      if(u!=head)                                               \
		L2hk+=IS_OUTEDGE(u,head);  /*Increment if there is an ISP*/ \
	    });                                                         \
	}                                                               \
        call_subroutine_path(hk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);


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
#define espRTP_change(subroutine_path, subroutine_focus)                \
  int L2th; /*Two-path counts for various edges*/                       \
  if(spcache) L2th = GETDDMUI(tail,head,spcache); else L2th=0;          \
  int htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/ \
  /* step through inedges of tail (k->t: k!=h,h->t,k<->h)*/             \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if(k!=head){                                                      \
        if(!spcache)                                                    \
          /*Do we have a t<->k<->h TP?  If so, add it to our count.*/   \
          L2th+=(IS_OUTEDGE(tail,k)&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)); \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2kt;                                                     \
          if(spcache) L2kt = GETUDMUI(k,tail,spcache);                  \
          else{                                                         \
            L2kt=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get (k,t)'s ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=tail)&&(IS_OUTEDGE(u,k)))                        \
                  L2kt+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(kt, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of tail (t->k: k!=h,h->t,k<->h)*/            \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk;                                                     \
          if(spcache) L2tk = GETUDMUI(tail,k,spcache);                  \
          else{                                                         \
            L2tk=0;                                                     \
            /*Now, count # u such that k<->u<->t (to get (tk)'s ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=tail)&&(IS_OUTEDGE(u,k)))                        \
                  L2tk+=(IS_OUTEDGE(u,tail)&&IS_OUTEDGE(tail,u));  /*k<->u<->t?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(tk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (k->h: k!=t,h->t,k<->t)*/             \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2kh;                                                     \
          if(spcache) L2kh = GETUDMUI(k,head,spcache);                  \
          else{                                                         \
            L2kh=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get k->h's ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=head)&&(IS_OUTEDGE(u,k)))                        \
                  L2kh+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (h->k: k!=t,h->t,k<->t)*/            \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2hk;                                                     \
          if(spcache) L2hk = GETUDMUI(head,k,spcache);                  \
          else{                                                         \
            L2hk=0;                                                     \
            /*Now, count # u such that k<->u<->h (to get h->k's ESP value)*/ \
            EXEC_THROUGH_OUTEDGES(k,f,u, {                              \
                if((u!=head)&&(IS_OUTEDGE(u,k)))                        \
                  L2hk+=(IS_OUTEDGE(u,head)&&IS_OUTEDGE(head,u));  /*k<->u<->h?*/ \
              });                                                       \
          }                                                             \
          call_subroutine_path(hk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);

#ifdef __cplusplus
#include "cpp/ergm_network.h"

namespace ergm {
inline namespace v1 {
namespace sp {

template<typename NetworkView>
inline int count_otp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.in_neighbors(head)) count += nw(tail, u);
  return count;
}

template<typename NetworkView>
inline int count_utp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.neighbors(head)) count += nw(u, tail);
  return count;
}

template<typename NetworkView>
inline int count_osp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.out_neighbors(head))
    if(u != tail) count += nw(tail, u);
  return count;
}

template<typename NetworkView>
inline int count_isp(NetworkView& nw, Vertex tail, Vertex head){
  int count = 0;
  for(auto u: nw.in_neighbors(head))
    if(u != tail) count += nw(u, tail);
  return count;
}

template<typename NetworkView>
inline int count_rtp(NetworkView& nw, Vertex tail, Vertex head, Vertex exclude1, Vertex exclude2){
  int count = 0;
  for(auto u: nw.out_neighbors(tail))
    if(u != exclude1 && u != exclude2 && nw(u, tail))
      count += nw(u, head) && nw(head, u);
  return count;
}

template<L2Type type, typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void dsp_change(Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus){
  if constexpr(type == L2UTP){
    for(auto u: nw.neighbors(head))
      if(u != tail)
        update_path(spcache ? GETUDMUI(tail, u, spcache) : count_utp(nw, tail, u));

    for(auto u: nw.neighbors(tail))
      if(u != head)
        update_path(spcache ? GETUDMUI(u, head, spcache) : count_utp(nw, u, head));
  }else if constexpr(type == L2OTP){
    for(auto k: nw.out_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    for(auto k: nw.in_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
  }else if constexpr(type == L2ITP){
    for(auto k: nw.out_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    for(auto k: nw.in_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
  }else if constexpr(type == L2RTP){
    if(nw(head, tail)){
      for(auto k: nw.out_neighbors(tail))
        if(k != head && nw(k, tail))
          update_path(spcache ? GETUDMUI(k, head, spcache) : count_rtp(nw, k, head, tail, head));

      for(auto k: nw.out_neighbors(head))
        if(k != tail && nw(k, head))
          update_path(spcache ? GETUDMUI(k, tail, spcache) : count_rtp(nw, k, tail, head, tail));
    }
  }else if constexpr(type == L2OSP){
    for(auto k: nw.in_neighbors(head))
      if(k != tail)
        update_path(spcache ? GETUDMUI(tail, k, spcache) : count_osp(nw, tail, k));
  }else if constexpr(type == L2ISP){
    for(auto k: nw.out_neighbors(tail))
      if(k != head)
        update_path(spcache ? GETUDMUI(k, head, spcache) : count_isp(nw, head, k));
  }
}

template<typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void dsp_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  switch(type){
  case L2UTP: dsp_change<L2UTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OTP: dsp_change<L2OTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ITP: dsp_change<L2ITP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2RTP: dsp_change<L2RTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OSP: dsp_change<L2OSP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ISP: dsp_change<L2ISP>(tail, head, nw, spcache, update_path, update_focus); break;
  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }
}

template<L2Type type, typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void esp_change(Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  if constexpr(type == L2UTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    for(auto u: nw.neighbors(head))
      if(nw(u, tail)){
        if(!spcache) L2th++;
        update_path(spcache ? GETUDMUI(tail, u, spcache) : count_utp(nw, tail, u));
        update_path(spcache ? GETUDMUI(u, head, spcache) : count_utp(nw, u, head));
      }

    update_focus(L2th);
  }else if constexpr(type == L2OTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;

    for(auto k: nw.out_neighbors(tail)){
      if(!spcache && k != head && nw(k, head)) L2th++;
      if(k != head && nw(head, k))
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));
    }

    for(auto k: nw.in_neighbors(head))
      if(k != tail && nw(k, tail))
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));

    update_focus(L2th);
  }else if constexpr(type == L2ITP){
    int L2th = spcache ? GETDDMUI(head, tail, spcache) : 0;

    for(auto k: nw.out_neighbors(head))
      if(k != tail && nw(k, tail)){
        if(!spcache) L2th++;
        update_path(spcache ? GETDDMUI(k, head, spcache) : count_otp(nw, k, head));
      }

    for(auto k: nw.in_neighbors(tail))
      if(k != head && nw(head, k))
        update_path(spcache ? GETDDMUI(tail, k, spcache) : count_otp(nw, tail, k));

    update_focus(L2th);
  }else if constexpr(type == L2RTP){
    int L2th = spcache ? GETDDMUI(tail, head, spcache) : 0;
    bool htedge = nw(head, tail);

    for(auto k: nw.in_neighbors(tail)){
      if(k != head){
        if(!spcache) L2th += nw(tail, k) && nw(head, k) && nw(k, head);
        if(htedge && nw(head, k) && nw(k, head))
          update_path(spcache ? GETUDMUI(k, tail, spcache) : count_rtp(nw, k, tail, tail, 0));
      }
    }

    for(auto k: nw.out_neighbors(tail))
      if(k != head && htedge && nw(head, k) && nw(k, head))
        update_path(spcache ? GETUDMUI(tail, k, spcache) : count_rtp(nw, k, tail, tail, 0));

    for(auto k: nw.in_neighbors(head))
      if(k != tail && htedge && nw(tail, k) && nw(k, tail))
        update_path(spcache ? GETUDMUI(k, head, spcache) : count_rtp(nw, k, head, head, 0));

    for(auto k: nw.out_neighbors(head))
      if(k != tail && htedge && nw(tail, k) && nw(k, tail))
        update_path(spcache ? GETUDMUI(head, k, spcache) : count_rtp(nw, k, head, head, 0));

    update_focus(L2th);
  }else if constexpr(type == L2OSP){
    int L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    for(auto k: nw.out_neighbors(tail))
      if(k != head){
        if(!spcache) L2th += nw(head, k);
        if(nw(k, head))
          update_path(spcache ? GETUDMUI(tail, k, spcache) : count_osp(nw, tail, k));
      }

    for(auto k: nw.in_neighbors(tail))
      if(k != head && nw(k, head))
        update_path(spcache ? GETUDMUI(k, tail, spcache) : count_osp(nw, tail, k));

    update_focus(L2th);
  }else if constexpr(type == L2ISP){
    int L2th = spcache ? GETUDMUI(tail, head, spcache) : 0;

    for(auto k: nw.in_neighbors(head))
      if(k != tail){
        if(!spcache) L2th += nw(k, tail);
        if(nw(tail, k))
          update_path(spcache ? GETUDMUI(k, head, spcache) : count_isp(nw, head, k));
      }

    for(auto k: nw.out_neighbors(head))
      if(k != tail && nw(tail, k))
        update_path(spcache ? GETUDMUI(head, k, spcache) : count_isp(nw, head, k));

    update_focus(L2th);
  }
}

template<typename NetworkView, typename UpdatePath, typename UpdateFocus>
inline void esp_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, StoreStrictDyadMapUInt *spcache, UpdatePath update_path, UpdateFocus update_focus){
  switch(type){
  case L2UTP: esp_change<L2UTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OTP: esp_change<L2OTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ITP: esp_change<L2ITP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2RTP: esp_change<L2RTP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2OSP: esp_change<L2OSP>(tail, head, nw, spcache, update_path, update_focus); break;
  case L2ISP: esp_change<L2ISP>(tail, head, nw, spcache, update_path, update_focus); break;
  default: error("In ergm shared partner helper, an unsupported type of triad: %d.", type);
  }
}

template<typename NetworkView>
inline int dsp_nonzero_change(L2Type type, Vertex tail, Vertex head, NetworkView& nw, Rboolean edgestate, StoreStrictDyadMapUInt *spcache){
  int echange = edgestate ? -1 : 1;
  int delta = 0;
  dsp_change(type, tail, head, nw, spcache,
             [&](int L2){ delta += (L2 + echange != 0) - (L2 != 0); },
             [&](int){});
  return delta;
}

} // namespace sp
} // namespace v1
} // namespace ergm
#endif

#endif // _CHANGESTATS_DGW_SP_H_
