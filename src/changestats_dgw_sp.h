/*  File src/changestats_dgw_sp.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2024 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_H_
#define _CHANGESTATS_DGW_SP_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_dyad_hashmap.h"

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5



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
        int L2tu = GETUDMUI(tail,u,spcache);                    \
        call_subroutine_path(tu, subroutine_path);              \
      }                                                         \
    });                                                         \
  EXEC_THROUGH_EDGES(tail,e,u, {                                \
      if (u!=head){                                             \
        int L2uh = GETUDMUI(u,head,spcache);                    \
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
        int L2tk = GETDDMUI(tail,k,spcache);                            \
        call_subroutine_path(tk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if (k!=head){ /*Only use contingent cases*/                       \
        int L2kh = GETDDMUI(k,head,spcache);                            \
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
        int L2kt = GETDDMUI(tail,k,spcache); /* spcache is an OTP cache. */ \
        call_subroutine_path(kt, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail, e, k, {                                    \
      if((k!=head)){ /*Only use contingent cases*/                      \
        int L2hk = GETDDMUI(k,head,spcache);                            \
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
        int L2tk = GETUDMUI(tail,k,spcache);                            \
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
      if(k!=head){                                                      \
        int L2kh = GETUDMUI(k,head,spcache);                            \
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
          int L2kh = GETUDMUI(k,head,spcache);                          \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      });                                                               \
    /* step through reciprocated outedges of tail (t->k: k!=h,k<-t)*/   \
    EXEC_THROUGH_OUTEDGES(head,e,k,{                                    \
        if(k!=tail&&IS_OUTEDGE(k,head)){                                \
          int L2kt = GETUDMUI(k,tail,spcache);                          \
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
  int L2th = GETDDMUI(tail,head,spcache);                       \
  /* step through outedges of head */                           \
  EXEC_THROUGH_EDGES(head,e,u, {                                \
      if (IS_UNDIRECTED_EDGE(u,tail) != 0){                     \
        int L2tu = GETUDMUI(tail,u,spcache);                    \
        int L2uh = GETUDMUI(u,head,spcache);                    \
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
  int L2th = GETDDMUI(tail,head,spcache);                               \
  /* step through outedges of tail (i.e., k: t->k)*/                    \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        int L2tk = GETDDMUI(tail,k,spcache);                            \
        call_subroutine_path(tk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (i.e., k: k->h)*/                     \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        int L2kh = GETDDMUI(k,head,spcache);                            \
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
  int L2th = GETDDMUI(head,tail,spcache);                               \
  /* step through outedges of head (i.e., k: h->k)*/                    \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if((k!=tail)&&(IS_OUTEDGE(k,tail))){ /*Only use contingent cases*/ \
        int L2hk = GETDDMUI(k,head,spcache);                            \
        call_subroutine_path(hk, subroutine_path);                      \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t)*/                     \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if((k!=head)&&(IS_OUTEDGE(head,k))){ /*Only use contingent cases*/ \
        int L2kt = GETDDMUI(tail,k,spcache);                            \
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
  int L2th = GETUDMUI(tail,head,spcache);                               \
  /* step through outedges of tail (i.e., k: t->k, k->h, k!=h)*/        \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
	if(IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk = GETUDMUI(tail,k,spcache);                  \
          call_subroutine_path(tk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of tail (i.e., k: k->t, k->h, k!=h)*/         \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if((k!=head)&&(IS_OUTEDGE(k,head))){ /*Only stats that could change*/ \
        int L2kt = GETUDMUI(k,tail,spcache);                            \
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
  int L2th = GETUDMUI(tail,head,spcache);                               \
  /* step through inedges of head (i.e., k: k->h, t->k, k!=t)*/         \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
	if(IS_OUTEDGE(tail,k)){ /*Only consider stats that could change*/ \
          int L2kh = GETUDMUI(k,head,spcache);                          \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (i.e., k: h->k, t->k, k!=t)*/        \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if((k!=tail)&&(IS_OUTEDGE(tail,k))){ /*Only stats that could change*/ \
        int L2hk = GETUDMUI(head,k,spcache);                            \
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
  int L2th = GETDDMUI(tail,head,spcache);                               \
  int htedge=IS_OUTEDGE(head,tail);  /*Is there an h->t (reciprocating) edge?*/ \
  /* step through inedges of tail (k->t: k!=h,h->t,k<->h)*/             \
  EXEC_THROUGH_INEDGES(tail,e,k, {                                      \
      if(k!=head){                                                      \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2kt = GETUDMUI(k,tail,spcache);                          \
          call_subroutine_path(kt, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of tail (t->k: k!=h,h->t,k<->h)*/            \
  EXEC_THROUGH_OUTEDGES(tail,e,k, {                                     \
      if(k!=head){                                                      \
        if(htedge&&IS_OUTEDGE(head,k)&&IS_OUTEDGE(k,head)){ /*Only consider stats that could change*/ \
          int L2tk = GETUDMUI(tail,k,spcache);                          \
          call_subroutine_path(tk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through inedges of head (k->h: k!=t,h->t,k<->t)*/             \
  EXEC_THROUGH_INEDGES(head,e,k, {                                      \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2kh = GETUDMUI(k,head,spcache);                          \
          call_subroutine_path(kh, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  /* step through outedges of head (h->k: k!=t,h->t,k<->t)*/            \
  EXEC_THROUGH_OUTEDGES(head,e,k, {                                     \
      if(k!=tail){                                                      \
        if(htedge&&IS_OUTEDGE(tail,k)&&IS_OUTEDGE(k,tail)){ /*Only consider stats that could change*/ \
          int L2hk = GETUDMUI(head,k,spcache);                          \
          call_subroutine_path(hk, subroutine_path);                    \
        }                                                               \
      }                                                                 \
    });                                                                 \
  call_subroutine_focus(th, subroutine_focus);

#endif // _CHANGESTATS_DGW_SP_H_
