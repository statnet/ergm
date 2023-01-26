/*  File inst/include/ergm_wtchangestat.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2023 Statnet Commons
 */
#ifndef _ERGM_WTCHANGESTAT_H_
#define _ERGM_WTCHANGESTAT_H_

#include "ergm_wtedgetree.h"

typedef struct WtModelTermstruct {
  void (*c_func)(Vertex, Vertex, double, struct WtModelTermstruct*, WtNetwork*, double);
  void (*d_func)(Edge, Vertex*, Vertex*, double*, struct WtModelTermstruct*, WtNetwork*);
  void (*i_func)(struct WtModelTermstruct*, WtNetwork*);
  void (*u_func)(Vertex, Vertex, double, struct WtModelTermstruct*, WtNetwork*, double);
  void (*f_func)(struct WtModelTermstruct*, WtNetwork*);
  void (*s_func)(struct WtModelTermstruct*, WtNetwork*);
  SEXP (*w_func)(struct WtModelTermstruct*, WtNetwork*);  
  void (*x_func)(unsigned int type, void *data, struct WtModelTermstruct*, WtNetwork*);
  void (*z_func)(struct WtModelTermstruct*, WtNetwork*, Rboolean);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int *iattrib; /* Ptr to vector of integer covariates (if necessary; generally unused) */
  int nstats;   /* Number of change statistics to be returned */
  unsigned int statspos; /* Position of this term's stats in the workspace vector. */ 
  double *dstats; /* ptr to change statistics returned */
  int ninputparams; /* Number of double input parameters passed to function */
  double *inputparams; /* ptr to double input parameters passed */
  int niinputparams; /* Number of integer input parameters passed to function */
  int *iinputparams; /* ptr to integer input parameters passed */
  double *statcache; /* vector of the same length as dstats */
  double *emptynwstats; /* vector of the same length as dstats or NULL*/
  void *storage; /* optional space for persistent storage */
  void **aux_storage; /* optional space for persistent public (auxiliary) storage */
  unsigned int n_aux;
  unsigned int *aux_slots;
  SEXP R; /* R term object. */
  SEXP ext_state; /* A place from which to read extended state. */
} WtModelTerm;

/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

#include "ergm_changestat_common.do_not_include_directly.h"

/* tell whether a particular edge exists */
#define _WtIS_OUTEDGE2(a,b) _WtIS_OUTEDGE3(a,b,nwp)
#define _WtIS_OUTEDGE3(a,b,nwp) (EdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define WtIS_OUTEDGE(...) _GET_OVERRIDE3(__VA_ARGS__, _WtIS_OUTEDGE3, _WtIS_OUTEDGE2,)(__VA_ARGS__)

#define _WtIS_INEDGE2(a,b) _WtIS_INEDGE3(a,b,nwp)
#define _WtIS_INEDGE3(a,b,nwp) (EdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define WTIS_INEDGE(...) _GET_OVERRIDE3(__VA_ARGS__, _WtIS_INEDGE3, _WtIS_INEDGE2,)(__VA_ARGS__)

#define _WtIS_UNDIRECTED_EDGE2(a,b) _WtIS_UNDIRECTED_EDGE3(a,b,nwp)
#define _WtIS_UNDIRECTED_EDGE3(a,b,nwp) (EdgetreeSearch(MIN((a),(b)),MAX((a),(b)),nwp->outedges)!=0?1:0)
#define WtIS_UNDIRECTED_EDGE(...) _GET_OVERRIDE3(__VA_ARGS__, _WtIS_UNDIRECTED_EDGE3, _WtIS_UNDIRECTED_EDGE2,)(__VA_ARGS__)

/* Return the Edge number of the smallest-labelled neighbor of the node 
   labelled "a".  Or, return the Edge number of the next-largest neighbor 
   starting from the pointer "e", which points to a node in an edgetree. 
   Mostly, these are utility macros used by the STEP_THROUGH_OUTEDGES 
   and STEP_THROUGH_INEDGES macros. */
#define WtMIN_OUTEDGE(a) (WtEdgetreeMinimum(nwp->outedges, (a)))
#define WtMIN_INEDGE(a) (WtEdgetreeMinimum(nwp->inedges, (a)))
#define WtNEXT_OUTEDGE(e) (WtEdgetreeSuccessor(nwp->outedges,(e)))
#define WtNEXT_INEDGE(e) (WtEdgetreeSuccessor(nwp->inedges,(e)))
/* As WtNEXT_*EDGE, but visits the parent nodes before the child
   nodes. */
#define WtNEXT_OUTEDGE_PRE(e) (WtEdgetreePreSuccessor(nwp->outedges,(e)))
#define WtNEXT_INEDGE_PRE(e) (WtEdgetreePreSuccessor(nwp->inedges,(e)))

/* The OUTWT and INWT macros give the weight of edge e, depending
   on whether it is an in-edge or an out-edge.  Presumably the first endnode
   of the edge is already known in this context. */
#define OUTWT(e) (nwp->outedges[(e)].weight)
#define INWT(e) (nwp->inedges[(e)].weight)

/* Return each of the out-neighbors or in-neighbors, one at a time,
   of node a.  At each iteration of the loop, the variable v gives the node 
   number of the corresponding neighbor.  The e variable, which should be
   initialized as type Edge, is merely the looping variable. */
#define WtSTEP_THROUGH_OUTEDGES(a,e,v) for((e)=WtMIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=WtNEXT_OUTEDGE(e))
#define WtSTEP_THROUGH_INEDGES(a,e,v) for((e)=WtMIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=WtNEXT_INEDGE(e))

/* As WtSTEP_THROUGH_*EDGES, but visit the parent nodes before the child
   nodes. This is useful for "copying" an edgetree. */
#define WtSTEP_THROUGH_OUTEDGES_PRE(a,e,v) for((e)=(a); ((v)=OUTVAL(e))!=0;WtNEXT_OUTEDGE_PRE(e))
#define WtSTEP_THROUGH_INEDGES_PRE(a,e,v) for((e)=(a);((v)=INVAL(e))!=0;(e)=WtNEXT_INEDGE_PRE(e))

// These are "declaring" versions of the above, optimized for use in EXEC_TROUGH_*EDGES macros.
#define WtSTEP_THROUGH_OUTEDGES_DECL(a,e,v) Vertex v; for(Edge e=WtMIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;e=WtNEXT_OUTEDGE(e))
#define WtSTEP_THROUGH_INEDGES_DECL(a,e,v) Vertex v; for(Edge e=WtMIN_INEDGE(a);((v)=INVAL(e))!=0;e=WtNEXT_INEDGE(e))
#define WtSTEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) Vertex v; for(Edge e=(a);((v)=OUTVAL(e))!=0;e=WtNEXT_OUTEDGE_PRE(e))
#define WtSTEP_THROUGH_INEDGES_PRE_DECL(a,e,v) Vertex v; for(Edge e=(a);((v)=INVAL(e))!=0;e=WtNEXT_INEDGE_PRE(e))

/* Instead of stepping through execute "subroutine" for each neighbor
   automatically adapting to undirected networks. w gets the weight of
   the edge in question. */
/* NOTE: For some reason, GCC complains when it encounters multiple
   declaration in a single statement inside the subroutine, so
   double v1, v2;
   might fail, while
   double v1;
   double v2;
   works.*/
#define WtEXEC_THROUGH_OUTEDGES(a,e,v,w,subroutine) {if(DIRECTED){ WtSTEP_THROUGH_OUTEDGES_DECL(a,e,v) {double w=OUTWT(e); subroutine} } else { WtEXEC_THROUGH_EDGES(a,e,v,w,subroutine) }}
#define WtEXEC_THROUGH_INEDGES(a,e,v,w,subroutine) {if(DIRECTED){ WtSTEP_THROUGH_INEDGES_DECL(a,e,v) {double w=INWT(e); subroutine} } else { WtEXEC_THROUGH_EDGES(a,e,v,w,subroutine) }}
#define WtEXEC_THROUGH_EDGES(a,e,v,w,subroutine) { {WtSTEP_THROUGH_OUTEDGES_DECL(a,e,v) {double w=OUTWT(e); subroutine}};  {WtSTEP_THROUGH_INEDGES_DECL(a,e,v) {Vertex v=INVAL(e); double w=INWT(e); subroutine}}; }
#define WtEXEC_THROUGH_OUTEDGES_PRE(a,e,v,w,subroutine) {if(DIRECTED){ WtSTEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {double w=OUTWT(e); subroutine} } else { WtEXEC_THROUGH_EDGES_PRE(a,e,v,w,subroutine) }}
#define WtEXEC_THROUGH_INEDGES_PRE(a,e,v,w,subroutine) {if(DIRECTED){ WtSTEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {double w=INWT(e); subroutine} } else { WtEXEC_THROUGH_EDGES_PRE(a,e,v,w,subroutine) }}
#define WtEXEC_THROUGH_EDGES_PRE(a,e,v,w,subroutine) { {WtSTEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {double w=OUTWT(e); subroutine}};  {WtSTEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {double w=INWT(e); subroutine}}; }

/* Non-adaptive versions of the above. (I.e. ForceOUT/INEDGES.) */
#define WtEXEC_THROUGH_FOUTEDGES(a,e,v,w,subroutine) WtSTEP_THROUGH_OUTEDGES_DECL(a,e,v) {double w=OUTWT(e); subroutine}
#define WtEXEC_THROUGH_FINEDGES(a,e,v,w,subroutine) WtSTEP_THROUGH_INEDGES_DECL(a,e,v) {double w=INWT(e); subroutine}
#define WtEXEC_THROUGH_FOUTEDGES_PRE(a,e,v,w,subroutine) WtSTEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {double w=OUTWT(e); subroutine}
#define WtEXEC_THROUGH_FINEDGES_PRE(a,e,v,w,subroutine) WtSTEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {double w=INWT(e); subroutine}

/* Exectute through all edges (nonzero values) in the network. */
#define WtEXEC_THROUGH_NET_EDGES(a,b,e,w,subroutine) {for(Vertex a=1; a <= N_NODES; a++){WtEXEC_THROUGH_FOUTEDGES(a, e, b, w, {subroutine});}}
#define WtEXEC_THROUGH_NET_EDGES_PRE(a,b,e,w,subroutine) {for(Vertex a=1; a <= N_NODES; a++){WtEXEC_THROUGH_FOUTEDGES_PRE(a, e, b, w, {subroutine});}}

/* Get and set the weight of the (a,b) edge. */
#define _WtGETWT2(a,b) _WtGETWT3(a,b,nwp)
#define _WtGETWT3(a,b,nwp) (WtGetEdge(a,b,nwp))
#define WtGETWT(...) _GET_OVERRIDE3(__VA_ARGS__, _WtGETWT3, _WtGETWT2,)(__VA_ARGS__)
#define _WtSETWT3(a,b,w) _WtSETWT4(a,b,w,nwp)
#define _WtSETWT4(a,b,w,nwp) (WtSetEdge(a,b,w,nwp))
#define WtSETWT(...) _GET_OVERRIDE4(__VA_ARGS__, _WtSETWT4, _WtSETWT3,)(__VA_ARGS__)

/* Cycle through all toggles proposed for the current step, then
   make the current toggle in case of more than one proposed toggle, then
   undo all of the toggles to reset the original network state.  */
#define WtFOR_EACH_TOGGLE for(unsigned int WtTOGGLEIND=0; WtTOGGLEIND<ntoggles; WtTOGGLEIND++)
/* The idea here is to essentially swap the contents of the proposed
   weights with the current weights, and then swap them back when
   done. */
#define WtTAIL (tail_var)
#define WtHEAD (head_var)
#define WtNEWWT (newwt_var)
#define WtOLDWT (oldwt_var)

#define WtGETOLDTOGGLEINFO() Vertex WtTAIL=tails[WtTOGGLEIND], WtHEAD=heads[WtTOGGLEIND]; double WtOLDWT=WtGETWT(WtTAIL,WtHEAD);
#define WtGETTOGGLEINFO() WtGETOLDTOGGLEINFO(); double WtNEWWT=weights[WtTOGGLEIND];
#define WtGETNEWTOGGLEINFO() Vertex WtTAIL=tails[WtTOGGLEIND], WtHEAD=heads[WtTOGGLEIND]; double WtNEWWT=weights[WtTOGGLEIND];

/* SETWT_WITH_BACKUP(a) must be called _after_ GETTOGGLEINFO! */
#define WtSETWT_WITH_BACKUP() {WtSETWT(WtTAIL,WtHEAD,WtNEWWT); weights[WtTOGGLEIND]=WtOLDWT;}
#define WtUNDO_SETWT() {WtGETOLDTOGGLEINFO(); WtSETWT(WtTAIL,WtHEAD,weights[WtTOGGLEIND]); weights[WtTOGGLEIND]=WtOLDWT;}
#define WtIF_MORE_TO_COME if(WtTOGGLEIND+1<ntoggles)
#define WtSETWT_IF_MORE_TO_COME() {WtIF_MORE_TO_COME{WtSETWT_WITH_BACKUP();}}
#define WtUNDO_PREVIOUS for(int TOGGLEIND=ntoggles-2; WtTOGGLEIND>=0; WtTOGGLEIND--)
#define WtUNDO_PREVIOUS_SETWTS() {WtUNDO_PREVIOUS{WtUNDO_SETWT();}}
/* Brings together the above operations:
   For each toggle:
      Get the current edge weight.
      Calculate the change.
      Back up the current edge weight by swapping weight[i] with current edge weight.
   For each toggle:
      Undo the changes by swapping them back. */
#define WtEXEC_THROUGH_TOGGLES(subroutine){WtFOR_EACH_TOGGLE{ WtGETTOGGLEINFO(); {subroutine}; WtSETWT_IF_MORE_TO_COME();}; WtUNDO_PREVIOUS_SETWTS();}

#define SAMEDYAD(a1,b1,a2,b2) (DIRECTED? a1==a2 && b1==b2 : MIN(a1,b1)==MIN(a2,b2) && MAX(a1,b1)==MAX(a2,b2))

#define WtGETOLDWT(a,b) (SAMEDYAD(WtTAIL,WtHEAD,a,b)?WtOLDWT:WtGETWT(a,b))
#define WtGETNEWWT(a,b) (SAMEDYAD(WtTAIL,WtHEAD,a,b)?WtNEWWT:WtGETWT(a,b))
#define WtGETNEWWTOLD(a,b,old) (SAMEDYAD(WtTAIL,WtHEAD,a,b)?WtNEWWT:(old))


/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtC_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate)
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, double *weights, WtModelTerm *mtp, WtNetwork *nwp)
#define WtI_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtU_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate)
#define WtF_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtW_CHANGESTAT_FN(a) SEXP (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtX_CHANGESTAT_FN(a) void (a) (unsigned int type, void *data, WtModelTerm *mtp, WtNetwork *nwp)
#define WtZ_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp, Rboolean skip_s)

/* This macro wraps two calls to an s_??? function with toggles
   between them. */
#define WtD_FROM_S							\
  {									\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));	\
    /* Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES. */	\
    WtFOR_EACH_TOGGLE{ WtGETTOGGLEINFO(); WtSETWT_WITH_BACKUP(); }		\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)			\
      mtp->dstats[i] -= mtp->statcache[i];				\
    WtFOR_EACH_TOGGLE{ WtUNDO_SETWT(); }					\
  }

/* This macro constructs a function that wraps D_FROM_S. */
#define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(a) WtD_FROM_S

#ifdef STUBFILE
#define STRICT_Wt_HEADERS
#endif

/* If STRICT_Wt_HEADERS is not set, give the terms more generic names. */
#ifndef STRICT_Wt_HEADERS

#define IS_OUTEDGE WtIS_OUTEDGE
#define IS_INEDGE WtIS_INEDGE
#define IS_UNDIRECTED_EDGE WtIS_UNDIRECTED_EDGE
#define MIN_OUTEDGE WtMIN_OUTEDGE
#define MIN_INEDGE WtMIN_INEDGE
#define NEXT_OUTEDGE WtNEXT_OUTEDGE
#define NEXT_INEDGE WtNEXT_INEDGE
#define NEXT_OUTEDGE_PRE WtNEXT_OUTEDGE_PRE
#define NEXT_INEDGE_PRE WtNEXT_INEDGE_PRE
#define STEP_THROUGH_OUTEDGES WtSTEP_THROUGH_OUTEDGES
#define STEP_THROUGH_INEDGES WtSTEP_THROUGH_INEDGES
#define STEP_THROUGH_OUTEDGES_PRE WtSTEP_THROUGH_OUTEDGES_PRE
#define STEP_THROUGH_INEDGES_PRE WtSTEP_THROUGH_INEDGES_PRE
#define STEP_THROUGH_OUTEDGES_DECL WtSTEP_THROUGH_OUTEDGES_DECL
#define STEP_THROUGH_INEDGES_DECL WtSTEP_THROUGH_INEDGES_DECL
#define STEP_THROUGH_OUTEDGES_PRE_DECL WtSTEP_THROUGH_OUTEDGES_PRE_DECL
#define STEP_THROUGH_INEDGES_PRE_DECL WtSTEP_THROUGH_INEDGES_PRE_DECL
#define EXEC_THROUGH_OUTEDGES WtEXEC_THROUGH_OUTEDGES
#define EXEC_THROUGH_INEDGES WtEXEC_THROUGH_INEDGES
#define EXEC_THROUGH_EDGES WtEXEC_THROUGH_EDGES
#define EXEC_THROUGH_OUTEDGES_PRE WtEXEC_THROUGH_OUTEDGES_PRE
#define EXEC_THROUGH_INEDGES_PRE WtEXEC_THROUGH_INEDGES_PRE
#define EXEC_THROUGH_EDGES_PRE WtEXEC_THROUGH_EDGES_PRE
#define EXEC_THROUGH_FOUTEDGES WtEXEC_THROUGH_FOUTEDGES
#define EXEC_THROUGH_FINEDGES WtEXEC_THROUGH_FINEDGES
#define EXEC_THROUGH_FOUTEDGES_PRE WtEXEC_THROUGH_FOUTEDGES_PRE
#define EXEC_THROUGH_FINEDGES_PRE WtEXEC_THROUGH_FINEDGES_PRE
#define EXEC_THROUGH_NET_EDGES WtEXEC_THROUGH_NET_EDGES
#define EXEC_THROUGH_NET_EDGES_PRE WtEXEC_THROUGH_NET_EDGES_PRE

#define GETWT WtGETWT
#define SETWT WtSETWT
#define FOR_EACH_TOGGLE WtFOR_EACH_TOGGLE
#define TAIL WtTAIL
#define HEAD WtHEAD
#define NEWWT WtNEWWT
#define OLDWT WtOLDWT

#define TOGGLEIND WtTOGGLEIND

#define GETOLDTOGGLEINFO WtGETOLDTOGGLEINFO
#define GETTOGGLEINFO WtGETTOGGLEINFO
#define GETNEWTOGGLEINFO WtGETNEWTOGGLEINFO
#define SETWT_WITH_BACKUP WtSETWT_WITH_BACKUP
#define UNDO_SETWT WtUNDO_SETWT
#define IF_MORE_TO_COME WtIF_MORE_TO_COME
#define SETWT_IF_MORE_TO_COME WtSETWT_IF_MORE_TO_COME
#define UNDO_PREVIOUS WtUNDO_PREVIOUS
#define UNDO_PREVIOUS_SETWTS WtUNDO_PREVIOUS_SETWTS
#define EXEC_THROUGH_TOGGLES WtEXEC_THROUGH_TOGGLES
#define GETOLDWT WtGETOLDWT
#define GETNEWWT WtGETNEWWT
#define GETNEWWTOLD WtGETNEWWTOLD

#define D_FROM_S WtD_FROM_S

#endif // STRICT_Wt_HEADERS

#endif
