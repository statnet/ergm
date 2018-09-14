/*  File src/ergm_changestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2017 Statnet Commons
 */
#ifndef _ERGM_CHANGESTAT_H_
#define _ERGM_CHANGESTAT_H_

#include "ergm_edgetree.h"
#include "ergm_edgelist.h"

typedef struct ModelTermstruct {
  void (*c_func)(Vertex, Vertex, struct ModelTermstruct*, Network*);
  void (*d_func)(Edge, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
  void (*i_func)(struct ModelTermstruct*, Network*);
  void (*u_func)(Vertex, Vertex, struct ModelTermstruct*, Network*);
  void (*f_func)(struct ModelTermstruct*, Network*);
  void (*s_func)(struct ModelTermstruct*, Network*);
  double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
  int nstats;   /* Number of change statistics to be returned */
  unsigned int statspos; /* Position of this term's stats in the workspace vector. */ 
  double *dstats; /* ptr to change statistics returned */
  int ninputparams; /* Number of input parameters passed to function */
  double *inputparams; /* ptr to input parameters passed */
  double *statcache; /* vector of the same length as dstats */
  void *storage; /* optional space for persistent storage */
  void **aux_storage; /* optional space for persistent public (auxiliary) storage */
} ModelTerm;

#include "ergm_changestat_common.do_not_include_directly.h"

/****************************************************
 Macros to make life easier when writing C code for change statistics:  */

/* return number of tail and head node in the directed node pair
   tail -> head of the selected toggle */
#define TAIL(a) (tails[(a)])
#define HEAD(a) (heads[(a)])

/* tell whether a particular edge exists */
#define IS_OUTEDGE(a,b) (EdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (EdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))

/* Return the Edge number of the smallest-labelled neighbor of the node 
   labelled "a".  Or, return the Edge number of the next-largest neighbor 
   starting from the pointer "e", which points to a node in an edgetree. 
   Mostly, these are utility macros used by the STEP_THROUGH_OUTEDGES 
   and STEP_THROUGH_INEDGES macros. */
#define MIN_OUTEDGE(a) (EdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (EdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (EdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (EdgetreeSuccessor(nwp->inedges,(e)))
/* As NEXT_*EDGE, but visits the parent nodes before the child
   nodes. */
#define NEXT_OUTEDGE_PRE(e) (EdgetreePreSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE_PRE(e) (EdgetreePreSuccessor(nwp->inedges,(e)))

/* Return each of the out-neighbors or in-neighbors, one at a time,
   of node a.  At each iteration of the loop, the variable v gives the node 
   number of the corresponding neighbor.  The e variable, which should be
   initialized as type Edge, is merely the looping variable. */
#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))

/* As STEP_THROUGH_*EDGES, but visit the parent nodes before the child
   nodes. This is useful for "copying" an edgetree. */
#define STEP_THROUGH_OUTEDGES_PRE(a,e,v) for((e)=(a); ((v)=OUTVAL(e))!=0;NEXT_OUTEDGE_PRE(e))
#define STEP_THROUGH_INEDGES_PRE(a,e,v) for((e)=(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE_PRE(e))


// These are "declaring" versions of the above, optimized for use in EXEC_TROUGH_*EDGES macros.
#define STEP_THROUGH_OUTEDGES_DECL(a,e,v) Vertex v; for(Edge e=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;e=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES_DECL(a,e,v) Vertex v; for(Edge e=MIN_INEDGE(a);((v)=INVAL(e))!=0;e=NEXT_INEDGE(e))
#define STEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) Vertex v; for(Edge e=(a);((v)=OUTVAL(e))!=0;e=NEXT_OUTEDGE_PRE(e))
#define STEP_THROUGH_INEDGES_PRE_DECL(a,e,v) Vertex v; for(Edge e=(a);((v)=INVAL(e))!=0;e=NEXT_INEDGE_PRE(e))

/* Instead of stepping through execute "subroutine" for each neighbor
   automatically adapting to undirected networks. */
/* NOTE: For some reason, GCC complains when it encounters multiple
   declaration in a single statement inside the subroutine, so
   double v1, v2;
   might fail, while
   double v1;
   double v2;
   works.*/
#define EXEC_THROUGH_OUTEDGES(a,e,v,subroutine) {if(DIRECTED){ STEP_THROUGH_OUTEDGES_DECL(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,subroutine) }}
#define EXEC_THROUGH_INEDGES(a,e,v,subroutine) {if(DIRECTED){ STEP_THROUGH_INEDGES_DECL(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES(a,e,v,subroutine) }}
#define EXEC_THROUGH_EDGES(a,e,v,subroutine) { {STEP_THROUGH_OUTEDGES_DECL(a,e,v) {subroutine}};  {STEP_THROUGH_INEDGES_DECL(a,e,v) {subroutine}}; }
#define EXEC_THROUGH_OUTEDGES_PRE(a,e,v,subroutine) {if(DIRECTED){ STEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES_PRE(a,e,v,subroutine) }}
#define EXEC_THROUGH_INEDGES_PRE(a,e,v,subroutine) {if(DIRECTED){ STEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {subroutine} } else { EXEC_THROUGH_EDGES_PRE(a,e,v,subroutine) }}
#define EXEC_THROUGH_EDGES_PRE(a,e,v,subroutine) { {STEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {subroutine}};  {STEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {subroutine}}; }

/* Non-adaptive versions of the above. (I.e. ForceOUT/INEDGES.) */
#define EXEC_THROUGH_FOUTEDGES(a,e,v,subroutine) STEP_THROUGH_OUTEDGES_DECL(a,e,v) {subroutine}
#define EXEC_THROUGH_FINEDGES(a,e,v,subroutine) STEP_THROUGH_INEDGES_DECL(a,e,v) {subroutine}
#define EXEC_THROUGH_FOUTEDGES_PRE(a,e,v,subroutine) STEP_THROUGH_OUTEDGES_PRE_DECL(a,e,v) {subroutine}
#define EXEC_THROUGH_FINEDGES_PRE(a,e,v,subroutine) STEP_THROUGH_INEDGES_PRE_DECL(a,e,v) {subroutine}

/* Exectute through all edges (nonzero values) in the network. */
#define EXEC_THROUGH_NET_EDGES(a,b,e,subroutine) {for(Vertex a=1; a <= N_NODES; a++){EXEC_THROUGH_FOUTEDGES(a, e, b, {subroutine});}}
#define EXEC_THROUGH_NET_EDGES_PRE(a,b,e,subroutine) {for(Vertex a=1; a <= N_NODES; a++){EXEC_THROUGH_FOUTEDGES_PRE(a, e, b, {subroutine});}}

/* Change the status of the (a,b) edge:  Add it if it's absent, or 
   delete it if it's present. */
#define TOGGLE(a,b) (ToggleEdge((a),(b),nwp));
#define TOGGLE_DISCORD(a,b) (ToggleEdge((a),(b),nwp+1));

/* Get and set the value (0 or 1) of the (a,b) edge. */
#define GETWT(a,b) (GetEdge(a,b,nwp))
#define SETWT(a,b,w) (SetEdge(a,b,w,nwp))

/* *** don't forget tail-> head, so these functions now toggle (tails, heads), instead of (heads, tails) */

#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
#define IF_MORE_TO_COME(a) if((a)+1<ntoggles)
#define TOGGLE_IF_MORE_TO_COME(a) IF_MORE_TO_COME(a){TOGGLE(tails[(a)],heads[(a)])}
#define TOGGLE_DISCORD_IF_MORE_TO_COME(a) IF_MORE_TO_COME(a){TOGGLE_DISCORD(tails[(a)],heads[(a)])}
#define UNDO_PREVIOUS(a) (a)--; while(--(a)>=0)
#define UNDO_PREVIOUS_TOGGLES(a)  UNDO_PREVIOUS(a){TOGGLE(tails[(a)],heads[(a)])}
#define UNDO_PREVIOUS_DISCORD_TOGGLES(a) UNDO_PREVIOUS(a){TOGGLE(tails[(a)],heads[(a)]); TOGGLE_DISCORD(tails[(a)],heads[(a)])}

/****************************************************/
/* changestat function prototypes */

/* *** don't forget tail -> head, so this prototype now accepts tails first, not heads first */

#define CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)

/* NB:  CHANGESTAT_FN is now deprecated (replaced by D_CHANGESTAT_FN) */
#define C_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp)
#define D_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp)
#define I_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define U_CHANGESTAT_FN(a) void (a) (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp)
#define F_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define S_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)

/* This macro wraps two calls to an s_??? function with toggles
   between them. */
#define D_FROM_S							\
  {									\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    memcpy(mtp->statcache,mtp->dstats,N_CHANGE_STATS*sizeof(double));	\
    /* Note: This cannot be abstracted into EXEC_THROUGH_TOGGLES. */	\
    int j;								\
    FOR_EACH_TOGGLE(j) TOGGLE(TAIL(j),HEAD(j));				\
    (*(mtp->s_func))(mtp, nwp);  /* Call s_??? function */		\
    for(unsigned int i=0; i<N_CHANGE_STATS; i++)			\
      mtp->dstats[i] -= mtp->statcache[i];				\
    FOR_EACH_TOGGLE(j) TOGGLE(TAIL(j),HEAD(j));				\
  }

/* This macro constructs a function that wraps D_FROM_S. */
#define D_FROM_S_FN(a) D_CHANGESTAT_FN(a) D_FROM_S

#endif
