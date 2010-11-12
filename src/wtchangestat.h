#ifndef WTCHANGESTAT_H
#define WTCHANGESTAT_H

#include "wtedgetree.h"

typedef struct WtModelTermstruct {
  void (*d_func)(Edge, Vertex*, Vertex*, double*, struct WtModelTermstruct*, WtNetwork*);
  	void (*s_func)(struct WtModelTermstruct*, WtNetwork*);
        void (*t_func)(struct WtModelTermstruct*, WtNetwork*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} WtModelTerm;


/****************************************************
 Macros to make life easier                         *
 Note:  These things still need to be documented    */ 

#define IS_OUTEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (WtEdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))
#define MIN_OUTEDGE(a) (WtEdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (WtEdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (WtEdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (WtEdgetreeSuccessor(nwp->inedges,(e)))
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)
//#define TOGGLE(a,b) (WtToggleEdge((a),(b),nwp));
//#define TOGGLE_DISCORD(a,b) (WtToggleEdge((a),(b),nwp+1));

#define GETWT(h,t) (WtGetEdge(h,t,nwp))
#define SETWT(h,t,w) (WtSetEdge(h,t,w,nwp))
#define N_NODES (nwp->nnodes)
#define N_DYADS (nwp->directed_flag?(nnodes*(nnodes-1)):nwp->bipartite?nwp->bipartite*(nnodes-nwp->bipartite):(nnodes*(nnodes-1)/2))
#define DIRECTED (nwp->directed_flag)
#define BIPARTITE (nwp->bipartite)

#define N_CHANGE_STATS (mtp->nstats)
#define INPUT_PARAM (mtp->inputparams)
#define CHANGE_STAT (mtp->dstats)
#define INPUT_ATTRIB (mtp->attrib)
#define N_INPUT_PARAMS (mtp->ninputparams)


// Original STEP_THROUGH_*EDGES macros
#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))

// Also execute for each edge, automatically adapting to undirected networks.
#define EXEC_THROUGH_OUTEDGES(a,e,v,subroutine) if(DIRECTED){ STEP_THROUGH_OUTEDGES(a,e,v) {subroutine} } else EXEC_THROUGH_EDGES(a,e,v,subroutine)
#define EXEC_THROUGH_INEDGES(a,e,v,subroutine) if(DIRECTED){ STEP_THROUGH_INEDGES(a,e,v) {subroutine} } else EXEC_THROUGH_EDGES(a,e,v,subroutine)
#define EXEC_THROUGH_EDGES(a,e,v,subroutine) STEP_THROUGH_OUTEDGES(a,e,v) {subroutine}  STEP_THROUGH_INEDGES(a,e,v) {subroutine} 

// Non-adaptive versions. (I.e. ForceOUT/INEDGES.)
#define EXEC_THROUGH_FOUTEDGES(a,e,v,subroutine) STEP_THROUGH_OUTEDGES(a,e,v) {subroutine}
#define EXEC_THROUGH_FINEDGES(a,e,v,subroutine) STEP_THROUGH_INEDGES(a,e,v) {subroutine}



#define ZERO_ALL_CHANGESTATS(a) for((a)=0; (a)<N_CHANGE_STATS; (a)++) CHANGE_STAT[(a)]=0.0
#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
// The idea here is to essentially swap the contents of the proposed
// weights with the current weights, and then swap them back when
// done.
#define OLDWT oldwt
#define GETOLDWT(a) OLDWT=GETWT(heads[(a)],tails[(a)])

// Yes, the next two macros are, in fact, identical:
// swap(x,y) is a self-inverse.
#define SETWT_WITH_BACKUP(a) {GETOLDWT(a); SETWT(heads[(a)],tails[(a)],weights[(a)]); weights[(a)]=OLDWT;}
#define UNDO_SETWT(a) {GETOLDWT(a); SETWT(heads[(a)],tails[(a)],weights[(a)]); weights[(a)]=OLDWT;}
#define SETWT_IF_MORE_TO_COME(a) {if((a)+1<ntoggles) {SETWT_WITH_BACKUP(a)}}
#define UNDO_PREVIOUS_SETWTS(a) (a)--; while(--(a)>=0){UNDO_SETWT(a)}

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, double *weights, WtModelTerm *mtp, WtNetwork *nwp)
#define WtT_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (WtModelTerm *mtp, WtNetwork *nwp)

WtD_CHANGESTAT_FN(d_from_s);
// This could be done more efficiently (saving a function call) 
// by assigning a function pointer as follows:
 #define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(*a)=d_from_s;
// However, it looks like it might confuse the function finding routines.
// In the future, it might be a good idea to have the initialization
// code autodetect when D_ function is not found, but S_ function is, and set it properly.
//#define WtD_FROM_S_FN(a) WtD_CHANGESTAT_FN(a){ d_from_s(ntoggles, heads, tails, weights, mtp, nwp); }

#endif
