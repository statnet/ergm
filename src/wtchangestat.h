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

#define ZERO_ALL_CHANGESTATS(a) for((a)=0; (a)<N_CHANGE_STATS; (a)++) CHANGE_STAT[(a)]=0.0
#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
// The idea here is to essentially swap the contents of the proposed
// weights with the current weights, and then swap them back when
// done.
#define GETOLDWT(a) oldwt=GETWT(heads[(a)],tails[(a)])
#define OLDWT oldwt;
#define SETWT_IF_MORE_TO_COME(a) {if((a)+1<ntoggles) {SETWT(heads[(a)],tails[(a)],weights[(a)]); weights[(a)]=oldwt;}}
#define UNDO_PREVIOUS_SETWTS(a) (a)--; while(--(a)>=0){oldwt=GETWT(heads[(a)],tails[(a)]); SETWT(heads[(a)],tails[(a)],weights[(a)]); weights[(a)]=oldwt; }

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define WtD_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, double *weights, double *oldweights,  WtModelTerm *mtp, WtNetwork *nwp)
#define WtT_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define WtS_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)              
#endif
