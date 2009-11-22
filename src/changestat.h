#ifndef CHANGESTAT_H
#define CHANGESTAT_H

#include "edgetree.h"

typedef struct ModelTermstruct {
	void (*d_func)(Edge, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
  	void (*s_func)(struct ModelTermstruct*, Network*);
        void (*t_func)(struct ModelTermstruct*, Network*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} ModelTerm;


/****************************************************
 Macros to make life easier                         *
 Note:  These things still need to be documented    */ 
#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 

#define IS_OUTEDGE(a,b) (EdgetreeSearch((a),(b),nwp->outedges)!=0?1:0)
#define IS_INEDGE(a,b) (EdgetreeSearch((a),(b),nwp->inedges)!=0?1:0)
#define IS_UNDIRECTED_EDGE(a,b) IS_OUTEDGE(MIN(a,b), MAX(a,b))
#define MIN_OUTEDGE(a) (EdgetreeMinimum(nwp->outedges, (a)))
#define MIN_INEDGE(a) (EdgetreeMinimum(nwp->inedges, (a)))
#define NEXT_OUTEDGE(e) (EdgetreeSuccessor(nwp->outedges,(e)))
#define NEXT_INEDGE(e) (EdgetreeSuccessor(nwp->inedges,(e)))
#define OUTVAL(e) (nwp->outedges[(e)].value)
#define INVAL(e) (nwp->inedges[(e)].value)
#define TOGGLE(a,b) (ToggleEdge((a),(b),nwp));
#define TOGGLE_DISCORD(a,b) (ToggleEdge((a),(b),nwp+1));

#define STEP_THROUGH_OUTEDGES(a,e,v) for((e)=MIN_OUTEDGE(a);((v)=OUTVAL(e))!=0;(e)=NEXT_OUTEDGE(e))
#define STEP_THROUGH_INEDGES(a,e,v) for((e)=MIN_INEDGE(a);((v)=INVAL(e))!=0;(e)=NEXT_INEDGE(e))

#define N_NODES (nwp->nnodes)
#define N_DYADS (nwp->directed_flag?(nnodes*(nnodes-1)):nwp->bipartite?nwp->bipartite*(nnodes-nwp->bipartite):(nnodes*(nnodes-1)/2))
#define OUT_DEG (nwp->outdegree)
#define IN_DEG (nwp->indegree)
#define DIRECTED (nwp->directed_flag)
#define BIPARTITE (nwp->bipartite)
#define N_EDGES (nwp->nedges)
#define NEXT_INEDGE_NUM (nwp->next_inedge)
#define NEXT_OUTEDGE_NUM (nwp->next_outedge)

#define N_CHANGE_STATS (mtp->nstats)
#define INPUT_PARAM (mtp->inputparams)
#define CHANGE_STAT (mtp->dstats)
#define INPUT_ATTRIB (mtp->attrib)
#define N_INPUT_PARAMS (mtp->ninputparams)

#define ZERO_ALL_CHANGESTATS(a) for((a)=0; (a)<N_CHANGE_STATS; (a)++) CHANGE_STAT[(a)]=0.0
#define FOR_EACH_TOGGLE(a) for((a)=0; (a)<ntoggles; (a)++)
#define TOGGLE_IF_MORE_TO_COME(a) if((a)+1<ntoggles) TOGGLE(heads[(a)],tails[(a)])
#define TOGGLE_DISCORD_IF_MORE_TO_COME(a) if((a)+1<ntoggles) TOGGLE_DISCORD(heads[(a)],tails[(a)])
#define UNDO_PREVIOUS_TOGGLES(a) (a)--; while(--(a)>=0) TOGGLE(heads[(a)],tails[(a)])
#define UNDO_PREVIOUS_DISCORD_TOGGLES(a) (a)--; while(--(a)>=0) {TOGGLE(heads[(a)],tails[(a)]); TOGGLE_DISCORD(heads[(a)],tails[(a)])}

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, ModelTerm *mtp, Network *nwp)
/* NB:  CHANGESTAT_FN is now deprecated (replaced by D_CHANGESTAT_FN) */
#define D_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, ModelTerm *mtp, Network *nwp)
#define T_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define S_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)              
#endif
