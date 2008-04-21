#ifndef CHANGESTATS_H
#define CHANGESTATS_H

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
#define D_CHANGESTAT_FN(a) void (a) (Edge ntoggles, Vertex *heads, Vertex *tails, ModelTerm *mtp, Network *nwp)
#define T_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
#define S_CHANGESTAT_FN(a) void (a) (ModelTerm *mtp, Network *nwp)
/********************  changestats:  A    ***********/
D_CHANGESTAT_FN(d_absdiff);
D_CHANGESTAT_FN(d_absdiffcat);
D_CHANGESTAT_FN(d_altkstar);
D_CHANGESTAT_FN(d_asymmetric);
/********************  changestats:  B    ***********/
D_CHANGESTAT_FN(d_b1concurrent);
D_CHANGESTAT_FN(d_b1concurrent_by_attr);
D_CHANGESTAT_FN(d_b1factor);
D_CHANGESTAT_FN(d_b1degree);
D_CHANGESTAT_FN(d_b1degree_by_attr);
D_CHANGESTAT_FN(d_b2concurrent);
D_CHANGESTAT_FN(d_b2concurrent_by_attr);
D_CHANGESTAT_FN(d_b2degree);
D_CHANGESTAT_FN(d_b2degree_by_attr);
D_CHANGESTAT_FN(d_b2factor);
D_CHANGESTAT_FN(d_balance);
D_CHANGESTAT_FN(d_boundeddegree);
D_CHANGESTAT_FN(d_boundedidegree);
D_CHANGESTAT_FN(d_boundedodegree);
D_CHANGESTAT_FN(d_boundedistar);
  double my_choose(double n, int r);
D_CHANGESTAT_FN(d_boundedkstar);
D_CHANGESTAT_FN(d_boundedostar);
D_CHANGESTAT_FN(d_boundedtriangle);
  Vertex CountTriangles (Vertex h, Vertex t, int outcount, 
                         int incount, Network *nwp);
/********************  changestats:  C    ***********/
D_CHANGESTAT_FN(d_concurrent);
D_CHANGESTAT_FN(d_concurrent_by_attr);
D_CHANGESTAT_FN(d_ctriple);
D_CHANGESTAT_FN(d_cycle);
  void edgewise_path_recurse(Network *g, Vertex dest, 
     Vertex curnode, Vertex *availnodes, long int availcount, 
     long int curlen, double *countv, long int maxlen);
  void edgewise_cycle_census(Network *g, Vertex t, Vertex h, 
     double *countv, long int maxlen);
/********************  changestats:  D    ***********/
D_CHANGESTAT_FN(d_degree);
D_CHANGESTAT_FN(d_degree_by_attr);
D_CHANGESTAT_FN(d_degree_w_homophily);
D_CHANGESTAT_FN(d_density);
D_CHANGESTAT_FN(d_dsp);
D_CHANGESTAT_FN(d_dyadcov);
/********************  changestats:  E    ***********/
D_CHANGESTAT_FN(d_edgecov);
D_CHANGESTAT_FN(d_edges);S_CHANGESTAT_FN(s_edges);
D_CHANGESTAT_FN(d_esp);
/********************  changestats:  F    ***********/
/********************  changestats:  G    ***********/
D_CHANGESTAT_FN(d_gwb1degree);
D_CHANGESTAT_FN(d_gwb1degree_by_attr);
D_CHANGESTAT_FN(d_gwdegree);
D_CHANGESTAT_FN(d_gwdegree_by_attr);
D_CHANGESTAT_FN(d_gwdsp);
D_CHANGESTAT_FN(d_gwb2degree);
D_CHANGESTAT_FN(d_gwb2degree_by_attr);
D_CHANGESTAT_FN(d_gwesp);
D_CHANGESTAT_FN(d_gwidegree);
D_CHANGESTAT_FN(d_gwidegree_by_attr);
D_CHANGESTAT_FN(d_gwodegree);
D_CHANGESTAT_FN(d_gwodegree_by_attr);
D_CHANGESTAT_FN(d_gwtdsp);
D_CHANGESTAT_FN(d_gwtesp);
/********************  changestats:   H    ***********/
D_CHANGESTAT_FN(d_hamming);
D_CHANGESTAT_FN(d_hamming_weighted);
D_CHANGESTAT_FN(d_hammingmix_constant);
D_CHANGESTAT_FN(d_hammingmix);
/********************  changestats:   I    ***********/
D_CHANGESTAT_FN(d_idegree);
D_CHANGESTAT_FN(d_idegree_by_attr);
D_CHANGESTAT_FN(d_idegree_w_homophily);
D_CHANGESTAT_FN(d_intransitive);
D_CHANGESTAT_FN(d_isolates);
D_CHANGESTAT_FN(d_istar);
/********************  changestats:   K    ***********/
D_CHANGESTAT_FN(d_kstar);
/********************  changestats:   L    ***********/
D_CHANGESTAT_FN(d_localtriangle);
/********************  changestats:   M    ***********/
D_CHANGESTAT_FN(d_m2star);
D_CHANGESTAT_FN(d_meandeg);
D_CHANGESTAT_FN(d_mix);
D_CHANGESTAT_FN(d_mutual);
/********************  changestats:   N    ***********/                       
D_CHANGESTAT_FN(d_nearsimmelian);
D_CHANGESTAT_FN(d_nodecov);
D_CHANGESTAT_FN(d_nodefactor);
D_CHANGESTAT_FN(d_nodeicov);
D_CHANGESTAT_FN(d_nodeifactor);
D_CHANGESTAT_FN(d_nodematch);
D_CHANGESTAT_FN(d_nodemix);
D_CHANGESTAT_FN(d_nodeocov);
D_CHANGESTAT_FN(d_nodeofactor);
/********************  changestats:   O    ***********/
D_CHANGESTAT_FN(d_odegree);
D_CHANGESTAT_FN(d_odegree_by_attr);
D_CHANGESTAT_FN(d_odegree_w_homophily);
D_CHANGESTAT_FN(d_ostar);
/********************  changestats:   R    ***********/
D_CHANGESTAT_FN(d_receiver);
/********************  changestats:   S    ***********/
D_CHANGESTAT_FN(d_sender);
D_CHANGESTAT_FN(d_simmelian);
D_CHANGESTAT_FN(d_simmelianties);
D_CHANGESTAT_FN(d_smalldiff);
D_CHANGESTAT_FN(d_sociality);
/********************  changestats:   T    ***********/
D_CHANGESTAT_FN(d_tdsp);
D_CHANGESTAT_FN(d_tesp);
D_CHANGESTAT_FN(d_transitive);
D_CHANGESTAT_FN(d_triadcensus);
D_CHANGESTAT_FN(d_triangle);
D_CHANGESTAT_FN(d_tripercent);
D_CHANGESTAT_FN(d_ttriple);

              
#endif
