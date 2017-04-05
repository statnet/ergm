#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgeTree.h"

typedef struct ModelTermstruct {
	void (*func)(int, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
	double *attrib; /* Ptr to vector of covariates (if necessary; generally unused) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} ModelTerm;

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
#define CHANGESTAT_FN(a) void (a) (int ntoggles, Vertex *heads, Vertex *tail, ModelTerm *mtp, Network *nwp);
/********************  changestats:  A    ***********/
CHANGESTAT_FN(d_absdiff)         
CHANGESTAT_FN(d_absdiffcat)
CHANGESTAT_FN(d_b1concurrent)
CHANGESTAT_FN(d_b1concurrent_by_attr)
CHANGESTAT_FN(d_b1factor)
CHANGESTAT_FN(d_b1degree)
CHANGESTAT_FN(d_b1degree_by_attr)
CHANGESTAT_FN(d_altkstar)
CHANGESTAT_FN(d_asymmetric)
/********************  changestats:  B    ***********/
CHANGESTAT_FN(d_balance)
CHANGESTAT_FN(d_boundeddegree)
CHANGESTAT_FN(d_boundedidegree)
CHANGESTAT_FN(d_boundedodegree)
CHANGESTAT_FN(d_boundedistar)
  double my_choose(double n, int r);
CHANGESTAT_FN(d_boundedkstar)
CHANGESTAT_FN(d_boundedostar)
CHANGESTAT_FN(d_boundedtriangle)
  Vertex CountTriangles (Vertex h, Vertex t, int outcount, 
                         int incount, Network *nwp);
/********************  changestats:  C    ***********/
CHANGESTAT_FN(d_concurrent) 
CHANGESTAT_FN(d_concurrent_by_attr)
CHANGESTAT_FN(d_ctriple)
CHANGESTAT_FN(d_cycle)
  void edgewise_path_recurse(Network *g, Vertex dest, 
     Vertex curnode, Vertex *availnodes, long int availcount, 
     long int curlen, double *countv, long int maxlen, 
     int directed);
  void edgewise_cycle_census(Network *g, Vertex t, Vertex h, 
     double *countv, long int maxlen, int directed);
/********************  changestats:  D    ***********/
CHANGESTAT_FN(d_degree)
CHANGESTAT_FN(d_degree_by_attr)
CHANGESTAT_FN(d_degree_w_homophily)
CHANGESTAT_FN(d_density)
CHANGESTAT_FN(d_dsp)
CHANGESTAT_FN(d_dyadcov)
/********************  changestats:  E    ***********/
CHANGESTAT_FN(d_b2concurrent)
CHANGESTAT_FN(d_b2concurrent_by_attr)
CHANGESTAT_FN(d_b2degree)
CHANGESTAT_FN(d_b2degree_by_attr)
CHANGESTAT_FN(d_edgecov)
CHANGESTAT_FN(d_edges)
CHANGESTAT_FN(d_esp)          
CHANGESTAT_FN(d_b2factor)
/********************  changestats:  F    ***********/
/********************  changestats:  G    ***********/
CHANGESTAT_FN(d_gwb1degree)
CHANGESTAT_FN(d_gwb1degree_by_attr)
CHANGESTAT_FN(d_gwdegree)
CHANGESTAT_FN(d_gwdegree_by_attr)
CHANGESTAT_FN(d_gwdsp)
CHANGESTAT_FN(d_gwb2degree)
CHANGESTAT_FN(d_gwb2degree_by_attr)
CHANGESTAT_FN(d_gwesp)
CHANGESTAT_FN(d_gwidegree)
CHANGESTAT_FN(d_gwidegree_by_attr)
CHANGESTAT_FN(d_gwodegree)
CHANGESTAT_FN(d_gwodegree_by_attr)
CHANGESTAT_FN(d_gwtdsp)
CHANGESTAT_FN(d_gwtesp)
/********************  changestats:   H    ***********/
CHANGESTAT_FN(d_hamming)
CHANGESTAT_FN(d_hamming_weighted)
CHANGESTAT_FN(d_hammingmix_constant)
CHANGESTAT_FN(d_hammingmix)
/********************  changestats:   I    ***********/
CHANGESTAT_FN(d_idegree)
CHANGESTAT_FN(d_idegree_by_attr)
CHANGESTAT_FN(d_idegree_w_homophily)
CHANGESTAT_FN(d_intransitive)
CHANGESTAT_FN(d_isolates)
CHANGESTAT_FN(d_istar)
/********************  changestats:   K    ***********/
CHANGESTAT_FN(d_kstar)
/********************  changestats:   L    ***********/
CHANGESTAT_FN(d_localtriangle)
/********************  changestats:   M    ***********/
CHANGESTAT_FN(d_m2star)
CHANGESTAT_FN(d_meandeg)
CHANGESTAT_FN(d_mix)
CHANGESTAT_FN(d_mutual)
/********************  changestats:   N    ***********/                       
CHANGESTAT_FN(d_nearsimmelian)
CHANGESTAT_FN(d_nodecov)
CHANGESTAT_FN(d_nodefactor)
CHANGESTAT_FN(d_nodeicov)
CHANGESTAT_FN(d_nodeifactor)
CHANGESTAT_FN(d_nodematch)
CHANGESTAT_FN(d_nodemix)
CHANGESTAT_FN(d_nodeocov)
CHANGESTAT_FN(d_nodeofactor)
/********************  changestats:   O    ***********/
CHANGESTAT_FN(d_odegree)
CHANGESTAT_FN(d_odegree_by_attr)
CHANGESTAT_FN(d_odegree_w_homophily)
CHANGESTAT_FN(d_ostar)
/********************  changestats:   R    ***********/
CHANGESTAT_FN(d_receiver)
/********************  changestats:   S    ***********/
CHANGESTAT_FN(d_sender)
CHANGESTAT_FN(d_simmelian)
CHANGESTAT_FN(d_simmelianties)
CHANGESTAT_FN(d_smalldiff)
CHANGESTAT_FN(d_sociality)
/********************  changestats:   T    ***********/
CHANGESTAT_FN(d_tdsp)
CHANGESTAT_FN(d_tesp)
CHANGESTAT_FN(d_transitive)
CHANGESTAT_FN(d_triadcensus)
CHANGESTAT_FN(d_triangle)
CHANGESTAT_FN(d_tripercent)
CHANGESTAT_FN(d_ttriple)

              
#endif
