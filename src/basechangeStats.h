#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgeTree.h"

double my_choose(double n, int r);

#define CHOOSE(n,r) ((n)<(r) ? (0) : (my_choose((double)(n),(int)(r)))) 

/* Node and dyad covariates are now passed as part of inputparams.  However,
   attrib can still be set to point to the start of these attributes if
   you want; see comments in InitErgm.r          Dave H  12/17/2003 */
typedef struct ModelTermstruct {
	void (*func)(int, Vertex*, Vertex*, struct ModelTermstruct*, Network*);
	double *attrib; /* Ptr to vector of covariates (if necessary) */
	int nstats;   /* Number of change statistics to be returned */
	double *dstats; /* ptr to change statistics returned */
	int ninputparams; /* Number of input parameters passed to function */
	double *inputparams; /* ptr to input parameters passed */
} ModelTerm;

void d_absdiff (int ntoggles, Vertex *heads, Vertex *tail, 
                ModelTerm *mtp, Network *nwp);
void d_boundedistar (int ntoggles, Vertex *heads, Vertex *tails, 
		     ModelTerm *mtp, Network *nwp);
void d_boundedostar (int ntoggles, Vertex *heads, Vertex *tails, 
		     ModelTerm *mtp, Network *nwp);
void d_boundedkstar (int ntoggles, Vertex *heads, Vertex *tails, 
		     ModelTerm *mtp, Network *nwp);
void d_boundedtriangle (int ntoggles, Vertex *heads, Vertex *tails, 
			ModelTerm *mtp, Network *nwp);
void d_concurrent (int ntoggles, Vertex *heads, Vertex *tails, 
			ModelTerm *mtp, Network *nwp);
void d_concurrent_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
			ModelTerm *mtp, Network *nwp);
void d_ctriad (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_degree (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_degree_w_homophily (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_degree_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	        ModelTerm *mtp, Network *nwp);
void d_degreep (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_degreep_w_homophily (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_degreep_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	        ModelTerm *mtp, Network *nwp);
void d_spartner (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp);
void d_wdegree (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_wspartner (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp);
void d_geodegree (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp);
void d_gwdegreealpha (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_altkstar (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_altistar (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_altostar (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_gwdegreelambda (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_gwdegree (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_gwdegree_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_gwidegree (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_gwidegree_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_gwodegree (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_gwodegree_by_attr (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_gwd (int ntoggles, Vertex *heads, Vertex *tails, 
	    ModelTerm *mtp, Network *nwp);
void d_geospartner (int ntoggles, Vertex *heads, Vertex *tails, 
		    ModelTerm *mtp, Network *nwp);
void d_geotwopath (int ntoggles, Vertex *heads, Vertex *tails, 
		   ModelTerm *mtp, Network *nwp);
void d_dyadcov (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_heideriandynamic (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_simmeliandynamic (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_edgecov (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_edges (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_meandeg (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_density (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_idegree (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_istar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_kstar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_nodematch (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp);
void d_receiver (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp);
void d_sender (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_sociality (int ntoggles, Vertex *heads, Vertex *tails, 
		   ModelTerm *mtp, Network *nwp);
void d_nearsimmian (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_simmian (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_mutual (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_asymmetric (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_nodecov (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_sendercov (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp);
void d_receivercov (int ntoggles, Vertex *heads, Vertex *tails, 
		    ModelTerm *mtp, Network *nwp);
void d_hamming (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_duration (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp);
void d_factor (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_nodefactor (int ntoggles, Vertex *heads, Vertex *tails, 
		   ModelTerm *mtp, Network *nwp);
void d_nodeifactor (int ntoggles, Vertex *heads, Vertex *tails, 
		   ModelTerm *mtp, Network *nwp);
void d_nodeofactor (int ntoggles, Vertex *heads, Vertex *tails, 
		   ModelTerm *mtp, Network *nwp);
void d_odegree (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_ostar (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_smalldiff (int ntoggles, Vertex *heads, Vertex *tails, 
		  ModelTerm *mtp, Network *nwp);
void d_localtriangle (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp);
void d_triangle (int ntoggles, Vertex *heads, Vertex *tails, 
		 ModelTerm *mtp, Network *nwp);
void d_tricorr (int ntoggles, Vertex *heads, Vertex *tails, 
		ModelTerm *mtp, Network *nwp);
void d_ttriad (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_kappa (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
void d_boundeddegree (int ntoggles, Vertex *heads, Vertex *tails, 
		      ModelTerm *mtp, Network *nwp);
void d_boundedodegree (int ntoggles, Vertex *heads, Vertex *tails, 
		       ModelTerm *mtp, Network *nwp);
void d_boundedidegree (int ntoggles, Vertex *heads, Vertex *tails, 
		       ModelTerm *mtp, Network *nwp);
void d_hiertriad (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
void d_hiertriaddegree (int ntoggles, Vertex *heads, Vertex *tails, 
	       ModelTerm *mtp, Network *nwp);
double numposthree (Vertex t, Network *nwp);
Vertex CountTriangles (Vertex h, Vertex t, int outcount, int incount, Network *nwp);

void edgewise_path_recurse(Network *g, Vertex dest, Vertex curnode, Vertex *availnodes, long int availcount, long int curlen, double *countv, long int maxlen, int directed);

void edgewise_cycle_census(Network *g, Vertex t, Vertex h, double *countv, long int maxlen, int directed);

void d_cycle (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);

void d_berninhom (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
              
void d_spatial (int ntoggles, Vertex *heads, Vertex *tails, 
	      ModelTerm *mtp, Network *nwp);
              
void d_nodemix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_mix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_hammingmix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_hammingfixmix (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_icvar (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
void d_idc (int ntoggles, Vertex *heads, Vertex *tails,
              ModelTerm *mtp, Network *nwp);
#endif
