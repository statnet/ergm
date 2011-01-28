#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgetree.h"
#include "changestat.h"

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
D_CHANGESTAT_FN(d_b1starmix);
D_CHANGESTAT_FN(d_b1starmixhomophily);
D_CHANGESTAT_FN(d_b1twostar);
D_CHANGESTAT_FN(d_b2concurrent);
D_CHANGESTAT_FN(d_b2concurrent_by_attr);
D_CHANGESTAT_FN(d_b2degree);
D_CHANGESTAT_FN(d_b2degree_by_attr);
D_CHANGESTAT_FN(d_b2factor);
D_CHANGESTAT_FN(d_b2starmix);
D_CHANGESTAT_FN(d_b2starmixhomophily);
D_CHANGESTAT_FN(d_b2twostar);
D_CHANGESTAT_FN(d_balance);
D_CHANGESTAT_FN(d_boundeddegree);
D_CHANGESTAT_FN(d_boundedidegree);
D_CHANGESTAT_FN(d_boundedodegree);
D_CHANGESTAT_FN(d_boundedistar);
  double my_choose(double n, int r);
D_CHANGESTAT_FN(d_boundedkstar);
D_CHANGESTAT_FN(d_boundedostar);
D_CHANGESTAT_FN(d_boundedtriangle);

/* *** don't forget tail-> head, so this function accepts tail first, not head */

  Vertex CountTriangles (Vertex tail, Vertex head, int outcount, 
                         int incount, Network *nwp);
/********************  changestats:  C    ***********/
D_CHANGESTAT_FN(d_concurrent);
D_CHANGESTAT_FN(d_concurrent_by_attr);
D_CHANGESTAT_FN(d_ctriple);
D_CHANGESTAT_FN(d_cycle);
  void edgewise_path_recurse(Network *g, Vertex dest, 
     Vertex curnode, Vertex *availnodes, long int availcount, 
     long int curlen, double *countv, long int maxlen);

/* *** I didn't swap heads and tails here, since these already 
   seem in line with the tails->heads naming convenction*/
       
  void edgewise_cycle_census(Network *g, Vertex tail, Vertex head, 
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
D_CHANGESTAT_FN(d_gwnsp);
D_CHANGESTAT_FN(d_gwodegree);
D_CHANGESTAT_FN(d_gwodegree_by_attr);
D_CHANGESTAT_FN(d_gwtdsp);
D_CHANGESTAT_FN(d_gwtesp);
D_CHANGESTAT_FN(d_gwtnsp);
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
D_CHANGESTAT_FN(d_nsp);
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
D_CHANGESTAT_FN(d_threepath);
D_CHANGESTAT_FN(d_tnsp);
D_CHANGESTAT_FN(d_transitive);
D_CHANGESTAT_FN(d_transitiveties); S_CHANGESTAT_FN(s_transitiveties);
D_CHANGESTAT_FN(d_triadcensus);
D_CHANGESTAT_FN(d_triangle);
D_CHANGESTAT_FN(d_tripercent);
D_CHANGESTAT_FN(d_ttriple);

              
#endif
