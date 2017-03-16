/*  File src/changestats.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgetree.h"
#include "changestat.h"

/********************  changestats:  A    ***********/
C_CHANGESTAT_FN(c_absdiff);
C_CHANGESTAT_FN(c_absdiffcat);
D_CHANGESTAT_FN(d_adegcor); S_CHANGESTAT_FN(s_adegcor);
C_CHANGESTAT_FN(c_altkstar);
C_CHANGESTAT_FN(c_asymmetric);
/********************  changestats:  B    ***********/
C_CHANGESTAT_FN(c_b1concurrent);
C_CHANGESTAT_FN(c_b1concurrent_by_attr);
C_CHANGESTAT_FN(c_b1factor);
C_CHANGESTAT_FN(c_b1degree);
C_CHANGESTAT_FN(c_b1degree_by_attr);
C_CHANGESTAT_FN(c_b1starmix);
C_CHANGESTAT_FN(c_b1starmixhomophily);
C_CHANGESTAT_FN(c_b1twostar);
C_CHANGESTAT_FN(c_b2concurrent);
C_CHANGESTAT_FN(c_b2concurrent_by_attr);
C_CHANGESTAT_FN(c_b2degree);
C_CHANGESTAT_FN(c_b2degree_by_attr);
C_CHANGESTAT_FN(c_b2factor);
C_CHANGESTAT_FN(c_b2starmix);
C_CHANGESTAT_FN(c_b2starmixhomophily);
C_CHANGESTAT_FN(c_b2twostar);
C_CHANGESTAT_FN(c_balance);
C_CHANGESTAT_FN(c_boundeddegree);
C_CHANGESTAT_FN(c_boundedidegree);
C_CHANGESTAT_FN(c_boundedodegree);
C_CHANGESTAT_FN(c_boundedistar);
C_CHANGESTAT_FN(c_boundedkstar);
C_CHANGESTAT_FN(c_boundedostar);
C_CHANGESTAT_FN(c_boundedtriangle);

/* *** don't forget tail-> head, so this function accepts tail first, not head */

  Vertex CountTriangles (Vertex tail, Vertex head, int outcount, 
                         int incount, Network *nwp);
/********************  changestats:  C    ***********/
C_CHANGESTAT_FN(c_concurrent);
C_CHANGESTAT_FN(c_concurrent_by_attr);
C_CHANGESTAT_FN(c_ctriple);
C_CHANGESTAT_FN(c_cycle);
  void edgewise_path_recurse(Network *g, Vertex dest, 
     Vertex curnode, Vertex *availnodes, long int availcount, 
     long int curlen, double *countv, long int maxlen);

/* *** I didn't swap heads and tails here, since these already 
   seem in line with the tails->heads naming convenction*/
       
  void edgewise_cycle_census(Network *g, Vertex tail, Vertex head, 
     double *countv, long int maxlen);
/********************  changestats:  D    ***********/
C_CHANGESTAT_FN(c_degcor); S_CHANGESTAT_FN(s_degcor);
C_CHANGESTAT_FN(c_degcrossprod);
C_CHANGESTAT_FN(c_degree);
C_CHANGESTAT_FN(c_degree_by_attr);
C_CHANGESTAT_FN(c_degree_w_homophily);
C_CHANGESTAT_FN(c_degreepopularity);
C_CHANGESTAT_FN(c_density);
C_CHANGESTAT_FN(c_diff);
C_CHANGESTAT_FN(c_dsp);
C_CHANGESTAT_FN(c_dyadcov);
/********************  changestats:  E    ***********/
C_CHANGESTAT_FN(c_edgecov);
C_CHANGESTAT_FN(c_edges);S_CHANGESTAT_FN(s_edges);
C_CHANGESTAT_FN(c_esp);
/********************  changestats:  F    ***********/
/********************  changestats:  G    ***********/
C_CHANGESTAT_FN(c_gwb1degree);
C_CHANGESTAT_FN(c_gwb1degree_by_attr);
C_CHANGESTAT_FN(c_gwdegree);
C_CHANGESTAT_FN(c_gwdegree_by_attr);
C_CHANGESTAT_FN(c_gwdsp);
C_CHANGESTAT_FN(c_gwb2degree);
C_CHANGESTAT_FN(c_gwb2degree_by_attr);
C_CHANGESTAT_FN(c_gwesp);
C_CHANGESTAT_FN(c_gwidegree);
C_CHANGESTAT_FN(c_gwidegree_by_attr);
C_CHANGESTAT_FN(c_gwnsp);
C_CHANGESTAT_FN(c_gwodegree);
C_CHANGESTAT_FN(c_gwodegree_by_attr);
C_CHANGESTAT_FN(c_gwtdsp);
C_CHANGESTAT_FN(c_gwtesp);
C_CHANGESTAT_FN(c_gwtnsp);
/********************  changestats:   H    ***********/
C_CHANGESTAT_FN(c_hamming);
C_CHANGESTAT_FN(c_hammingmix);
/********************  changestats:   I    ***********/
C_CHANGESTAT_FN(c_idegree);
C_CHANGESTAT_FN(c_idegree_by_attr);
C_CHANGESTAT_FN(c_idegree_w_homophily);
C_CHANGESTAT_FN(c_idegreepopularity);
C_CHANGESTAT_FN(c_intransitive);
C_CHANGESTAT_FN(c_isolates);
S_CHANGESTAT_FN(s_isolates);
C_CHANGESTAT_FN(c_istar);
/********************  changestats:   K    ***********/
C_CHANGESTAT_FN(c_kstar);
/********************  changestats:   L    ***********/
C_CHANGESTAT_FN(c_localtriangle);
/********************  changestats:   M    ***********/
C_CHANGESTAT_FN(c_m2star);
C_CHANGESTAT_FN(c_meandeg);
C_CHANGESTAT_FN(c_mix);
C_CHANGESTAT_FN(c_mutual);
C_CHANGESTAT_FN(c_mutual_by_attr);
/********************  changestats:   N    ***********/                       
C_CHANGESTAT_FN(c_nearsimmelian);
C_CHANGESTAT_FN(c_nodecov);
C_CHANGESTAT_FN(c_nodefactor);
C_CHANGESTAT_FN(c_nodeicov);
C_CHANGESTAT_FN(c_nodeifactor);
C_CHANGESTAT_FN(c_nodematch);
C_CHANGESTAT_FN(c_nodemix);
C_CHANGESTAT_FN(c_nodeocov);
C_CHANGESTAT_FN(c_nodeofactor);
C_CHANGESTAT_FN(c_nsp);
/********************  changestats:   O    ***********/
C_CHANGESTAT_FN(c_odegree);
C_CHANGESTAT_FN(c_odegree_by_attr);
C_CHANGESTAT_FN(c_odegree_w_homophily);
C_CHANGESTAT_FN(c_opentriad);
C_CHANGESTAT_FN(c_ostar);
C_CHANGESTAT_FN(c_odegreepopularity);
/********************  changestats:   P    ***********/
D_CHANGESTAT_FN(d_pdegcor); S_CHANGESTAT_FN(s_pdegcor);
/********************  changestats:   R    ***********/
D_CHANGESTAT_FN(d_rdegcor); S_CHANGESTAT_FN(s_rdegcor);
C_CHANGESTAT_FN(c_receiver);
/********************  changestats:   S    ***********/
C_CHANGESTAT_FN(c_sender);
C_CHANGESTAT_FN(c_simmelian);
C_CHANGESTAT_FN(c_simmelianties);
C_CHANGESTAT_FN(c_smalldiff);
C_CHANGESTAT_FN(c_sociality);
/********************  changestats:   T    ***********/
C_CHANGESTAT_FN(c_tdsp);
C_CHANGESTAT_FN(c_tesp);
C_CHANGESTAT_FN(c_threetrail);
C_CHANGESTAT_FN(c_tnsp);
C_CHANGESTAT_FN(c_transitive);
C_CHANGESTAT_FN(c_transitiveties); S_CHANGESTAT_FN(s_transitiveties);
C_CHANGESTAT_FN(c_triadcensus);
C_CHANGESTAT_FN(c_triangle);
C_CHANGESTAT_FN(c_tripercent);
C_CHANGESTAT_FN(c_ttriple);
C_CHANGESTAT_FN(c_transitiveties);



              
#endif
