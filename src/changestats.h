/*  File src/changestats.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/********************  changestats:  A    ***********/
D_CHANGESTAT_FN(d_absdiff);
D_CHANGESTAT_FN(d_absdiffcat);
D_CHANGESTAT_FN(d_adegcor); S_CHANGESTAT_FN(s_adegcor);
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

void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode, 
     Vertex *visited, long int curlen, double *countv, long int maxlen, int semi);

void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head, 
                           double *countv, long int maxlen, int semi);
/********************  changestats:  D    ***********/
D_CHANGESTAT_FN(d_degcor); S_CHANGESTAT_FN(s_degcor);
D_CHANGESTAT_FN(d_degcrossprod);
D_CHANGESTAT_FN(d_degree);
D_CHANGESTAT_FN(d_degree_by_attr);
D_CHANGESTAT_FN(d_degree_w_homophily);
D_CHANGESTAT_FN(d_degreepopularity);
D_CHANGESTAT_FN(d_density);
D_CHANGESTAT_FN(d_diff);
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
D_CHANGESTAT_FN(d_hammingmix);
/********************  changestats:   I    ***********/
D_CHANGESTAT_FN(d_idegree);
D_CHANGESTAT_FN(d_idegree_by_attr);
D_CHANGESTAT_FN(d_idegree_w_homophily);
D_CHANGESTAT_FN(d_idegreepopularity);
D_CHANGESTAT_FN(d_intransitive);
D_CHANGESTAT_FN(d_isolatededges);
D_CHANGESTAT_FN(d_isolates);
S_CHANGESTAT_FN(s_isolates);
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
D_CHANGESTAT_FN(d_mutual_by_attr);
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
D_CHANGESTAT_FN(d_opentriad);
D_CHANGESTAT_FN(d_ostar);
D_CHANGESTAT_FN(d_odegreepopularity);
/********************  changestats:   P    ***********/
D_CHANGESTAT_FN(d_pdegcor); S_CHANGESTAT_FN(s_pdegcor);
/********************  changestats:   R    ***********/
D_CHANGESTAT_FN(d_rdegcor); S_CHANGESTAT_FN(s_rdegcor);
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
D_CHANGESTAT_FN(d_threetrail);
D_CHANGESTAT_FN(d_tnsp);
D_CHANGESTAT_FN(d_transitive);
D_CHANGESTAT_FN(d_transitiveties); S_CHANGESTAT_FN(s_transitiveties);
D_CHANGESTAT_FN(d_triadcensus);
D_CHANGESTAT_FN(d_triangle);
D_CHANGESTAT_FN(d_tripercent);
D_CHANGESTAT_FN(d_ttriple);
D_CHANGESTAT_FN(d_transitiveties);



              
#endif
