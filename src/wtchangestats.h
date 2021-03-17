/*  File src/wtchangestats.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef WTCHANGESTATS_H
#define WTCHANGESTATS_H

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_storage.h"
#include "ergm_Rutil.h"


/********************  Utility macros    ***********/

// 0 = untransformed
// 1 = square root
#define TRANSFORM_DYADVAL(y, transcode) (transcode==0? y : transcode==1? sqrt(y): 0)

/********************  changestats:   A    ***********/
WtD_CHANGESTAT_FN(d_atleast);

/********************  changestats:   C    ***********/
WtD_CHANGESTAT_FN(d_cyclicalweighs); WtS_CHANGESTAT_FN(s_cyclicalweights);
WtD_CHANGESTAT_FN(d_cyclicalweights_threshold); WtS_CHANGESTAT_FN(s_cyclicalweights_threshold);

/********************  changestats:   E    ***********/
WtD_CHANGESTAT_FN(d_edgecov_nonzero); WtD_CHANGESTAT_FN(d_edgecov_sum);

/********************  changestats:   D    ***********/
WtD_CHANGESTAT_FN(d_diff_nonzero);
WtD_CHANGESTAT_FN(d_diff_sum);

/********************  changestats:   G    ***********/
WtD_CHANGESTAT_FN(d_greaterthan);

/********************  changestats:   I    ***********/
WtD_CHANGESTAT_FN(d_ininterval);

/********************  changestats:   M    ***********/
WtD_CHANGESTAT_FN(d_mutual_wt_product);
WtD_CHANGESTAT_FN(d_mutual_wt_geom_mean);
WtD_CHANGESTAT_FN(d_mutual_wt_min); 
WtD_CHANGESTAT_FN(d_mutual_wt_nabsdiff);
WtD_CHANGESTAT_FN(d_mutual_wt_threshold);

/********************  changestats:   N    ***********/
WtD_CHANGESTAT_FN(d_nodecorr);
WtD_CHANGESTAT_FN(d_nodecov_nonzero);
WtD_CHANGESTAT_FN(d_nodecov_sum);
WtD_CHANGESTAT_FN(d_nodefactor_nonzero);
WtD_CHANGESTAT_FN(d_nodefactor_sum);
WtD_CHANGESTAT_FN(d_nodeicorr);
WtD_CHANGESTAT_FN(d_nodeifactor_nonzero);
WtD_CHANGESTAT_FN(d_nodeifactor_sum);
WtD_CHANGESTAT_FN(d_nodematch_nonzero);
WtD_CHANGESTAT_FN(d_nodematch_sum);
WtI_CHANGESTAT_FN(i_nodemix_nonzero);
WtI_CHANGESTAT_FN(i_nodemix_sum);
WtC_CHANGESTAT_FN(c_nodemix_nonzero);
WtC_CHANGESTAT_FN(c_nodemix_sum);
WtF_CHANGESTAT_FN(f_nodemix_nonzero);
WtF_CHANGESTAT_FN(f_nodemix_sum);
WtD_CHANGESTAT_FN(d_nodeocorr);
WtD_CHANGESTAT_FN(d_nodeofactor_nonzero);
WtD_CHANGESTAT_FN(d_nodeofactor_sum);
WtD_CHANGESTAT_FN(d_nonzero);

/********************  changestats:   S    ***********/
WtD_CHANGESTAT_FN(d_sum);
WtD_CHANGESTAT_FN(d_sum_pow);

/********************  changestats:   T    ***********/
WtD_CHANGESTAT_FN(d_transitiveweights_threshold); WtS_CHANGESTAT_FN(s_transitiveweights_threshold);
WtD_CHANGESTAT_FN(d_transitiveweighs); WtS_CHANGESTAT_FN(s_transitiveweights);

              
#endif
