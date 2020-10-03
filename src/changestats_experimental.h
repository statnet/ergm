/*  File src/changestats_experimental.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2020 Statnet Commons
 */
#ifndef CHANGESTATS_EXPERIMENTAL_H
#define CHANGESTATS_EXPERIMENTAL_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
/********************  changestats:  A    ***********/
D_CHANGESTAT_FN(d_b1kappa);
D_CHANGESTAT_FN(d_b1share);
D_CHANGESTAT_FN(d_altistar);
D_CHANGESTAT_FN(d_altostar);
/********************  changestats:  B    ***********/
D_CHANGESTAT_FN(d_berninhom);
D_CHANGESTAT_FN(d_biduration);
D_CHANGESTAT_FN(d_bimix);
D_CHANGESTAT_FN(d_bkappa);
/********************  changestats:  C    ***********/
/********************  changestats:  D    ***********/
D_CHANGESTAT_FN(d_degreep);
D_CHANGESTAT_FN(d_degreep_by_attr);
D_CHANGESTAT_FN(d_degreep_w_homophily);
D_CHANGESTAT_FN(d_dissolve);
D_CHANGESTAT_FN(d_duration);
/********************  changestats:  E    ***********/
D_CHANGESTAT_FN(d_b2kappa);
/********************  changestats:  F    ***********/
D_CHANGESTAT_FN(d_factor);
/********************  changestats:  G    ***********/
D_CHANGESTAT_FN(d_geodegree);
D_CHANGESTAT_FN(d_geospartner);
D_CHANGESTAT_FN(d_gwb1);
D_CHANGESTAT_FN(d_gwd);
D_CHANGESTAT_FN(d_gwdegree706);
D_CHANGESTAT_FN(d_gwdegreealpha);
D_CHANGESTAT_FN(d_gwdegreelambda);
D_CHANGESTAT_FN(d_gwb2);
D_CHANGESTAT_FN(d_gwb1share);
D_CHANGESTAT_FN(d_gwb2share);
/********************  changestats:   H    ***********/
D_CHANGESTAT_FN(d_heideriandynamic);
D_CHANGESTAT_FN(d_hiertriad);
  double numposthree (Vertex t, Network *nwp);
D_CHANGESTAT_FN(d_hiertriaddegree);
/********************  changestats:   I    ***********/
D_CHANGESTAT_FN(d_icvar);
D_CHANGESTAT_FN(d_idc);
D_CHANGESTAT_FN(d_intransitivedynamic);
D_CHANGESTAT_FN(d_intransitivity);
/********************  changestats:   K    ***********/
D_CHANGESTAT_FN(d_kappa);
/********************  changestats:   L    ***********/
/********************  changestats:   M    ***********/
D_CHANGESTAT_FN(d_monopolymixmat);
/********************  changestats:   N    ***********/
/********************  changestats:   O    ***********/
/********************  changestats:   R    ***********/
/********************  changestats:   S    ***********/
D_CHANGESTAT_FN(d_simmeliandynamic);
D_CHANGESTAT_FN(d_spatial);
/********************  changestats:   T    ***********/
D_CHANGESTAT_FN(d_transitivedynamic);
D_CHANGESTAT_FN(d_transitivity);
       
#endif
