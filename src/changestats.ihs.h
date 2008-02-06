#ifndef CHANGESTATS_IHS_H
#define CHANGESTATS_IHS_H

#include "edgetree.ihs.h"
#include "changestats.h"

/****************************************************/
/* changestat function prototypes, 
   plus a few supporting function prototypes */
/********************  changestats:  A    ***********/
CHANGESTAT_FN(d_b1kappa)
CHANGESTAT_FN(d_b1share)
CHANGESTAT_FN(d_altistar)
CHANGESTAT_FN(d_altostar)
/********************  changestats:  B    ***********/
CHANGESTAT_FN(d_berninhom)
CHANGESTAT_FN(d_biduration)
CHANGESTAT_FN(d_bimix)
CHANGESTAT_FN(d_bkappa)
/********************  changestats:  C    ***********/
/********************  changestats:  D    ***********/
CHANGESTAT_FN(d_degreep)
CHANGESTAT_FN(d_degreep_by_attr)
CHANGESTAT_FN(d_degreep_w_homophily)
CHANGESTAT_FN(d_dissolve)
CHANGESTAT_FN(d_duration)
/********************  changestats:  E    ***********/
CHANGESTAT_FN(d_b2kappa)
/********************  changestats:  F    ***********/
CHANGESTAT_FN(d_factor)
CHANGESTAT_FN(d_formation)
/********************  changestats:  G    ***********/
CHANGESTAT_FN(d_geodegree)
CHANGESTAT_FN(d_geospartner)
CHANGESTAT_FN(d_gwb1)
CHANGESTAT_FN(d_gwd)
CHANGESTAT_FN(d_gwdegree706)
CHANGESTAT_FN(d_gwdegreealpha)
CHANGESTAT_FN(d_gwdegreelambda)
CHANGESTAT_FN(d_gwb2)
CHANGESTAT_FN(d_gwb1share)
CHANGESTAT_FN(d_gwb2share)
/********************  changestats:   H    ***********/
CHANGESTAT_FN(d_heideriandynamic)
CHANGESTAT_FN(d_hiertriad)
  double numposthree (Vertex t, Network *nwp);
CHANGESTAT_FN(d_hiertriaddegree)
/********************  changestats:   I    ***********/
CHANGESTAT_FN(d_icvar)
CHANGESTAT_FN(d_idc)
CHANGESTAT_FN(d_intransitivedynamic)
CHANGESTAT_FN(d_intransitivity)
/********************  changestats:   K    ***********/
CHANGESTAT_FN(d_kappa)
/********************  changestats:   L    ***********/
/********************  changestats:   M    ***********/
CHANGESTAT_FN(d_monopolymixmat)
/********************  changestats:   N    ***********/
/********************  changestats:   O    ***********/
/********************  changestats:   R    ***********/
/********************  changestats:   S    ***********/
CHANGESTAT_FN(d_simmeliandynamic)
CHANGESTAT_FN(d_spatial)
/********************  changestats:   T    ***********/
CHANGESTAT_FN(d_transitivedynamic)
CHANGESTAT_FN(d_transitivity)
       
#endif
