/*  File src/changestats_dgw_sp.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5

/*DSP calculation functions*/
void dspUTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void dspOTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void dspITP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void dspRTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void dspOSP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void dspISP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);


/*ESP calculation functions*/
void espUTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espOTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espITP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espRTP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espOSP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);
void espISP_calc(Edge ntoggles, Vertex *tails, Vertex *heads, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs);

/*Changescore functions*/
D_CHANGESTAT_FN(d_desp);
D_CHANGESTAT_FN(d_dgwesp);

/*Changescore functions*/
D_CHANGESTAT_FN(d_ddsp);
D_CHANGESTAT_FN(d_dgwdsp);

/*Changescore functions*/
D_CHANGESTAT_FN(d_dnsp);
D_CHANGESTAT_FN(d_dgwnsp);

#endif
