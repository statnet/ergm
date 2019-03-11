/*  File src/changestats_dgw_sp.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution
 *
 *  Copyright 2003-2019 Statnet Commons
 */
#ifndef _CHANGESTATS_DGW_SP_H_
#define _CHANGESTATS_DGW_SP_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5

/* /\*DSP calculation functions*\/ */
/* static inline void dspUTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void dspOTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void dspITP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void dspRTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void dspOSP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void dspISP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */


/* /\*ESP calculation functions*\/ */
/* static inline void espUTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void espOTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void espITP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void espRTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void espOSP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */
/* static inline void espISP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs); */

/*Changescore functions*/
C_CHANGESTAT_FN(c_desp);
I_CHANGESTAT_FN(i_dgwesp);
C_CHANGESTAT_FN(c_dgwesp);

/*Changescore functions*/
C_CHANGESTAT_FN(c_ddsp);
I_CHANGESTAT_FN(i_dgwdsp);
C_CHANGESTAT_FN(c_dgwdsp);

/*Changescore functions*/
I_CHANGESTAT_FN(i_dnsp);
C_CHANGESTAT_FN(c_dnsp);
I_CHANGESTAT_FN(i_dgwnsp);
C_CHANGESTAT_FN(c_dgwnsp);

#endif // _CHANGESTATS_DGW_SP_H_
