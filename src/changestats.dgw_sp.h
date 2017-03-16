#ifndef CHANGESTATS_H
#define CHANGESTATS_H

#include "edgetree.h"
#include "changestat.h"

#define ESPUTP 0
#define ESPOTP 1
#define ESPITP 2
#define ESPRTP 3
#define ESPOSP 4
#define ESPISP 5

/*DSP calculation functions*/
static inline void dspUTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void dspOTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void dspITP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void dspRTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void dspOSP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void dspISP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);


/*ESP calculation functions*/
static inline void espUTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void espOTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void espITP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void espRTP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void espOSP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);
static inline void espISP_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, int nd, double *dvec, double *cs, unsigned int cur_weight);

/*Changescore functions*/
C_CHANGESTAT_FN(c_desp);
C_CHANGESTAT_FN(c_dgwesp);

/*Changescore functions*/
C_CHANGESTAT_FN(c_ddsp);
C_CHANGESTAT_FN(c_dgwdsp);

/*Changescore functions*/
C_CHANGESTAT_FN(c_dnsp);
C_CHANGESTAT_FN(c_dgwnsp);

#endif
