#ifndef _CHANGESTATS_DGW_SP_ML_H_
#define _CHANGESTATS_DGW_SP_ML_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_storage.h"
#include "ergm_changestat_multilayer.h"
#include "changestats_dgw_sp.h"

static inline int ergm_Check2Path(unsigned int e11, unsigned int e12, unsigned int e21, unsigned int e22, unsigned int any_order){
  if(any_order) return((e11&&e22) || (e12&&e21));
  else return(e11&&e22);
}
// FIXME: Optimize
static inline int ergm_c_LayerLogic2Path(Vertex tail1, Vertex head1, Vertex tail2, Vertex head2,
					 StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order,
					 int c11, int c12, int c21, int c22){
  if(!ll1->onwp->directed_flag) any_order = TRUE;
  
  unsigned int e11o = c11? c11==-1 : ML_GETWT(ll1, tail1, head1), e22o = c22? c22==-1 : ML_GETWT(ll2, tail2, head2), e12o=0, e21o=0;
  unsigned int e11n = c11? c11==+1 : e11o, e22n = c22? c22==+1 : e22o, e12n=0, e21n=0;
  
  if(any_order){
    e12o = c12? c12==-1 : ML_GETWT(ll2, tail1, head1);
    e21o = c21? c21==-1 : ML_GETWT(ll1, tail2, head2);
    e12n = c12? c12==+1 : e12o;
    e21n = c21? c21==+1 : e21o;
  }

  return ergm_Check2Path(e11n, e12n, e21n, e22n, any_order)
    - ergm_Check2Path(e11o, e12o, e21o, e22o, any_order);
}


static inline unsigned int ergm_LayerLogic2Path(Vertex tail1, Vertex head1, Vertex tail2, Vertex head2, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order){
  if(!ll1->onwp->directed_flag) any_order = TRUE;

  unsigned int e11 = ML_GETWT(ll1, tail1, head1), e22 = ML_GETWT(ll2, tail2, head2), e12, e21;
  if(any_order){
    e12 = ML_GETWT(ll2, tail1, head1);
    e21 = ML_GETWT(ll1, tail2, head2);
  }
  return ergm_Check2Path(e11, e12, e21, e22, any_order);
}


/* /\*DSP calculation functions*\/ */
/* static inline void dspUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* /\* static inline void dspRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); *\/ */
/* static inline void dspOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void dspISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, unsigned int any_order, int nd, double *dvec, double *cs); */


/* /\*ESP calculation functions*\/ */
/* static inline void espUTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espOTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espITP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* /\* static inline void espRTP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); *\/ */
/* static inline void espOSP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */
/* static inline void espISP_ML_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, StoreLayerLogic *ll0, StoreLayerLogic *ll1, StoreLayerLogic *ll2, StoreLayerLogic *ll3, unsigned int any_order, int nd, double *dvec, double *cs); */

/* /\*Changescore functions*\/ */
/* C_CHANGESTAT_FN(c_desp_ML); */
/* I_CHANGESTAT_FN(i_dgwesp_ML); */
/* C_CHANGESTAT_FN(c_dgwesp_ML); */

/* /\*Changescore functions*\/ */
/* C_CHANGESTAT_FN(c_ddsp_ML); */
/* I_CHANGESTAT_FN(i_dgwdsp_ML); */
/* C_CHANGESTAT_FN(c_dgwdsp_ML); */

/* /\*Changescore functions*\/ */
/* I_CHANGESTAT_FN(i_dnsp_ML); */
/* C_CHANGESTAT_FN(c_dnsp_ML); */
/* I_CHANGESTAT_FN(i_dgwnsp_ML); */
/* C_CHANGESTAT_FN(c_dgwnsp_ML); */

#endif // _CHANGESTATS_DGW_SP_ML_H_
