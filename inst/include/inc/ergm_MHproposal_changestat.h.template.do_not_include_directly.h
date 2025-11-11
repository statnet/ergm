#include "ergm_Rutil.h"

static inline bool ETYPE(ChangeStats1_changed)(ETYPE(MHProposal) *MHp, ETYPE(Network) *nwp, ETYPE(Model) *m, EWTTYPE edgestate) {
  ETYPE(ChangeStats1)(*MHp->toggletail, *MHp->togglehead, IFEWT(*MHp->toggleweight,) nwp, m, edgestate);
  for (unsigned int i = 0; i < m->n_stats; i++) if (m->workspace[i] != 0) return(TRUE);
  return(FALSE);
}

static inline bool ETYPE(ChangeStats_changed)(ETYPE(MHProposal) *MHp, ETYPE(Network) *nwp, ETYPE(Model) *m) {
  ETYPE(ChangeStats)(MHp->ntoggles, MHp->toggletail, MHp->togglehead, IFEWT(MHp->toggleweight,) nwp, m);
  for (unsigned int i = 0; i < m->n_stats; i++) if (m->workspace[i] != 0) return(TRUE);
  return(FALSE);
}

#define GET_CHANGESTATS_MODEL(into)                             \
  {                                                             \
    SEXP slot = getListElement(MHp->R, "ChangeStat_pos");       \
    if (slot == R_NilValue) into = NULL;                        \
    else into = MH_AUX_STORAGE_NUM(Rf_asInteger(slot));         \
  }
