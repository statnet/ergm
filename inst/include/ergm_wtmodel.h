/*  File inst/include/ergm_wtmodel.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_WTMODEL_H_
#define _ERGM_WTMODEL_H_

#include "ergm_wtedgetree.h"
#include "ergm_wtchangestat.h"
#include "ergm_wtMHproposal.h"

#include "ergm_edgetype_set_double.h"
#include "inc/ergm_model.h.template.do_not_include_directly.h"
#include "ergm_edgetype_unset.h"

#define WtFOR_EACH_TERM(m) for(WtModelTerm *mtp = (m)->termarray; mtp < (m)->termarray + (m)->n_terms; mtp++)

#define WtEXEC_THROUGH_TERMS(m, subroutine)				\
  WtFOR_EACH_TERM(m){							\
    subroutine;								\
  }

#define WtFOR_EACH_TERM_INREVERSE(m) for(WtModelTerm *mtp = (m)->termarray + (m)->n_terms - 1; mtp >= (m)->termarray; mtp--)

#define WtEXEC_THROUGH_TERMS_INREVERSE(m, subroutine)			\
  WtFOR_EACH_TERM_INREVERSE(m){						\
    subroutine;								\
  }

#define WtEXEC_THROUGH_TERMS_INTO(m, output, subroutine)		\
  WtFOR_EACH_TERM(m){							\
    double *dstats = output + mtp->statspos;				\
    subroutine;								\
  }

#define WtSEND_X_SIGNAL(nwp, m, MHp, type, data)                        \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  WtEXEC_THROUGH_TERMS((m), {                                           \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });

#define WtSEND_X_SIGNAL_INREVERSE(nwp, m, MHp, type, data)              \
  WtEXEC_THROUGH_TERMS_INREVERSE((m), {                                 \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });                                                                 \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp));

#define WtSEND_X_SIGNAL_INTO(nwp, m, MHp, output, type, data)           \
  if((MHp) && ((WtMHProposal*)(MHp))->x_func) ((WtMHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  WtEXEC_THROUGH_TERMS_INTO((m), output, {                              \
      if(mtp->x_func){                                                  \
        mtp->dstats = dstats;                                           \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
      }                                                                 \
    });

#endif

