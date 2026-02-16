/*  File inst/include/ergm_model.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_MODEL_H_
#define _ERGM_MODEL_H_

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_MHproposal.h"

#include "ergm_edgetype_set_binary.h"
#include "inc/ergm_model.h.template.do_not_include_directly.h"
#include "ergm_edgetype_unset.h"

#define FOR_EACH_TERM(m) for(ModelTerm *mtp = (m)->termarray; mtp < (m)->termarray + (m)->n_terms; mtp++)


#define EXEC_THROUGH_TERMS(m, subroutine)				\
  FOR_EACH_TERM(m){							\
    subroutine;								\
  }

#define FOR_EACH_TERM_INREVERSE(m) for(ModelTerm *mtp = (m)->termarray + (m)->n_terms - 1; mtp >= (m)->termarray; mtp--)

#define EXEC_THROUGH_TERMS_INREVERSE(m, subroutine)			\
  FOR_EACH_TERM_INREVERSE(m){						\
    subroutine;								\
  }

#define EXEC_THROUGH_TERMS_INTO(m, output, subroutine)			\
  FOR_EACH_TERM(m){							\
    double *dstats = output + mtp->statspos;				\
    subroutine;								\
  }

#define SEND_X_SIGNAL(nwp, m, MHp, type, data)                          \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  EXEC_THROUGH_TERMS((m), {                                             \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });

#define SEND_X_SIGNAL_INREVERSE(nwp, m, MHp, type, data)                 \
  EXEC_THROUGH_TERMS_INREVERSE((m), {                                   \
      if(mtp->x_func)                                                   \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
    });                                                                 \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp));

#define SEND_X_SIGNAL_INTO(nwp, m, MHp, output, type, data)              \
  if((MHp) && ((MHProposal*)(MHp))->x_func) ((MHProposal*)(MHp))->x_func((type), (data), (MHp), (nwp)); \
  EXEC_THROUGH_TERMS_INTO((m), output, {                                \
      if(mtp->x_func){                                                  \
        mtp->dstats = dstats;                                           \
        (*(mtp->x_func))((type), (data), (mtp), (nwp));                 \
      }                                                                 \
    });

#endif

