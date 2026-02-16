/*  File inst/include/cpp/ergm_wtproposal.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#pragma once

#include "ergm_proposal_base.h"
#include "ergm_wtMHproposal.h"

namespace ergm {
inline namespace v1 {

template<typename StorageType = void>
class ErgmCppWtProposal : public ErgmCppProposalBase<WtMHProposal, StorageType> {
public:
  ErgmCppWtProposal(WtMHProposal* mhp)
    : ErgmCppProposalBase<WtMHProposal, StorageType>(mhp),
      weight((mhp->toggleweight) {}

  double* weight;
};

} // namespace v1
} // namespace ergm
