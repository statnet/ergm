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
