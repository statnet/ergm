#pragma once

#include "ergm_proposal_base.h"
#include "ergm_wtMHproposal.h"

template<typename StorageType = void>
class ErgmCppWtProposal : public ErgmCppProposalBase<WtMHProposal, StorageType> {
public:
  ErgmCppWtProposal(WtMHProposal* mhp)
    : ErgmCppProposalBase<WtMHProposal, StorageType>(mhp),
      weight((mhp->toggleweight) {}

  double* weight;
};
