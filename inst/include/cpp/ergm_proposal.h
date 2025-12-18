#pragma once

#include "ergm_proposal_base.h"
#include "ergm_MHproposal.h"

namespace ergm {
inline namespace v1 {

template<typename StorageType = void>
using ErgmCppProposal = ErgmCppProposalBase<MHProposal, StorageType>;

} // namespace v1
} // namespace ergm
