#pragma once

#include "ergm_proposal_base.h"
#include "ergm_MHproposal.h"

template<typename StorageType = void>
using ErgmCppProposal = ErgmCppProposalBase<MHProposal, StorageType>;
