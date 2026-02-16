/*  File inst/include/cpp/ergm_proposal_base.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */

#pragma once
#include <cstddef>

#include "FixedArray.h"
#include "ergm_auxstorage_proxy.h"
#include "ergm_R_proxy.h"

// Template C++ wrapper for MHProposal and WtMHProposal structs
// Usage: ErgmCppProposalBase<MHProposal> or ErgmCppProposalBase<WtMHProposal>

namespace ergm {
inline namespace v1 {

template<typename ProposalType, typename StorageType = void>
class ErgmCppProposalBase {
public:
  ErgmCppProposalBase(ProposalType* mhp)
    : ptr(mhp),
      dinput(mhp->inputs, mhp->ninputs),
      iinput(mhp->iinputs, mhp->niinputs),
      storage(reinterpret_cast<StorageType*&>(mhp->storage)),
      aux_storage(mhp),
      R(mhp->R) {}

  ProposalType* ptr;

  // Direct member access
  Edge& size = ptr->ntoggles;
  Vertex* tail = ptr->toggletail;
  Vertex* head = ptr->togglehead;
  double& logratio = ptr->logratio;
  const FixedArray<double> dinput;
  const FixedArray<int> iinput;
  StorageType*& storage;

  // Aux storage proxy
  AuxStorageProxy<ProposalType> aux_storage;

  // R list and attribute proxy
  RListProxy R;
};

} // namespace v1
} // namespace ergm
