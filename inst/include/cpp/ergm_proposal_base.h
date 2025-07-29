
#pragma once
#include <cstddef>

#include "ergm_auxstorage_proxy.h"

// Template C++ wrapper for MHProposal and WtMHProposal structs
// Usage: ErgmCppProposalBase<MHProposal> or ErgmCppProposalBase<WtMHProposal>

template<typename ProposalType, typename StorageType = void>
class ErgmCppProposalBase {
public:
  ErgmCppProposalBase(ProposalType* mhp)
    : ptr(mhp),
      dinput(mhp->inputs, mhp->ninputs),
      iinput(mhp->iinputs, mhp->niinputs),
      storage(reinterpret_cast<StorageType*&>(mhp->storage)),
      aux_storage(mhp) {}

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
};
