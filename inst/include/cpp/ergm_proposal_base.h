
#pragma once
#include <cstddef>

#include "FixedArray.h"
#include "ergm_auxstorage_proxy.h"
#include "ergm_R_proxy.h"

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
