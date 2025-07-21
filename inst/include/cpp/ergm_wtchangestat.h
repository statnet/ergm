#pragma once

#include "FixedArray.h"

extern "C" {
#include "../ergm_wtchangestat.h"
}

#define WtC_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void c_ ## name (Vertex tail, Vertex head, double weight, WtModelTerm *mtp, WtNetwork *nwp, double edgestate) { \
    ErgmCppWtNetwork nw(nwp);                                             \
    ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                                           \
    {impl;}                                                             \
  }

#define WtS_CHANGESTAT_CPP(name, impl, ...)                                    \
  extern "C" void s_ ## name (WtModelTerm *mtp, WtNetwork *nwp) { \
    ErgmCppWtNetwork nw(nwp);                                             \
    ErgmCppWtModelTerm<__VA_ARGS__> mt(mtp);                                           \
    {impl;}                                                             \
  }


template<typename StorageType = void>
class ErgmCppWtModelTerm {
public:
  ErgmCppWtModelTerm(WtModelTerm* mtp)
    : mtp_(mtp),
      stat(mtp->dstats, mtp->nstats),
      dinput(mtp->inputparams, mtp->ninputparams),
      iinput(mtp->iinputparams, mtp->niinputparams),
      dattrib(mtp->attrib, (mtp->attrib && mtp->inputparams) ? (mtp->ninputparams - (mtp->attrib - mtp->inputparams)) : 0),
      iattrib(mtp->iattrib, (mtp->iattrib && mtp->iinputparams) ? (mtp->niinputparams - (mtp->iattrib - mtp->iinputparams)) : 0),
      storage(static_cast<StorageType*>(mtp->storage))
  {}

  FixedArray<double> stat;
  const FixedArray<double> dinput;
  const FixedArray<int> iinput;
  const FixedArray<double> dattrib;
  const FixedArray<int> iattrib;
  StorageType* storage;

  void set_storage(StorageType* ptr) {
    storage = ptr;
    mtp_->storage = static_cast<void*>(ptr);
  }

  // Encapsulate aux_storage access
  void* aux_storage(size_t i) const {
    return mtp_->aux_storage[mtp_->aux_slots[i]];
  }
  void* operator[](size_t i) const {
    return aux_storage(i);
  }

private:
  WtModelTerm* mtp_;
};
