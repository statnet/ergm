#pragma once

#include "FixedArray.h"

template<typename ModelTermType, typename StorageType = void>
class ErgmCppModelTermBase {
public:
  ErgmCppModelTermBase(ModelTermType* mtp)
    : mtp_(mtp),
      stat(mtp->dstats, mtp->nstats),
      dinput(mtp->inputparams, mtp->ninputparams),
      iinput(mtp->iinputparams, mtp->niinputparams),
      dattrib(mtp->attrib, (mtp->attrib && mtp->inputparams) ? (mtp->ninputparams - (mtp->attrib - mtp->inputparams)) : 0),
      iattrib(mtp->iattrib, (mtp->iattrib && mtp->iinputparams) ? (mtp->niinputparams - (mtp->iattrib - mtp->iinputparams)) : 0),
      storage(static_cast<StorageType*>(mtp->storage)),
      aux_storage(mtp)
  {}

  FixedArray<double> stat;
  const FixedArray<double> dinput;
  const FixedArray<int> iinput;
  const FixedArray<double> dattrib;
  const FixedArray<int> iattrib;
  StorageType* storage;

  // aux_storage access via proxy
  struct AuxStorageProxy {
    ModelTermType* mtp_;
    AuxStorageProxy(ModelTermType* mtp) : mtp_(mtp) {}
    void* operator[](size_t i) const {
      return mtp_->aux_storage[mtp_->aux_slots[i]];
    }
  } aux_storage;

  void set_storage(StorageType* ptr) {
    storage = ptr;
    mtp_->storage = static_cast<void*>(ptr);
  }

private:
  ModelTermType* mtp_;
};
