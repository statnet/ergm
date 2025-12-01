#pragma once

#include "FixedArray.h"
#include "ergm_auxstorage_proxy.h"
#include "ergm_R_proxy.h"

template<typename ModelTermType, typename StorageType = void>
class ErgmCppModelTermBase {
public:

  ErgmCppModelTermBase(ModelTermType* mtp)
    : stat(mtp->dstats, mtp->nstats),
      dinput(mtp->inputparams, mtp->ninputparams),
      iinput(mtp->iinputparams, mtp->niinputparams),
      dattrib(mtp->attrib, (mtp->attrib && mtp->inputparams) ? (mtp->ninputparams - (mtp->attrib - mtp->inputparams)) : 0),
      iattrib(mtp->iattrib, (mtp->iattrib && mtp->iinputparams) ? (mtp->niinputparams - (mtp->iattrib - mtp->iinputparams)) : 0),
      storage(reinterpret_cast<StorageType*&>(mtp->storage)),
      aux_storage(mtp),
      R(mtp),
      mtp_(mtp)
  {}

  FixedArray<double> stat;
  const FixedArray<double> dinput;
  const FixedArray<int> iinput;
  const FixedArray<double> dattrib;
  const FixedArray<int> iattrib;
  StorageType*& storage;

  // aux_storage access via proxy
  AuxStorageProxy<ModelTermType> aux_storage;

  RListProxy<ModelTermType> R;

private:
  ModelTermType* mtp_;
};
