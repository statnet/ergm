/*  File inst/include/cpp/ergm_modelterm_base.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#pragma once

#include "FixedArray.h"
#include "ergm_auxstorage_proxy.h"
#include "ergm_R_proxy.h"

namespace ergm {
inline namespace v1 {

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
      R(mtp->R),
      ext_state(mtp->ext_state),
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

  RListProxy R;
  RListProxy ext_state;

private:
  ModelTermType* mtp_;
};

} // namespace v1
} // namespace ergm
