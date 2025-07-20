#pragma once

#include "FixedArray.h"

extern "C" {
#include "../ergm_changestat.h"
}

#define C_CHANGESTAT_CPP(name, impl)                                    \
  extern "C" void c_ ## name (Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm mt(mtp);                                           \
    {impl;}                                                             \
  }

#define S_CHANGESTAT_CPP(name, impl)                                    \
  extern "C" void s_ ## name (ModelTerm *mtp, Network *nwp) { \
    ErgmCppNetwork nw(nwp);                                             \
    ErgmCppModelTerm mt(mtp);                                           \
    {impl;}                                                             \
  }

class ErgmCppModelTerm {
public:
  ErgmCppModelTerm(ModelTerm* mtp)
    : mtp_(mtp),
      stat(mtp->dstats, mtp->nstats),
      dinput(mtp->inputparams, mtp->ninputparams),
      iinput(mtp->iinputparams, mtp->niinputparams),
      dattrib(mtp->attrib, (mtp->attrib && mtp->inputparams) ? (mtp->ninputparams - (mtp->attrib - mtp->inputparams)) : 0),
      iattrib(mtp->iattrib, (mtp->iattrib && mtp->iinputparams) ? (mtp->niinputparams - (mtp->iattrib - mtp->iinputparams)) : 0)
  {}

  FixedArray<double> stat;
  const FixedArray<double> dinput;
  const FixedArray<int> iinput;
  const FixedArray<double> dattrib;
  const FixedArray<int> iattrib;

private:
  ModelTerm* mtp_;
};
