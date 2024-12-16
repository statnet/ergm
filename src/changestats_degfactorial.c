#include "ergm_changestat.h"
#include "ergm_storage.h"

#define mk_c_sum_lfactorial_degreem1(_LABEL_, _OUT_, _IN_)              \
                                                                        \
  C_CHANGESTAT_FN(c_sum_lfactorial_ ## _LABEL_ ## degreem1){            \
    int echange=edgestate ? -1:+1;                                      \
                                                                        \
    Vertex otd = 0, ntd = 0, ohd = 0, nhd = 0;                          \
                                                                        \
    if(_OUT_){otd = OUT_DEG[tail] + ( _IN_ ? IN_DEG[tail] : 0); ntd = otd + echange;} \
    if(_IN_){ohd = ( _OUT_ ? OUT_DEG[head] : 0) + IN_DEG[head]; nhd = ohd + echange;} \
                                                                        \
    CHANGE_STAT[0] += (ntd?lgammafn(ntd):0) - (otd?lgammafn(otd):0)     \
      + (nhd?lgammafn(nhd):0) - (ohd?lgammafn(ohd):0);                  \
  }

mk_c_sum_lfactorial_degreem1(, TRUE, TRUE)
mk_c_sum_lfactorial_degreem1(o, TRUE, FALSE)
mk_c_sum_lfactorial_degreem1(i, FALSE, TRUE)
