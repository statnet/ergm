#ifndef _ERGM_CHANGESTAT_OPERATOR_H_
#define _ERGM_CHANGESTAT_OPERATOR_H_

#include "ergm_model.h"
#include "ergm_storage.h"

Model *unpack_Model_as_double(double **x);

I_CHANGESTAT_FN(i_OnAuxnet);
C_CHANGESTAT_FN(c_OnAuxnet);
U_CHANGESTAT_FN(u_OnAuxnet);
F_CHANGESTAT_FN(f_OnAuxnet);

#endif
