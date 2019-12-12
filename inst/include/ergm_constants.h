#ifndef _ERGM_CONSTANTS_H_
#define _ERGM_CONSTANTS_H_

typedef enum MCMCStatus_enum {
  MCMC_OK = 0,
  MCMC_TOO_MANY_EDGES = 1,
  MCMC_MH_FAILED = 2
} MCMCStatus;

typedef enum ErgmStateExtFlag_enum {
  ERGM_STATE_R_CHANGED = -1,
  ERGM_STATE_C_CHANGED = +1,
  ERGM_STATE_RECONCILED = 0
} ErgmStateExtFlag;

#endif // _ERGM_CONSTANTS_H_
