#ifndef _ERGM_CONSTANTS_H_
#define _ERGM_CONSTANTS_H_

// Macro indicating the version of the C API.
#define ERGM_API 4.0

typedef enum MCMCStatus_enum {
  MCMC_OK = 0,
  MCMC_TOO_MANY_EDGES = 1,
  MCMC_MH_FAILED = 2
} MCMCStatus;

#define ERGM_STATE_EMPTY_NET 1u
#define ERGM_STATE_NO_INIT_S 2u
#define ERGM_STATE_NO_INIT_PROP 4u

typedef enum ErgmStateExtFlag_enum {
  ERGM_STATE_R_CHANGED = -1,
  ERGM_STATE_C_CHANGED = +1,
  ERGM_STATE_RECONCILED = 0
} ErgmStateExtFlag;

#endif // _ERGM_CONSTANTS_H_
