#ifndef _ERGM_CONSTANTS_H_
#define _ERGM_CONSTANTS_H_

typedef enum MCMCStatus_enum {
  MCMC_OK = 0,
  MCMC_TOO_MANY_EDGES = 1,
  MCMC_MH_FAILED = 2
} MCMCStatus;

#endif // _ERGM_CONSTANTS_H_
