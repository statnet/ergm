/*  File inst/include/ergm_constants.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_CONSTANTS_H_
#define _ERGM_CONSTANTS_H_

/* Macros indicating the version of the C API. This is used primarily
   if the client package wants to conditionally compile based on
   'ergm' version. They should be updated with every minor version
   update. */
#define ERGM_API_MAJOR 4
#define ERGM_API_MINOR 9

/* ABI version: this should be updated only when the ABI changes, even
   if the change is 100% source-compatible. This includes, for
   example, adding an element to one of the structs.

   Then, it should be set to a signed int with value major*1e6 +
   minor. Make sure to remove any leading zeros!
*/
#define ERGM_ABI_VERSION 4000009

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
