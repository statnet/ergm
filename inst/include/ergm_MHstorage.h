/*  File inst/include/ergm_MHstorage.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

#ifndef _ERGM_MHSTORAGE_H_
#define _ERGM_MHSTORAGE_H_

/* Storage utilities for MH proposal clients */
#define MH_STORAGE (/* (stored_type *) */ MHp->storage)
#define MH_N_AUX (MHp->n_aux)
#define MH_ALLOC_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (MH_STORAGE = R_Calloc(nmemb, stored_type));
#define MH_GET_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) MH_STORAGE;

#define MH_AUX_STORAGE (/* (stored_type *) */ MHp->aux_storage)
#define _MH_GET_AUX_STORAGE2(stored_type, store_into) stored_type *store_into = MH_AUX_STORAGE[MHp->aux_slots[0]];
#define MH_AUX_STORAGE_NUM(ind) /* (stored_type *) */ MH_AUX_STORAGE[MHp->aux_slots[ind]]
#define _MH_GET_AUX_STORAGE3(ind, stored_type, store_into) stored_type *store_into = MH_AUX_STORAGE_NUM(ind);
// For backwards-compatibility. TODO: Delete around 3.10 release.
#define MH_GET_AUX_STORAGE_NUM(stored_type, store_into, ind) _MH_GET_AUX_STORAGE3(ind, stored_type, store_into)
#define MH_GET_AUX_STORAGE(...) _GET_OVERRIDE3(__VA_ARGS__, _MH_GET_AUX_STORAGE3, _MH_GET_AUX_STORAGE2,,)(__VA_ARGS__)


#ifdef STUBFILE
#define STRICT_MH_HEADERS
#endif

#ifndef STRICT_MH_HEADERS

#define STORAGE MH_STORAGE
#define N_AUX MH_N_AUX
#define ALLOC_STORAGE MH_ALLOC_STORAGE
#define GET_STORAGE MH_GET_STORAGE
#define AUX_STORAGE MH_AUX_STORAGE
#define GET_AUX_STORAGE MH_GET_AUX_STORAGE
#define AUX_STORAGE_NUM MH_AUX_STORAGE_NUM
#define GET_AUX_STORAGE_NUM MH_GET_AUX_STORAGE_NUM

#endif // STRICT_MH_HEADERS

#endif // _ERGM_MHSTORAGE_H_
