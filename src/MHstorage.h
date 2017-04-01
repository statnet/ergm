/*  File src/changestat.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */

/* Storage utilities for MH proposal clients */
#define MH_STORAGE (/* (stored_type *) */ MHp->storage)
#define MH_ALLOC_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (MH_STORAGE = calloc(nmemb, sizeof(stored_type)));
#define MH_GET_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) MH_STORAGE;

#define MH_AUX_STORAGE (/* (stored_type *) */ MHp->aux_storage[(unsigned int) MH_INPUTS[0]])
#define MH_GET_AUX_STORAGE(stored_type, store_into) stored_type *store_into = MH_AUX_STORAGE;
#define MH_AUX_STORAGE_NUM(ind) (/* (stored_type *) */ MHp->aux_storage[(unsigned int) MH_INPUTS[ind]])
#define MH_GET_AUX_STORAGE_NUM(stored_type, store_into, ind) stored_type *store_into = MH_AUX_STORAGE_NUM(ind);

#ifndef STRICT_MH_HEADERS

#define STORAGE MH_STORAGE
#define ALLOC_STORAGE MH_ALLOC_STORAGE
#define GET_STORAGE MH_GET_STORAGE
#define AUX_STORAGE MH_AUX_STORAGE
#define GET_AUX_STORAGE MH_GET_AUX_STORAGE
#define AUX_STORAGE_NUM MH_AUX_STORAGE_NUM
#define GET_AUX_STORAGE_NUM MH_GET_AUX_STORAGE_NUM

#endif // STRICT_MH_HEADERS
