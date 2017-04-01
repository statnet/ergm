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
#define STORAGE (/* (stored_type *) */ MHp->storage)
#define ALLOC_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (STORAGE = calloc(nmemb, sizeof(stored_type)));
#define GET_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) STORAGE;

#define AUX_STORAGE (/* (stored_type *) */ MHp->aux_storage[(unsigned int) MH_INPUTS[0]])
#define GET_AUX_STORAGE(stored_type, store_into) stored_type *store_into = AUX_STORAGE;
#define AUX_STORAGE_NUM(ind) (/* (stored_type *) */ MHp->aux_storage[(unsigned int) MH_INPUTS[ind]])
#define GET_AUX_STORAGE_NUM(stored_type, store_into, ind) stored_type *store_into = AUX_STORAGE_NUM(ind);
