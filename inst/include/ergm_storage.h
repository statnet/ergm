/*  File inst/include/ergm_storage.h in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */

#ifndef _ERGM_STORAGE_H_
#define _ERGM_STORAGE_H_

#include "ergm_variadic_macros.h"

/*** Storage utilities ***/

/** Private storage **/

// Pointer to a term's private storage
#define STORAGE (mtp->storage)

// 1. Declares a stored_type *store_into.
// 2. Allocates a vector of `nmemb` elements of type `stored_type`.
// 3. Saves its pointer to private storage.
// 4. Also assigns its pointer to store_into.
#define ALLOC_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (STORAGE = R_Calloc(nmemb, stored_type));

// 1. Declares a stored_type *store_into.
// 2. Assigns pointer to private storage to store_into.
#define GET_STORAGE(stored_type, store_into) stored_type *store_into = (stored_type *) STORAGE;


/** Public (auxiliary) storage **/

// Pointer to an auxiliary term's assigned storage slot or to the storage slot of a statistic's first requested auxiliary.
#define AUX_STORAGE (mtp->aux_storage[mtp->aux_slots[0]])
#define N_AUX (mtp->n_aux)

// Should be used by auxiliary terms' initialization or updating function only.
// 1. Declares a stored_type *store_into.
// 2. Allocates a vector of `nmemb` elements of type `stored_type`.
// 3. Saves its pointer to its assigned public storage slot.
// 4. Also assigns its pointer to store_into.
#define ALLOC_AUX_STORAGE(nmemb, stored_type, store_into) stored_type *store_into = (stored_type *) (AUX_STORAGE = R_Calloc(nmemb, stored_type));
// 1. Declares a stored_type *store_into.
// 2. Assigns pointer to its assigned auxiliary storage slot (or, for a statistic, its first requested auxiliary) to store_into.
#define _GET_AUX_STORAGE2(stored_type, store_into) stored_type *store_into = (stored_type *) AUX_STORAGE;
// Pointer to the storage slot of a statistic's ind'th requested auxiliary.
#define AUX_STORAGE_NUM(ind) (mtp->aux_storage[mtp->aux_slots[ind]])
// 1. Declares a stored_type *store_into.
// 2. Assigns pointer to its ind'th requested auxiliary to store_into.
#define _GET_AUX_STORAGE3(ind, stored_type, store_into) stored_type *store_into = (stored_type *) AUX_STORAGE_NUM(ind);
// For backwards-compatibility. TODO: Delete around 4.10 release.
#define GET_AUX_STORAGE_NUM(stored_type, store_into, ind) _GET_AUX_STORAGE3(ind, stored_type, store_into)
// This version takes 2 or 3 arguments; if 3 arguments, the first argument is the slot number.
#define GET_AUX_STORAGE(...) _GET_OVERRIDE3(__VA_ARGS__, _GET_AUX_STORAGE3, _GET_AUX_STORAGE2,,)(__VA_ARGS__)

/* Allocate a sociomatrix as auxiliary storage. */
#define ALLOC_AUX_SOCIOMATRIX(stored_type, store_into)			\
  /* Note: the following code first sets up a 2D array indexed from 0, then shifts all pointers by -1 so that sm[t][h] would work for vertex IDs. */ \
  ALLOC_AUX_STORAGE(N_TAILS, stored_type*, store_into);			\
  Dyad sm_size = BIPARTITE? N_TAILS*N_HEADS : DIRECTED ? N_NODES*N_NODES : N_NODES*(N_NODES+1)/2; /* For consistency, and possible future capabilities, include diagonal: */ \
  ALLOC_STORAGE(sm_size, stored_type, data); /* A stored_type* to data. */ \
  Dyad pos = 0;	  /* Start of the next row's data in the data vector. */ \
  for(Vertex t=0; t<N_TAILS; t++){                                      \
  /* First set up the pointer to the right location in the data vector, */ \
  if(BIPARTITE){							\
  store_into[t] = data+pos - N_TAILS; /* This is so that store_into[t][h=BIPARTITE] would address the 0th element of that row. */ \
  pos += N_HEADS;							\
  }else if(DIRECTED){							\
    store_into[t] = data+pos;						\
    pos += N_HEADS;							\
  }else{ /* Undirected. */						\
    store_into[t] = data+pos - t; /* tail <= head, so this is so that store_into[t][h=t] would address the 0th element of that row. */ \
    pos += N_HEADS-t+1; /* Each row has N_NODES - t + 1 elements (including diagonal). */ \
  }									\
  store_into[t]--; /* Now, shift the pointer by -1. */			\
  }									\
									\
  store_into--; /* Shift the pointer array by -1. */			\
  AUX_STORAGE = store_into; /* This is needed to make sure the pointer array itself is updated. */

/* R_Free a sociomatrix in auxiliary storage. */
/* If we hadn't shifted the pointers by -1, this would not have been
   necessary. We need to shift the array back into place so that it's
   automatically deallocated. */
#define FREE_AUX_SOCIOMATRIX AUX_STORAGE = (void **)AUX_STORAGE + 1;

#endif // _ERGM_STORAGE_H_
