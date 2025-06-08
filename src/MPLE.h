/*  File src/MPLE.h in package ergm, part of the Statnet suite of packages for
 *  network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef MPLE_H
#define MPLE_H

#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_model.h"
#include "ergm_rlebdm.h"
#include "ergm_state.h"
#include "ergm_khash.h"
#include "ergm_kvec.h"

typedef struct{
  int edges, nonedges;
} ENE;

static inline unsigned int kh_hash_DVec(double *newRow, unsigned int rowLength){
  /* Cast all pointers to unsigned char pointers, since data need to
     be fed to the hash function one byte at a time. */
  unsigned char *cnewRow = (unsigned char *) newRow;
  unsigned int crowLength = rowLength * sizeof(double);

  unsigned int hash=0;

#define HASH_LOOP(hash, keybyte){ hash+=keybyte; hash+=(hash<<10); hash^= (hash>>6); }
  for(unsigned int i=0; i<crowLength; i++) HASH_LOOP(hash, cnewRow[i]);
#undef HASH_LOOP

  hash += (hash<<3);
  hash ^= (hash>>11);
  hash += (hash<<15);

  return(hash);
}

#define kh_DVec_hash_func(key) (khint32_t)(kh_hash_DVec(key, h->l))
#define kh_DVec_hash_equal(a,b) !memcmp(a, b, h->l*sizeof(double))

KHASH_INIT(DVecMapENE, double*, ENE, true, kh_DVec_hash_func, kh_DVec_hash_equal, size_t l;)
typedef khash_t(DVecMapENE) StoreDVecMapENE;

StoreDVecMapENE *MpleInit_hash_wl_RLE(ErgmState *s, RLEBDM1D *wl, Edge maxNumDyads);

#endif
