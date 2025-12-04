/*  File src/changestats_dyad_ind.c.template.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */

// This file has been converted to idiomatic C++ using the
// ergm_changestat.h and ergm_wtchangestat.h macro definitions.
// Each function is wrapped using C_CHANGESTAT_CPP or WtC_CHANGESTAT_CPP.

#include <cmath>
#include "ergm_storage.h"

using namespace std;

// Ugly but functional: otherwise, *C_CHANGESTAT_CPP will clobber SVARIANT().
#define c_SVARIANT(name) SVARIANT(c_ ## name)
#define i_SVARIANT(name) SVARIANT(i_ ## name)
#define u_SVARIANT(name) SVARIANT(u_ ## name)
#define f_SVARIANT(name) SVARIANT(f_ ## name)

/********************  changestats:  A    ***********/

/*
 changestat: c_absdiff
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(absdiff), {
  double p = mt.dinput[0];
  if(p == 1.0) {
    mt.stat[0] = ECHANGE(fabs(mt.dinput[tail] - mt.dinput[head]));
  } else {
    mt.stat[0] = ECHANGE(pow(fabs(mt.dinput[tail] - mt.dinput[head]), p));
  }
})

/*
 changestat: c_absdiffcat
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(absdiffcat), {
  double tailval = mt.dattrib[tail - 1];
  double headval = mt.dattrib[head - 1];
  double absdiff = fabs(tailval - headval);
  if (absdiff > 0) {
    for (unsigned int j = 0; j < mt.stat.size(); j++) {
      mt.stat[j] += (absdiff == mt.dinput[j]) ? ECHANGE1 : 0.0;
    }
  }
})

/*
 changestat: attrcov - with storage
*/
struct SVARIANT(AttrcovStorage) {
  int* nodecov;
  double** mat;
};

ETYPE(I_CHANGESTAT_CPP)(SVARIANT(attrcov), {
  ALLOC_STORAGE(1, SVARIANT(AttrcovStorage), sto);
  sto->nodecov = INTEGER(mt.R["nodecov"]);

  int nr = asInteger(mt.R["nr"]);
  int nc = asInteger(mt.R["nc"]);

  sto->mat = R_Calloc(nc, double*);
  sto->mat[0] = REAL(mt.R["mat"]);
  for(int i = 1; i < nc; i++) {
    sto->mat[i] = sto->mat[i - 1] + nr;
  }
}, SVARIANT(AttrcovStorage))

ETYPE(C_CHANGESTAT_CPP)(SVARIANT(attrcov), {
  mt.stat[0] += ECHANGE(mt.storage->mat[mt.storage->nodecov[head]][mt.storage->nodecov[tail]]);
}, SVARIANT(AttrcovStorage))

ETYPE(F_CHANGESTAT_CPP)(SVARIANT(attrcov), {
  R_Free(mt.storage->mat);
}, SVARIANT(AttrcovStorage))

/*
 changestat: c_b2cov
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(b2cov), {
  unsigned int oshift = mt.dinput.size() / mt.stat.size();
  Vertex nb1 = nw.bip;

  for(unsigned int j = 0, o = 0; j < mt.stat.size(); j++, o += oshift) {
    double sum = mt.dattrib[head - nb1 + o - 1];
    mt.stat[j] += ECHANGE(sum);
  }
})

/*
 changestat: c_b2factor
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(b2factor), {
  int headpos = mt.dattrib[head - 1 - nw.bip];
  if (headpos != -1) {
    mt.stat[headpos] += ECHANGE1;
  }
})

/********************  changestats:  D    ***********/

/*
 changestat: c_density
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(density), {
  Dyad ndyads = nw.n * (nw.n - 1) / (nw.dir ? 1 : 2);
  mt.stat[0] += ECHANGE(1.0 / ndyads);
})

/*
 changestat: c_diff
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(diff), {
  double p = mt.dinput[0];
  int mul = mt.iinput[0];
  int sign_code = mt.iinput[1];

  double change = (mt.dinput[tail] - mt.dinput[head]) * mul;

  switch(sign_code) {
    case 1:  // identity
      break;
    case 2:  // abs
      change = fabs(change);
      break;
    case 3:  // positive only
      change = change < 0 ? 0 : change;
      break;
    case 4:  // negative only
      change = change > 0 ? 0 : change;
      break;
    default:
      error("Invalid sign action code passed to c_diff.");
      break;
  }

  if(p == 0.0) {
    change = (change > 0) ? 1.0 : ((change < 0) ? -1.0 : 0.0);
  } else if(p != 1.0) {
    change = pow(change, p);
  }

  mt.stat[0] += ECHANGE(change);
})

/********************  changestats:  E    ***********/

/*
 changestat: c_edgecov
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(edgecov), {
  int noffset = nw.bip;
  int nrow = (noffset > 0) ? noffset : (int)mt.dinput[0];

  double val = mt.dattrib[(head - 1 - noffset) * nrow + (tail - 1)];
  mt.stat[0] += ECHANGE(val);
})

/*
 changestat: c_edges
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(edges), {
  mt.stat[0] = ECHANGE1;
})

/********************  changestats:  L    ***********/

/*
 changestat: c_meandeg
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(meandeg), {
  if(nw.dir) {
    mt.stat[0] = ECHANGE(1.0 / nw.n);
  } else {
    mt.stat[0] = ECHANGE(2.0 / nw.n);
  }
})

/*
 changestat: c_mixmat - General mixing matrix (mm) implementation
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(mixmat), {
  unsigned int symm = mt.iinput[0] & 1;
  unsigned int marg = mt.iinput[0] & 2;
  const int* tx = mt.iinput.data();
  const int* hx = nw.bip ? mt.iinput.data() : mt.iinput.data() + nw.n;
  const int* cells = nw.bip ? mt.iinput.data() + nw.n + 1 : mt.iinput.data() + nw.n * 2 + 1;

  unsigned int diag = tx[tail] == tx[head] && hx[tail] == hx[head];
  for(unsigned int j = 0; j < mt.stat.size(); j++) {
    unsigned int thmatch = tx[tail] == cells[j*2] && hx[head] == cells[j*2+1];
    unsigned int htmatch = tx[head] == cells[j*2] && hx[tail] == cells[j*2+1];

    int w = (nw.dir || nw.bip) ? thmatch :
      (symm ? (thmatch || htmatch) : (thmatch + htmatch)) * (symm && marg && diag ? 2 : 1);
    if(w) {
      mt.stat[j] += ECHANGE(w);
    }
  }
})

/*
 changestat: c_nodecov
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodecov), {
  unsigned int oshift = mt.dinput.size() / mt.stat.size();

  for(unsigned int j = 0, o = 0; j < mt.stat.size(); j++, o += oshift) {
    double sum = mt.dattrib[tail + o - 1] + mt.dattrib[head + o - 1];
    mt.stat[j] += ECHANGE(sum);
  }
})

/*
 changestat: c_nodefactor
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodefactor), {
  int tailpos = mt.iattrib[tail - 1];
  int headpos = mt.iattrib[head - 1];
  if (tailpos != -1) {
    mt.stat[tailpos] += ECHANGE1;
  }
  if (headpos != -1) {
    mt.stat[headpos] += ECHANGE1;
  }
})

/*
 changestat: c_nodeicov
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodeicov), {
  unsigned int oshift = mt.dinput.size() / mt.stat.size();

  for(unsigned int j = 0, o = 0; j < mt.stat.size(); j++, o += oshift) {
    double sum = mt.dattrib[head + o - 1];
    mt.stat[j] += ECHANGE(sum);
  }
})

/*
 changestat: c_nodeifactor
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodeifactor), {
  int headpos = mt.dattrib[head - 1];
  if (headpos != -1) {
    mt.stat[headpos] += ECHANGE1;
  }
})

/*
 changestat: c_nodematch
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodematch), {
  Vertex ninputs = mt.dinput.size() - nw.n;

  double matchval = mt.dinput[tail + ninputs - 1];
  if (matchval == mt.dinput[head + ninputs - 1]) {
    if (ninputs == 0) {
      mt.stat[0] += ECHANGE1;
    } else {
      for (unsigned int j = 0; j < ninputs; j++) {
        if (matchval == mt.dinput[j]) {
          mt.stat[j] += ECHANGE1;
        }
      }
    }
  }
})

/*
 changestat: nodemix - with storage
*/
struct SVARIANT(NodemixStorage) {
  int* nodecov;
  int** indmat;
};

ETYPE(I_CHANGESTAT_CPP)(SVARIANT(nodemix), {
  ALLOC_STORAGE(1, SVARIANT(NodemixStorage), sto);
  sto->nodecov = INTEGER(mt.R["nodecov"]);

  int nr = asInteger(mt.R["nr"]);
  int nc = asInteger(mt.R["nc"]);

  sto->indmat = R_Calloc(nr, int*);
  sto->indmat[0] = INTEGER(mt.R["indmat"]);
  for(int i = 1; i < nr; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + nc;
  }
}, SVARIANT(NodemixStorage))

ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodemix), {
  int index = mt.storage->indmat[mt.storage->nodecov[tail]][mt.storage->nodecov[head]];
  if(index >= 0) {
    mt.stat[index] += ECHANGE1;
  }
}, SVARIANT(NodemixStorage))

ETYPE(F_CHANGESTAT_CPP)(SVARIANT(nodemix), {
  R_Free(mt.storage->indmat);
}, SVARIANT(NodemixStorage))

/*
 changestat: c_nodeocov
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodeocov), {
  unsigned int oshift = mt.dinput.size() / mt.stat.size();

  for(unsigned int j = 0, o = 0; j < mt.stat.size(); j++, o += oshift) {
    double sum = mt.dattrib[tail + o - 1];
    mt.stat[j] += ECHANGE(sum);
  }
})

/*
 changestat: c_nodeofactor
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(nodeofactor), {
  int tailpos = mt.dattrib[tail - 1];
  if (tailpos != -1) {
    mt.stat[tailpos] += ECHANGE1;
  }
})

/********************  changestats:  O    ***********/

/*
 changestat: c_receiver
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(receiver), {
  Vertex deg;

  for(unsigned int j = 0; j < mt.stat.size(); j++) {
    deg = (Vertex)mt.dinput[j];
    if(deg == head) {
      mt.stat[j] += ECHANGE1;
      break;
    }
  }
})

/********************  changestats:  S    ***********/

/*
 changestat: c_sender
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(sender), {
  Vertex deg;

  for(unsigned int j = 0; j < mt.stat.size(); j++) {
    deg = (Vertex)mt.dinput[j];
    if(deg == tail) {
      mt.stat[j] += ECHANGE1;
      break;
    }
  }
})

/*
 changestat: c_smalldiff
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(smalldiff), {
  mt.stat[0] += (fabs(mt.dattrib[tail - 1] - mt.dattrib[head - 1]) > mt.dinput[0])
    ? 0.0 : ECHANGE1;
})

/*
 changestat: c_sociality
*/
ETYPE(C_CHANGESTAT_CPP)(SVARIANT(sociality), {
  unsigned int ninputs = mt.dinput.size();
  unsigned int nstats = mt.stat.size();

  if(ninputs > nstats + 1) {
    // match on attributes
    double tailattr = mt.dattrib[tail - 1 + nstats + 1];
    double headattr = mt.dattrib[head - 1 + nstats + 1];

    if(tailattr == headattr) {
      // Find tail
      for(unsigned int j = 0; j < nstats; j++) {
        if((Vertex)mt.dinput[j] == tail) {
          mt.stat[j] += ECHANGE1;
          break;
        }
      }
      // Find head
      for(unsigned int j = 0; j < nstats; j++) {
        if((Vertex)mt.dinput[j] == head) {
          mt.stat[j] += ECHANGE1;
          break;
        }
      }
    }
  } else {
    // Find tail
    for(unsigned int j = 0; j < nstats; j++) {
      if((Vertex)mt.dinput[j] == tail) {
        mt.stat[j] += ECHANGE1;
        break;
      }
    }
    // Find head
    for(unsigned int j = 0; j < nstats; j++) {
      if((Vertex)mt.dinput[j] == head) {
        mt.stat[j] += ECHANGE1;
        break;
      }
    }
  }
})

