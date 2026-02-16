/*  File src/changestats_diversity.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_changestat.h"
#include "ergm_storage.h"

#define mk_c_nodecovrange(_LABEL_, _OUT_, _IN_, _BSHIFT_)               \
  C_CHANGESTAT_FN(c_ ## _LABEL_ ## covrange) {                          \
    if( _OUT_ ) {                                                       \
      double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf; \
                                                                        \
      EXEC_THROUGH_OUTEDGES(tail, e, u, {                               \
          oldmin = MIN(oldmin, INPUT_PARAM[u-1 _BSHIFT_]);              \
          oldmax = MAX(oldmax, INPUT_PARAM[u-1 _BSHIFT_]);              \
                                                                        \
          if(!edgestate || u!=head){                                    \
            newmin = MIN(newmin, INPUT_PARAM[u-1 _BSHIFT_]);            \
            newmax = MAX(newmax, INPUT_PARAM[u-1 _BSHIFT_]);            \
          }                                                             \
        });                                                             \
                                                                        \
      if(!edgestate){                                                   \
        newmin = MIN(newmin, INPUT_PARAM[head-1 _BSHIFT_]);             \
        newmax = MAX(newmax, INPUT_PARAM[head-1 _BSHIFT_]);             \
      }                                                                 \
                                                                        \
      CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin); \
    }                                                                   \
                                                                        \
    if( _IN_ ) {                                                        \
      double oldmin = R_PosInf, oldmax = R_NegInf, newmin = R_PosInf, newmax = R_NegInf; \
                                                                        \
      EXEC_THROUGH_INEDGES(head, e, u, {                                \
          oldmin = MIN(oldmin, INPUT_PARAM[u-1]);                       \
          oldmax = MAX(oldmax, INPUT_PARAM[u-1]);                       \
                                                                        \
          if(!edgestate || u!=tail){                                    \
            newmin = MIN(newmin, INPUT_PARAM[u-1]);                     \
            newmax = MAX(newmax, INPUT_PARAM[u-1]);                     \
          }                                                             \
        });                                                             \
                                                                        \
      if(!edgestate){                                                   \
        newmin = MIN(newmin, INPUT_PARAM[tail-1]);                      \
        newmax = MAX(newmax, INPUT_PARAM[tail-1]);                      \
      }                                                                 \
                                                                        \
      CHANGE_STAT[0] += (isinf(newmax) ? 0 : newmax-newmin) - (isinf(oldmax) ? 0 : oldmax-oldmin); \
    }                                                                   \
  }

mk_c_nodecovrange(node, TRUE, TRUE,)
mk_c_nodecovrange(nodeo, TRUE, FALSE,)
mk_c_nodecovrange(nodei, FALSE, TRUE,)
mk_c_nodecovrange(b1, TRUE, FALSE, - BIPARTITE)


#define mk_nodefactordistinct(_LABEL_, _EFF_NODES_, _OUT_, _IN_, _BSHIFT_) \
                                                                        \
  I_CHANGESTAT_FN(i_ ## _LABEL_ ## factordistinct) {                    \
    unsigned int ncats = IINPUT_PARAM[0];                               \
                                                                        \
    ALLOC_STORAGE(ncats*(_EFF_NODES_), unsigned int, freqs);            \
                                                                        \
    EXEC_THROUGH_NET_EDGES(tail, head, e, {                             \
        if(_OUT_ && IINPUT_PARAM[head _BSHIFT_]) freqs[(tail-1)*ncats + IINPUT_PARAM[head _BSHIFT_] - 1]++; \
        if(_IN_ && IINPUT_PARAM[tail]) freqs[(head-1 _BSHIFT_)*ncats + IINPUT_PARAM[tail] - 1]++; \
      });                                                               \
  }                                                                     \
                                                                        \
  C_CHANGESTAT_FN(c_ ## _LABEL_ ## factordistinct) {                    \
    unsigned int ncats = IINPUT_PARAM[0];                               \
    GET_STORAGE(unsigned int, freqs);                                   \
    int change = edgestate ? -1 : +1;                                   \
    if(_OUT_ && IINPUT_PARAM[head _BSHIFT_]) {                                   \
      unsigned int oldfreq = freqs[(tail-1)*ncats + IINPUT_PARAM[head _BSHIFT_] - 1]; \
      unsigned int newfreq = oldfreq + change;                          \
      CHANGE_STAT[0] += (newfreq!=0) - (oldfreq!=0);                    \
    }                                                                   \
                                                                        \
    if(_IN_ && IINPUT_PARAM[tail]) {                                    \
      unsigned int oldfreq = freqs[(head-1 _BSHIFT_)*ncats + IINPUT_PARAM[tail] - 1]; \
      unsigned int newfreq = oldfreq + change;                          \
      CHANGE_STAT[0] += (newfreq!=0) - (oldfreq!=0);                    \
    }                                                                   \
  }                                                                     \
                                                                        \
                                                                        \
  U_CHANGESTAT_FN(u_ ## _LABEL_ ## factordistinct) {                    \
    unsigned int ncats = IINPUT_PARAM[0];                               \
    GET_STORAGE(unsigned int, freqs);                                   \
    int change = edgestate ? -1 : +1;                                   \
    if(_OUT_ && IINPUT_PARAM[head _BSHIFT_])                            \
      freqs[(tail-1)*ncats + IINPUT_PARAM[head _BSHIFT_] - 1] += change;    \
    if(_IN_ && IINPUT_PARAM[tail])                                      \
      freqs[(head-1 _BSHIFT_)*ncats + IINPUT_PARAM[tail] - 1] += change;    \
  }

mk_nodefactordistinct(node, N_NODES, TRUE, TRUE,)
mk_nodefactordistinct(nodeo, N_NODES, TRUE, FALSE,)
mk_nodefactordistinct(nodei, N_NODES, FALSE, TRUE,)
mk_nodefactordistinct(b1, BIPARTITE, TRUE, FALSE, - BIPARTITE)
mk_nodefactordistinct(b2, N_NODES-BIPARTITE, FALSE, TRUE, - BIPARTITE)
