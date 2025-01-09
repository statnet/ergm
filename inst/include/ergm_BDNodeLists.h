/*  File inst/include/ergm_BDNodeLists.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2025 Statnet Commons
 */
#ifndef _ERGM_BDNODELISTS_H_
#define _ERGM_BDNODELISTS_H_

#include "ergm_nodelist.h"

/*
   The BDNodeLists data structure stores NodeLists with a two-fold stratification
   based on the ego node's combined_vattr value and the alter node's bd_vattr
   value. For a bd_vattr value i and a combined_vattr value j, tails[i][j] (resp.
   heads[i][j], boths[i][j]) is the NodeList containing the set of egos with
   combined_vattr value j that are of strictly submaximal degree as tails (resp.
   heads, both tails and heads) to alters of bd_vattr value i. For undirected
   networks, tails[i][j], heads[i][j], and boths[i][j] all coincide. The notion
   of strict submaximality of degree is defined in terms of the supplied maxout,
   maxin degree bounds and outdegree, indegree degree states. (maxout, maxin,
   outdegree, and indegree are arrays with two indices, the first being the
   bd_vattr value of the alter and the second being the 1-based nodal index
   of the ego.)
*/

typedef struct {
  NodeList ***tails;
  NodeList ***heads;
  NodeList ***boths;

  int **tailpos;
  int **headpos;
  int **bothpos;

  int *combined_vattr;
  int *bd_vattr;
  
  int combined_nlevels;
  int bd_nlevels;

  int directed;
} BDNodeLists;

static inline BDNodeLists *BDNodeListsInitialize(int **maxout,
                                                 int **maxin,
                                                 int **outdegree,
                                                 int **indegree,
                                                 int *combined_vattr,
                                                 int combined_nlevels,
                                                 int *bd_vattr,
                                                 int bd_nlevels,
                                                 int *combined_vattr_counts,
                                                 Network *nwp) {
  BDNodeLists *lists = R_Calloc(1, BDNodeLists);

  // do some copying
  lists->directed = DIRECTED;

  lists->combined_vattr = combined_vattr;
  lists->bd_vattr = bd_vattr;

  lists->combined_nlevels = combined_nlevels;
  lists->bd_nlevels = bd_nlevels;

  // set up node lists
  lists->bothpos = R_Calloc(bd_nlevels, int *);
  lists->tailpos = DIRECTED ? R_Calloc(bd_nlevels, int *) : lists->bothpos;
  lists->headpos = DIRECTED ? R_Calloc(bd_nlevels, int *) : lists->bothpos;

  lists->boths = R_Calloc(bd_nlevels, NodeList **);
  lists->tails = DIRECTED ? R_Calloc(bd_nlevels, NodeList **) : lists->boths;
  lists->heads = DIRECTED ? R_Calloc(bd_nlevels, NodeList **) : lists->boths;

  for(int i = 0; i < bd_nlevels; i++) {
    lists->bothpos[i] = R_Calloc(N_NODES + 1, int);
    if(DIRECTED) {
      lists->tailpos[i] = R_Calloc(N_NODES + 1, int);
      lists->headpos[i] = R_Calloc(N_NODES + 1, int);
    }

    lists->boths[i] = R_Calloc(combined_nlevels, NodeList *);
    if(DIRECTED) {
      lists->tails[i] = R_Calloc(combined_nlevels, NodeList *);
      lists->heads[i] = R_Calloc(combined_nlevels, NodeList *);
    }

    for(int j = 0; j < combined_nlevels; j++) {
      lists->boths[i][j] = NodeListInitialize(combined_vattr_counts[j], lists->bothpos[i]);
      if(DIRECTED) {
        lists->tails[i][j] = NodeListInitialize(combined_vattr_counts[j], lists->tailpos[i]);
        lists->heads[i][j] = NodeListInitialize(combined_vattr_counts[j], lists->headpos[i]);
      }
    }
  }

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int attr_val = combined_vattr[vertex];

    for(int i = 0; i < bd_nlevels; i++) {
      if(DIRECTED) {
        if(indegree[i][vertex] < maxin[i][vertex] && outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(lists->boths[i][attr_val], vertex);
        } else if(outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(lists->tails[i][attr_val], vertex);
        } else if(indegree[i][vertex] < maxin[i][vertex]) {
          NodeListInsert(lists->heads[i][attr_val], vertex);
        }
      } else if(indegree[i][vertex] + outdegree[i][vertex] < maxout[i][vertex]) {
        NodeListInsert(lists->boths[i][attr_val], vertex);
      }
    }
  }

  return lists;
}

static inline void BDNodeListsDestroy(BDNodeLists *lists) {
  for(int i = 0; i < lists->bd_nlevels; i++) {
    for(int j = 0; j < lists->combined_nlevels; j++) {
      NodeListDestroy(lists->boths[i][j]);
      if(lists->directed) {
        NodeListDestroy(lists->tails[i][j]);
        NodeListDestroy(lists->heads[i][j]);
      }
    }

    R_Free(lists->boths[i]);
    if(lists->directed) {
      R_Free(lists->tails[i]);
      R_Free(lists->heads[i]);
    }

    R_Free(lists->bothpos[i]);
    if(lists->directed) {
      R_Free(lists->tailpos[i]);
      R_Free(lists->headpos[i]);
    }
  }

  R_Free(lists->boths);
  if(lists->directed) {
    R_Free(lists->tails);
    R_Free(lists->heads);
  }

  R_Free(lists->bothpos);
  if(lists->directed) {
    R_Free(lists->tailpos);
    R_Free(lists->headpos);
  }

  R_Free(lists);
}

// update NodeLists as appropriate for a toggle with tail, head maximality
// changes given by tailcondition and headcondition flags
static inline void BDNodeListsToggleIf(Vertex tail, Vertex head, BDNodeLists *lists, int tailcondition, int headcondition) {
  if(tailcondition) {
    int tailattr = lists->combined_vattr[tail];
    int headattr = lists->bd_vattr[head];
    if(lists->directed && (lists->bothpos[headattr][tail] || lists->headpos[headattr][tail])) {
      NodeListToggle(lists->boths[headattr][tailattr], tail);
      NodeListToggle(lists->heads[headattr][tailattr], tail);
    } else {
      NodeListToggle(lists->tails[headattr][tailattr], tail);
    }
  }

  if(headcondition) {
    int tailattr = lists->bd_vattr[tail];
    int headattr = lists->combined_vattr[head];
    if(lists->directed && (lists->bothpos[tailattr][head] || lists->tailpos[tailattr][head])) {
      NodeListToggle(lists->boths[tailattr][headattr], head);
      NodeListToggle(lists->tails[tailattr][headattr], head);
    } else {
      NodeListToggle(lists->heads[tailattr][headattr], head);
    }
  }
}

// count of egos with the same combined_vattr value as tail that are of strictly
// submaximal degree as tails to alters with the same bd_vattr value as head
static inline int BDNodeListsTailCount(Vertex tail, Vertex head, BDNodeLists *lists) {
  int tailattr = lists->combined_vattr[tail];
  int headattr = lists->bd_vattr[head];
  return lists->tails[headattr][tailattr]->length + lists->directed*lists->boths[headattr][tailattr]->length;
}

// count of egos with the same combined_vattr value as head that are of strictly
// submaximal degree as heads to alters with the same bd_vattr value as tail
static inline int BDNodeListsHeadCount(Vertex tail, Vertex head, BDNodeLists *lists) {
  int tailattr = lists->bd_vattr[tail];
  int headattr = lists->combined_vattr[head];
  return lists->heads[tailattr][headattr]->length + lists->directed*lists->boths[tailattr][headattr]->length;
}

#endif
