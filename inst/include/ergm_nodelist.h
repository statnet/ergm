/*  File inst/include/ergm_nodelist.h in package ergm, part of the Statnet
 *  suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#ifndef _ERGM_NODELIST_H_
#define _ERGM_NODELIST_H_

#include "ergm_edgetree_types.h"

/*
   This data structure stores a collection of nodes. The nodepos field allows
   for lookups: if node is present in the nodelist, nodepos[node] gives its
   (1-based) index in the nodes field. nodepos should be an array of length
   N_NODES + 1, which is passed in to the constructor by the client code to
   allow sharing of nodepos across different NodeLists, as when partitioning
   the set of all nodes based on e.g. a nodal attribute value.
*/

typedef struct {
  Vertex *nodes;
  int length;
  int *nodepos;
} NodeList;

static inline NodeList *NodeListInitialize(int nnodes, int *nodepos) {
  NodeList *nodelist = R_Calloc(1, NodeList);
  nodelist->nodes = R_Calloc(nnodes + 1, Vertex);
  nodelist->nodepos = nodepos;
  return nodelist;
}

static inline void NodeListDestroy(NodeList *nodelist) {
  R_Free(nodelist->nodes);
  R_Free(nodelist);
}

static inline void NodeListInsert(NodeList *nodelist, Vertex node) {
  nodelist->length++;
  nodelist->nodes[nodelist->length] = node;
  nodelist->nodepos[node] = nodelist->length;
}

static inline void NodeListDelete(NodeList *nodelist, Vertex node) {
  nodelist->nodes[nodelist->nodepos[node]] = nodelist->nodes[nodelist->length];
  nodelist->nodepos[nodelist->nodes[nodelist->length]] = nodelist->nodepos[node];
  nodelist->nodepos[node] = 0;
  nodelist->length--;
}

static inline void NodeListToggle(NodeList *nodelist, Vertex node) {
  if(nodelist->nodepos[node]) {
    NodeListDelete(nodelist, node);
  } else {
    NodeListInsert(nodelist, node);
  }
}

#endif
