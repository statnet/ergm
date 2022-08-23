/*  File inst/include/ergm_nodelist.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_NODELIST_H_
#define _ERGM_NODELIST_H_

#include "ergm_MHproposal.h"
#include "ergm_edgetree_types.h"
#include "ergm_changestat.h"

typedef struct {
  Vertex *nodes;
  int length;
  int *nodepos;
} NodeList;

static inline NodeList *NodeListInitialize(int nnodes, int *nodepos) {
  NodeList *nodelist = Calloc(1, NodeList);
  nodelist->nodes = Calloc(nnodes + 1, Vertex);
  nodelist->nodepos = nodepos;
  return nodelist;
}

static inline void NodeListDestroy(NodeList *nodelist) {
  Free(nodelist->nodes);
  Free(nodelist);
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

static inline void NodeListToggleKnown(NodeList *nodelist, Vertex node, int nodeflag) {
  if(nodeflag) {
    NodeListDelete(nodelist, node);
  } else {
    NodeListInsert(nodelist, node);
  }
}

static inline void NodeListToggle(NodeList *nodelist, Vertex node) {
  if(nodelist->nodepos[node]) {
    NodeListDelete(nodelist, node);      
  } else {
    NodeListInsert(nodelist, node);      
  }
}

#endif
