/*  File inst/include/ergm_BDNodeLists.h in package ergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2022 Statnet Commons
 */
#ifndef _ERGM_BDNODELISTS_H_
#define _ERGM_BDNODELISTS_H_

#include "ergm_nodelist.h"

typedef struct {  
  NodeList ***tails_by_attr;
  NodeList ***heads_by_attr;
  NodeList ***boths_by_attr;
  
  NodeList ***tails;
  NodeList ***heads;
  NodeList ***boths;

  int **tailpos;
  int **headpos;
  int **bothpos;
  
  int *bd_vattr;
    
  int directed;

  int total_nlevels;
  int bd_nlevels;
} BDNodeLists;

static inline BDNodeLists *BDNodeListsInitialize(int **maxout, 
                                                 int **maxin, 
                                                 int **indegree,
                                                 int **outdegree,
                                                 int *strat_vattr, 
                                                 int strat_nlevels, 
                                                 int *blocks_vattr, 
                                                 int blocks_nlevels, 
                                                 int *bd_vattr, 
                                                 int bd_nlevels, 
                                                 int *jointattrcounts, 
                                                 Network *nwp) {
  BDNodeLists *lists = Calloc(1, BDNodeLists);

  // do some copying    
  lists->directed = DIRECTED;
  
  lists->bd_nlevels = bd_nlevels;
  lists->total_nlevels = strat_nlevels*blocks_nlevels*bd_nlevels;

  lists->bd_vattr = bd_vattr;

  // set up node lists
  lists->bothpos = Calloc(bd_nlevels, int *);
  lists->tailpos = DIRECTED ? Calloc(bd_nlevels, int *) : lists->bothpos;  
  lists->headpos = DIRECTED ? Calloc(bd_nlevels, int *) : lists->bothpos;  
  
  for(int i = 0; i < bd_nlevels; i++) {
    lists->bothpos[i] = Calloc(N_NODES + 1, int);
    if(DIRECTED) {
      lists->tailpos[i] = Calloc(N_NODES + 1, int);  
      lists->headpos[i] = Calloc(N_NODES + 1, int); 
    }
  }
    
  lists->boths_by_attr = Calloc(lists->total_nlevels, NodeList **);
  lists->tails_by_attr = DIRECTED ? Calloc(lists->total_nlevels, NodeList **) : lists->boths_by_attr;
  lists->heads_by_attr = DIRECTED ? Calloc(lists->total_nlevels, NodeList **) : lists->boths_by_attr;
  
  for(int i = 0; i < lists->total_nlevels; i++) {
    lists->boths_by_attr[i] = Calloc(bd_nlevels, NodeList *);
    if(DIRECTED) {
      lists->tails_by_attr[i] = Calloc(bd_nlevels, NodeList *);
      lists->heads_by_attr[i] = Calloc(bd_nlevels, NodeList *);
    }
    
    for(int j = 0; j < bd_nlevels; j++) {
      lists->boths_by_attr[i][j] = NodeListInitialize(jointattrcounts[i], lists->bothpos[j]);
      if(DIRECTED) {
        lists->tails_by_attr[i][j] = NodeListInitialize(jointattrcounts[i], lists->tailpos[j]);
        lists->heads_by_attr[i][j] = NodeListInitialize(jointattrcounts[i], lists->headpos[j]);
      }
    }
  }

  lists->boths = Calloc(N_NODES + 1, NodeList **);
  lists->tails = DIRECTED ? Calloc(N_NODES + 1, NodeList **) : lists->boths;
  lists->heads = DIRECTED ? Calloc(N_NODES + 1, NodeList **) : lists->boths;

  for(Vertex vertex = 1; vertex <= N_NODES; vertex++) {
    int attr_val = strat_vattr[vertex]*blocks_nlevels*bd_nlevels + blocks_vattr[vertex]*bd_nlevels + bd_vattr[vertex];
    
    for(int i = 0; i < bd_nlevels; i++) {
      if(DIRECTED) {
        if(indegree[i][vertex] < maxin[i][vertex] && outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(lists->boths_by_attr[attr_val][i], vertex);        
        } else if(outdegree[i][vertex] < maxout[i][vertex]) {
          NodeListInsert(lists->tails_by_attr[attr_val][i], vertex);
        } else if(indegree[i][vertex] < maxin[i][vertex]) {
          NodeListInsert(lists->heads_by_attr[attr_val][i], vertex);
        }
      } else if(indegree[i][vertex] + outdegree[i][vertex] < maxout[i][vertex]) {
        NodeListInsert(lists->boths_by_attr[attr_val][i], vertex);        
      }
    }
    
    lists->boths[vertex] = lists->boths_by_attr[attr_val];
    if(DIRECTED) {
      lists->tails[vertex] = lists->tails_by_attr[attr_val];
      lists->heads[vertex] = lists->heads_by_attr[attr_val];
    }
  }

  return lists;
}

static inline void BDNodeListsDestroy(BDNodeLists *lists) {
  for(int i = 0; i < lists->total_nlevels; i++) {
    for(int j = 0; j < lists->bd_nlevels; j++) {
      NodeListDestroy(lists->boths_by_attr[i][j]);
      if(lists->directed) {
        NodeListDestroy(lists->tails_by_attr[i][j]);
        NodeListDestroy(lists->heads_by_attr[i][j]);
      }
    }
    Free(lists->boths_by_attr[i]);
    if(lists->directed) {
      Free(lists->tails_by_attr[i]);
      Free(lists->heads_by_attr[i]);
    }
  }
  Free(lists->boths_by_attr);
  if(lists->directed) {
    Free(lists->tails_by_attr);
    Free(lists->heads_by_attr);
  }
  Free(lists->boths);
  if(lists->directed) {
    Free(lists->tails);
    Free(lists->heads);
  }
  
  for(int i = 0; i < lists->bd_nlevels; i++) {
    Free(lists->bothpos[i]);
    if(lists->directed) {
      Free(lists->tailpos[i]);
      Free(lists->headpos[i]);
    }
  }
  Free(lists->bothpos);
  if(lists->directed) {
    Free(lists->tailpos);
    Free(lists->headpos);
  }
  
  Free(lists);
}

static inline void BDNodeListsToggleIf(Vertex tail, Vertex head, BDNodeLists *lists, int tailcondition, int headcondition) {
  if(tailcondition) {
    if(lists->directed && (lists->bothpos[lists->bd_vattr[head]][tail] || lists->headpos[lists->bd_vattr[head]][tail])) {
      NodeListToggle(lists->boths[tail][lists->bd_vattr[head]], tail);
      NodeListToggle(lists->heads[tail][lists->bd_vattr[head]], tail);
    } else {
      NodeListToggle(lists->tails[tail][lists->bd_vattr[head]], tail);      
    }
  }

  if(headcondition) {
    if(lists->directed && (lists->bothpos[lists->bd_vattr[tail]][head] || lists->tailpos[lists->bd_vattr[tail]][head])) {
      NodeListToggle(lists->boths[head][lists->bd_vattr[tail]], head);
      NodeListToggle(lists->tails[head][lists->bd_vattr[tail]], head);
    } else {
      NodeListToggle(lists->heads[head][lists->bd_vattr[tail]], head);        
    }
  }
}

#endif
