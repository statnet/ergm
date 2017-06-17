#ifndef _ERGM_CHANGESTAT_LAYERS_H_
#define _ERGM_CHANGESTAT_LAYERS_H_

#include "ergm_edgetree.h"

typedef struct {
  unsigned int nl;
  Network *nwp;
  double *lid;
  double *lmap;
} StoreNetsAndLIDAndLMapAndNL;

typedef struct {
  unsigned int nl;
  Network *inwp, *onwp;
  double *lid;
  double *lmap;
  double *commands;
  unsigned int *stacks;
} StoreLayerLogic;

#define GLOBAL_TAIL(t, l, ll) ((ll)->inwp->bipartite? (t) + ((l)-1)*(ll)->onwp->bipartite : (t) + ((l)-1)*(ll)->onwp->nnodes)
#define GLOBAL_HEAD(h, l, ll) ((ll)->inwp->bipartite? (h) + (ll)->inwp->bipartite - (ll)->onwp->bipartite + ((l)-1)*((ll)->onwp->nnodes-(ll)->onwp->bipartite) : (h) + ((l)-1)*(ll)->onwp->nnodes)

/*
  Network logic language:
  
  com > 0: reference to layer; look up dyad value and push
  
  com == -1: negation; pop, invert, push
  
  com == -2: logical and; pop 2, take conjunction, push

  com == -3: logical or; pop 2, take disjunction, push

  com == -4: test for equality; pop 2, check equality, push

  com == -5: test for inequality (xor); pop 2, check inequality, push

  change = 0 a.k.a. FALSE: Just evaluate as is.

  change = 1 a.k.a. TRUE: Difference between a hypothetical toggle of
    (tail,head) and its current state.
   
  change = 2: "encode": Instead of the difference return: a binary "encoding" of both:
    asis*1 + toggled*2

  change = 3: Return the post-toggle network.

 */

static inline int ergm_LayerLogic(Vertex tail, Vertex head, // Dyad to toggle on LHS network.
				  StoreLayerLogic *ll, // Layer Logic
				  unsigned int change
				  ){
  double *commands = ll->commands;
  unsigned int ncom = *(commands++);
  unsigned int *stack0=ll->stacks-1, *stack1=change? ll->stacks+ncom-1 : NULL; // stack0 and stack1 always point to the top element (if any)
  Vertex lt = ll->lmap[tail], lh = ll->lmap[head], tl = ll->lid[tail];

  for(unsigned int i=0; i<ncom; i++){
    int com = *(commands++);
    switch(com){
    case -1:{
      unsigned int x0 = *(stack0--);
      *(++stack0) = !x0;
      if(stack1){
	unsigned int x1 = *(stack1--);
	*(++stack1) = !x1;
      }
      break;}
    case -2:{
      unsigned int x0 = *(stack0--);
      unsigned int y0 = *(stack0--);
      *(++stack0) = x0 && y0;
      if(stack1){
	unsigned int x1 = *(stack1--);
	unsigned int y1 = *(stack1--);
	*(++stack1) = x1 && y1;
      }
      break;}
    case -3:{
      unsigned int x0 = *(stack0--);
      unsigned int y0 = *(stack0--);
      *(++stack0) = x0 || y0;
      if(stack1){
	unsigned int x1 = *(stack1--);
	unsigned int y1 = *(stack1--);
	*(++stack1) = x1 || y1;
      }
      break;}
    case -4:{
      unsigned int x0 = *(stack0--);
      unsigned int y0 = *(stack0--);
      *(++stack0) = x0 == y0;
      if(stack1){
	unsigned int x1 = *(stack1--);
	unsigned int y1 = *(stack1--);
	*(++stack1) = x1 == y1;
      }
      break;}
    case -5:{
      unsigned int x0 = *(stack0--);
      unsigned int y0 = *(stack0--);
      *(++stack0) = x0 != y0;
      if(stack1){
	unsigned int x1 = *(stack1--);
	unsigned int y1 = *(stack1--);
	*(++stack1) = x1 != y1;
      }
      break;}
    default:{
      Vertex l = com; 
      unsigned int x0 = GetEdge(GLOBAL_TAIL(lt, l, ll), GLOBAL_HEAD(lh, l, ll), ll->inwp);
      *(++stack0) = x0;
      if(stack1){
	unsigned int x1 = tl==l? !x0 : x0;
	*(++stack1) = x1;
      }
      break;}
    }
  }

  switch(change){
  case 1: return (int)*stack1 - (int)*stack0;
  case 2: return *stack0 | (*stack1<<1);
  case 3: return *stack1;
  default: return *stack0;
  }
}

#endif // _ERGM_CHANGESTAT_LAYERS_H_
