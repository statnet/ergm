#include "ergm_edgetree.h"
#include "ergm_changestat.h"
#include "ergm_changestat_operator.h"

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

  if(stack1) return (int)*stack1 - (int)*stack0;
  else return *stack0;
}

I_CHANGESTAT_FN(i__layer_net){
  double *inputs = INPUT_PARAM;
  ALLOC_AUX_STORAGE(1, StoreLayerLogic, ll); inputs++;
  ll->nl = *(inputs++);
  ll->inwp = nwp;
  ll->onwp = Calloc(1, Network);
  ll->lid = inputs - 1; // The -1 is because Vertex IDs count from 1.
  inputs += N_NODES;
  ll->lmap = inputs - 1;
  inputs += N_NODES;
  ll->commands = inputs;
  ll->stacks = Calloc(2*ll->commands[0], unsigned int);

  Vertex lnnodes = N_NODES/ll->nl, lbip = BIPARTITE/ll->nl;
  ll->onwp[0] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      if(ergm_LayerLogic(t, h, ll, TRUE)){
	ToggleEdge(ll->lmap[t], ll->lmap[h], ll->onwp);
      }
    });
}

U_CHANGESTAT_FN(u__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  if(ergm_LayerLogic(tail, head, ll, TRUE)){
    ToggleEdge(ll->lmap[tail], ll->lmap[head], ll->onwp);
  }
}

F_CHANGESTAT_FN(f__layer_net){ 
  GET_AUX_STORAGE(StoreLayerLogic, ll);
  NetworkDestroy(ll->onwp);
  Free(ll->onwp);
  Free(ll->stacks);
}

I_CHANGESTAT_FN(i__layer_nets){
  ALLOC_AUX_STORAGE(1, StoreNetsAndLIDAndLMapAndNL, li);
  li->nl = INPUT_PARAM[1];
  Vertex lnnodes = N_NODES/li->nl, lbip = BIPARTITE/li->nl;
  li->nwp = Calloc(li->nl+1, Network);
  li->lid = INPUT_PARAM+2 -1; // The -1 is because Vertex IDs count from 1.
  li->lmap = INPUT_PARAM+2+N_NODES -1;
  for(unsigned int l = 1; l <= li->nl; l++){
    li->nwp[l] = NetworkInitialize(NULL, NULL, 0, lnnodes, DIRECTED, lbip, 0, 0, NULL);
  }
  
  EXEC_THROUGH_NET_EDGES(t, h, e, {
      ToggleEdge(li->lmap[t], li->lmap[h], li->nwp + (int)li->lid[t]);
    });
}

U_CHANGESTAT_FN(u__layer_nets){ 
  GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li);
  ToggleEdge(li->lmap[tail], li->lmap[head], li->nwp + (int)li->lid[tail]);
}

F_CHANGESTAT_FN(f__layer_nets){ 
  GET_AUX_STORAGE(StoreNetsAndLIDAndLMapAndNL, li);
  for(unsigned int l = 1; l <= li->nl; l++){
    NetworkDestroy(li->nwp + l);
  }
  Free(li->nwp);
}

I_CHANGESTAT_FN(i_OnLayer){
  
  unsigned int nml = *INPUT_ATTRIB; // Number of layers *in the term*. Already shifted past the auxiliaries.
  
  ALLOC_STORAGE(nml, Model*, ms);

  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    double *inputs = INPUT_ATTRIB+1; // Rewind to the start of model spec.
    ms[ml] = unpack_Model_as_double(&inputs);
    InitStats(ll->onwp, ms[ml]);
  }
}

C_CHANGESTAT_FN(c_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
      Vertex lt = ll->lmap[tail], lh = ll->lmap[head];
      ChangeStats(1, &lt, &lh, ll->onwp, ms[ml]);
      for(unsigned int i=0; i<N_CHANGE_STATS; i++)
	CHANGE_STAT[i] += ms[ml]->workspace[i];
    }
  }
}

U_CHANGESTAT_FN(u_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;

  // Find the affected models.
  for(unsigned int ml=0; ml < nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    if(ergm_LayerLogic(tail, head, ll, TRUE)){ // network affected
      Vertex lt = ll->lmap[tail], lh = ll->lmap[head];
      Model *m = ms[ml];
      UPDATE_STORAGE(lt, lh, ll->onwp, ms[ml], NULL);
    }
  }
}

F_CHANGESTAT_FN(f_OnLayer){
  GET_STORAGE(Model*, ms);
  unsigned int nml = *INPUT_ATTRIB;
  for(unsigned int ml=0; ml<nml; ml++){
    GET_AUX_STORAGE_NUM(StoreLayerLogic, ll, ml);
    ModelDestroy(ll->onwp, ms[ml]);
  }
}
