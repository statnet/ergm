#ifndef _ERGM_CHANGESTAT_MULTILAYER_H_
#define _ERGM_CHANGESTAT_MULTILAYER_H_

#include "ergm_edgetree.h"

#include "changestat_multilayer_common.inc"

/* layer-aware macros eponymous to ergm_changestat.h */
#define ML_IS_OUTEDGE(ll, a,b) (EdgetreeSearch((a),(b),(ll)->onwp->outedges)!=0?1:0)
#define ML_IS_INEDGE(ll, a,b) (EdgetreeSearch((a),(b),(ll)->onwp->inedges)!=0?1:0)
#define ML_IS_UNDIRECTED_EDGE(ll, a,b) ML_IS_OUTEDGE((ll), MIN(a,b), MAX(a,b))
#define ML_MIN_OUTEDGE(ll, a) (EdgetreeMinimum((ll)->onwp->outedges, (a)))
#define ML_MIN_INEDGE(ll, a) (EdgetreeMinimum((ll)->onwp->inedges, (a)))
#define ML_NEXT_OUTEDGE(ll, e) (EdgetreeSuccessor((ll)->onwp->outedges,(e)))
#define ML_NEXT_INEDGE(ll, e) (EdgetreeSuccessor((ll)->onwp->inedges,(e)))
#define ML_STEP_THROUGH_OUTEDGES(ll, a,e,v) for((e)=ML_MIN_OUTEDGE((ll), a);((v)=ML_OUTVAL((ll), e))!=0;(e)=ML_NEXT_OUTEDGE((ll), e))
#define ML_STEP_THROUGH_INEDGES(ll, a,e,v) for((e)=ML_MIN_INEDGE((ll), a);((v)=ML_INVAL((ll), e))!=0;(e)=ML_NEXT_INEDGE((ll), e))
#define ML_STEP_THROUGH_OUTEDGES_DECL(ll, a,e,v) for(Edge e=ML_MIN_OUTEDGE((ll), a);ML_OUTVAL((ll), e)!=0;e=ML_NEXT_OUTEDGE((ll), e))
#define ML_STEP_THROUGH_INEDGES_DECL(ll, a,e,v) for(Edge e=ML_MIN_INEDGE((ll), a);ML_INVAL((ll), e)!=0;e=ML_NEXT_INEDGE((ll), e))
#define ML_EXEC_THROUGH_OUTEDGES(ll, a,e,v,subroutine) if(ML_DIRECTED((ll))){ ML_STEP_THROUGH_OUTEDGES_DECL((ll), a,e,v) {Vertex v=ML_OUTVAL((ll), e); subroutine} } else { ML_EXEC_THROUGH_EDGES((ll), a,e,v,subroutine) }
#define ML_EXEC_THROUGH_INEDGES(ll, a,e,v,subroutine) if(ML_DIRECTED((ll))){ ML_STEP_THROUGH_INEDGES_DECL((ll), a,e,v) {Vertex v=ML_INVAL((ll), e); subroutine} } else { ML_EXEC_THROUGH_EDGES((ll), a,e,v,subroutine) }
#define ML_EXEC_THROUGH_EDGES(ll, a,e,v,subroutine) { ML_STEP_THROUGH_OUTEDGES_DECL((ll), a,e,v) {Vertex v=ML_OUTVAL((ll), e); subroutine}  ML_STEP_THROUGH_INEDGES_DECL((ll), a,e,v) {Vertex v=ML_INVAL((ll), e); subroutine} }
#define ML_EXEC_THROUGH_FOUTEDGES(ll, a,e,v,subroutine) ML_STEP_THROUGH_OUTEDGES_DECL((ll), a,e,v) {Vertex v=ML_OUTVAL((ll), e); subroutine}
#define ML_EXEC_THROUGH_FINEDGES(ll, a,e,v,subroutine) ML_STEP_THROUGH_INEDGES_DECL((ll), a,e,v) {Vertex v=ML_INVAL((ll), e); subroutine}
#define ML_EXEC_THROUGH_NET_EDGES(ll, a,b,e,subroutine) for(Vertex a=1; a <= N_NODES; a++)  ML_EXEC_THROUGH_FOUTEDGES((ll), a, e, b, {subroutine});
#define ML_TOGGLE(ll, a,b) (ToggleEdge((a),(b),(ll)->onwp));
#define ML_TOGGLE_DISCORD(ll, a,b) (ToggleEdge((a),(b),(ll)->onwp+1));
#define ML_GETWT(ll, a,b) (GetEdge(a,b,(ll)->onwp))
#define ML_SETWT(ll, a,b,w) (SetEdge(a,b,w,(ll)->onwp))

typedef struct {
  unsigned int nl;
  Network *inwp, *onwp;
  double *lid;
  double *lmap;
  double *commands;
  unsigned int *stacks;
} StoreLayerLogic;

#define ML_OI_TAIL(ll, l, t) ((ll)->inwp->bipartite? (t) + ((l)-1)*(ll)->onwp->bipartite : (t) + ((l)-1)*(ll)->onwp->nnodes)
#define ML_OI_HEAD(ll, l, h) ((ll)->inwp->bipartite? (h) + (ll)->inwp->bipartite - (ll)->onwp->bipartite + ((l)-1)*((ll)->onwp->nnodes-(ll)->onwp->bipartite) : (h) + ((l)-1)*(ll)->onwp->nnodes)

#define ML_IO_TAIL(ll, t) ((ll)->lmap[t])
#define ML_IO_HEAD(ll, h) ((ll)->lmap[h])
#define ML_LID_TAIL(ll, t) ((ll)->lid[t])
#define ML_LID_HEAD(ll, h) ((ll)->lid[h])

#define ML_IGETWT(ll, l,a,b) (GetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b), ll->inwp))
#define ML_ISETWT(ll, l,a,b,w) (SetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b),w,(ll)->inwp))


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
  Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head), tl = ML_LID_TAIL(ll, tail);

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
      unsigned int x0 = ML_IGETWT(ll, l, lt, lh);
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

#endif // _ERGM_CHANGESTAT_MULTILAYER_H_
