#ifndef _ERGM_CHANGESTAT_MULTILAYER_H_
#define _ERGM_CHANGESTAT_MULTILAYER_H_

#include "ergm_changestat.h"
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
  double *stacks;
} StoreLayerLogic;

#define ML_OI_TAIL(ll, l, t) ((Vertex) ((ll)->inwp->bipartite? (t) + ((l)-1)*(ll)->onwp->bipartite : (t) + ((l)-1)*(ll)->onwp->nnodes))
#define ML_OI_HEAD(ll, l, h) ((Vertex) ((ll)->inwp->bipartite? (h) + (ll)->inwp->bipartite - (ll)->onwp->bipartite + ((l)-1)*((ll)->onwp->nnodes-(ll)->onwp->bipartite) : (h) + ((l)-1)*(ll)->onwp->nnodes))

#define ML_IO_TAIL(ll, t) ((Vertex) ((ll)->lmap[t]))
#define ML_IO_HEAD(ll, h) ((Vertex) (ll)->lmap[h])
#define ML_LID_TAIL(ll, t) ((Vertex) ((ll)->lid[t]))
#define ML_LID_HEAD(ll, h) ((Vertex) ((ll)->lid[h]))

#define ML_IGETWT(ll, l,a,b) (GetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b), ll->inwp))
#define ML_ISETWT(ll, l,a,b,w) (SetEdge(ML_OI_TAIL(ll, l, a), ML_OI_HEAD(ll, l, b),w,(ll)->inwp))


/*
  Network logic language:
  
  com > 0: reference to layer; look up dyad value and push

  com == 0: numeric literal; push value of next command

  see lookup table in R/InitErgmTerm.multilayer.R
    pack.LayerLogic_formula_as_double() for com < 0

  change = 0 a.k.a. FALSE: Just evaluate as is.

  change = 1 a.k.a. TRUE: Difference between a hypothetical toggle of
    (tail,head) and its current state.
   
  change = 2: "encode": Instead of the difference return: a binary
    "encoding" of both: asis*1 + toggled*2

  change = 3: Return the post-toggle network.

 */


#define ergm_UNOP(op)				\
  {						\
    double x0 = *(stack0--);			\
    *(++stack0) = (op x0);			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      *(++stack1) = (op x1);			\
    }						\
    break;}

#define ergm_UNFUN(fun)				\
  {						\
    double x0 = *(stack0--);			\
    *(++stack0) = fun(x0);			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      *(++stack1) = fun(x1);			\
    }						\
    break;}


#define ergm_BINOP(op)				\
  {						\
    double x0 = *(stack0--);			\
    double y0 = *(stack0--);			\
    *(++stack0) = (x0 op y0);			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      double y1 = *(stack1--);			\
      *(++stack1) = (x1 op y1);			\
    }						\
    break;}

#define ergm_BINFUN(fun)			\
  {						\
    double x0 = *(stack0--);			\
    double y0 = *(stack0--);			\
    *(++stack0) = fun(x0, y0);			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      double y1 = *(stack1--);			\
      *(++stack1) = fun(x1, y1);		\
    }						\
    break;}

#define ergm_LUNOP(op)				\
  {						\
    double x0 = *(stack0--);			\
    *(++stack0) = (op (x0!=0));			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      *(++stack1) = (op (x1!=0));		\
    }						\
    break;}

#define ergm_LUNFUN(fun)			\
  {						\
    double x0 = *(stack0--);			\
    *(++stack0) = fun(x0!=0);			\
    if(stack1){					\
      double x1 = *(stack1--);			\
      *(++stack1) = fun(x1!=0);			\
    }						\
    break;}


#define ergm_LBINOP(op)				\
  {						\
    double x0 = *(stack0--);			\
    double y0 = *(stack0--);			\
    *(++stack0) = ((x0!=0) op (y0!=0));		\
    if(stack1){					\
      double x1 = *(stack1--);			\
      double y1 = *(stack1--);			\
      *(++stack1) = ((x1!=0) op (y1!=0));	\
    }						\
    break;}

#define ergm_LBINFUN(fun)			\
  {						\
    double x0 = *(stack0--);			\
    double y0 = *(stack0--);			\
    *(++stack0) = fun(x0!=0, y0!=0);		\
    if(stack1){					\
      double x1 = *(stack1--);			\
      double y1 = *(stack1--);			\
      *(++stack1) = fun(x1!=0, y1!=0);		\
    }						\
    break;}


#define ergm_FLOORDIV(x,y) floor(x/y)

#define ergm_FROUND(x) fround(x,0)

static inline int ergm_LayerLogic(Vertex tail, Vertex head, // Dyad to toggle on LHS network.
				  StoreLayerLogic *ll, // Layer Logic
				  unsigned int change
				  ){
  double *commands = ll->commands;
  unsigned int ncom = *(commands++);
  double *stack0=ll->stacks-1, *stack1=change? ll->stacks+ncom-1 : NULL; // stack0 and stack1 always point to the top element (if any)
  Vertex lt = ML_IO_TAIL(ll, tail), lh = ML_IO_HEAD(ll, head), tl = ML_LID_TAIL(ll, tail);

  for(unsigned int i=0; i<ncom; i++){
    int com = *(commands++);
    switch(com){
    case 0:{
      double x0 = *(commands++);
      *(++stack0) = x0;
      if(stack1){
	*(++stack1) = x0;
      }
      break;}
    case -1:ergm_LUNOP(!)
    case -2:ergm_LBINOP(&&)
    case -3:ergm_LBINOP(||)
    case -4:ergm_LBINFUN(XOR)
    case -5:ergm_BINOP(==)
    case -6:ergm_BINOP(!=)
    case -7:ergm_BINOP(<)
    case -8:ergm_BINOP(>)
    case -9:ergm_BINOP(<=)
    case -10:ergm_BINOP(>=)
    case -11:ergm_BINOP(+)
    case -12:ergm_BINOP(-)
    case -13:ergm_BINOP(*)
    case -14:ergm_BINOP(/)
    case -15:ergm_BINFUN(fmod)
    case -16:ergm_UNOP(-)
    case -17:ergm_UNFUN(fabs)
    case -18:ergm_BINFUN(pow)
    case -19:ergm_BINFUN(ergm_FLOORDIV)
    case -20:ergm_UNFUN(ergm_FROUND)
    case -21:ergm_BINFUN(fround)
    case -22:ergm_UNFUN(sign)
    default:{
      Vertex l = com; 
      unsigned int x0 = ML_IGETWT(ll, l, lt, lh);
      *(++stack0) = x0;
      if(stack1){
	unsigned int x1 = tl==l? !(x0!=0) : x0;
	*(++stack1) = x1;
      }
      break;}
    }
  }

  switch(change){
  case 1: return (int)(*stack1!=0) - (int)(*stack0!=0);
  case 2: return (*stack0!=0) | ((*stack1!=0)<<1);
  case 3: return (*stack1!=0);
  default: return (*stack0!=0);
  }
}

#undef ergm_UNOP
#undef ergm_UNFUN
#undef ergm_BINOP
#undef ergm_BINFUN
#undef ergm_FLOORDIV

#endif // _ERGM_CHANGESTAT_MULTILAYER_H_
