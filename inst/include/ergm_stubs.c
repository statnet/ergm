
#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_changestat_common.do_not_include_directly.h"
double my_choose(double n, int r){
static double (*fun)(double,int) = NULL;
if(fun==NULL) fun = (double (*)(double,int)) R_FindSymbol("my_choose", "ergm", NULL);
return fun(n,r);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_changestats_auxnet.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_dyad_hashmap.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_dyad_hashmap_utils.h"
void PrintDyadMapUInt(StoreDyadMapUInt *h){
static void (*fun)(StoreDyadMapUInt *) = NULL;
if(fun==NULL) fun = (void (*)(StoreDyadMapUInt *)) R_FindSymbol("PrintDyadMapUInt", "ergm", NULL);
fun(h);
}
void PrintDyadSet(StoreDyadSet *h){
static void (*fun)(StoreDyadSet *) = NULL;
if(fun==NULL) fun = (void (*)(StoreDyadSet *)) R_FindSymbol("PrintDyadSet", "ergm", NULL);
fun(h);
}
StoreDyadSet * NetworkToDyadSet(Network *nwp){
static StoreDyadSet * (*fun)(Network *) = NULL;
if(fun==NULL) fun = (StoreDyadSet * (*)(Network *)) R_FindSymbol("NetworkToDyadSet", "ergm", NULL);
return fun(nwp);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_edgelist.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_edgetree.h"
Network * NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network * (*fun)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network * (*)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitialize", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
void NetworkDestroy(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("NetworkDestroy", "ergm", NULL);
fun(nwp);
}
Network * NetworkInitializeD(double *tails, double *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network * (*fun)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network * (*)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitializeD", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
Network * NetworkCopy(Network *src){
static Network * (*fun)(Network *) = NULL;
if(fun==NULL) fun = (Network * (*)(Network *)) R_FindSymbol("NetworkCopy", "ergm", NULL);
return fun(src);
}
SEXP Network2Redgelist(Network *nwp){
static SEXP (*fun)(Network *) = NULL;
if(fun==NULL) fun = (SEXP (*)(Network *)) R_FindSymbol("Network2Redgelist", "ergm", NULL);
return fun(nwp);
}
Network * Redgelist2Network(SEXP elR, Rboolean empty){
static Network * (*fun)(SEXP,Rboolean) = NULL;
if(fun==NULL) fun = (Network * (*)(SEXP,Rboolean)) R_FindSymbol("Redgelist2Network", "ergm", NULL);
return fun(elR,empty);
}
void SetEdge(Vertex tail, Vertex head, unsigned int weight, Network *nwp){
static void (*fun)(Vertex,Vertex,unsigned int,Network *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,unsigned int,Network *)) R_FindSymbol("SetEdge", "ergm", NULL);
fun(tail,head,weight,nwp);
}
int ToggleEdge(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("ToggleEdge", "ergm", NULL);
return fun(tail,head,nwp);
}
int ToggleKnownEdge(Vertex tail, Vertex head, Network *nwp, Rboolean edgeflag){
static int (*fun)(Vertex,Vertex,Network *,Rboolean) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *,Rboolean)) R_FindSymbol("ToggleKnownEdge", "ergm", NULL);
return fun(tail,head,nwp,edgeflag);
}
int AddEdgeToTrees(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("AddEdgeToTrees", "ergm", NULL);
return fun(tail,head,nwp);
}
void AddHalfedgeToTree(Vertex a, Vertex b, TreeNode *edges, Edge *last_edge){
static void (*fun)(Vertex,Vertex,TreeNode *,Edge *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,TreeNode *,Edge *)) R_FindSymbol("AddHalfedgeToTree", "ergm", NULL);
fun(a,b,edges,last_edge);
}
void CheckEdgetreeFull(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("CheckEdgetreeFull", "ergm", NULL);
fun(nwp);
}
int DeleteEdgeFromTrees(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("DeleteEdgeFromTrees", "ergm", NULL);
return fun(tail,head,nwp);
}
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,Edge *last_edge){
static int (*fun)(Vertex,Vertex,TreeNode *,Edge *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,TreeNode *,Edge *)) R_FindSymbol("DeleteHalfedgeFromTree", "ergm", NULL);
return fun(a,b,edges,last_edge);
}
void RelocateHalfedge(Edge from, Edge to, TreeNode *edges){
static void (*fun)(Edge,Edge,TreeNode *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,Edge,TreeNode *)) R_FindSymbol("RelocateHalfedge", "ergm", NULL);
fun(from,to,edges);
}
void AddOnNetworkToggle(Network *nwp, OnNetworkToggle callback, void *payload, unsigned int pos){
static void (*fun)(Network *,OnNetworkToggle,void *,unsigned int) = NULL;
if(fun==NULL) fun = (void (*)(Network *,OnNetworkToggle,void *,unsigned int)) R_FindSymbol("AddOnNetworkToggle", "ergm", NULL);
fun(nwp,callback,payload,pos);
}
void DeleteOnNetworkToggle(Network *nwp, OnNetworkToggle callback, void *payload){
static void (*fun)(Network *,OnNetworkToggle,void *) = NULL;
if(fun==NULL) fun = (void (*)(Network *,OnNetworkToggle,void *)) R_FindSymbol("DeleteOnNetworkToggle", "ergm", NULL);
fun(nwp,callback,payload);
}
int FindithEdge(Vertex *tail, Vertex *head, Edge i, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Edge,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Edge,Network *)) R_FindSymbol("FindithEdge", "ergm", NULL);
return fun(tail,head,i,nwp);
}
int GetRandEdge(Vertex *tail, Vertex *head, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Network *)) R_FindSymbol("GetRandEdge", "ergm", NULL);
return fun(tail,head,nwp);
}
int FindithNondge(Vertex *tail, Vertex *head, Dyad i, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Dyad,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Dyad,Network *)) R_FindSymbol("FindithNondge", "ergm", NULL);
return fun(tail,head,i,nwp);
}
int GetRandNonedge(Vertex *tail, Vertex *head, Network *nwp){
static int (*fun)(Vertex *,Vertex *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Network *)) R_FindSymbol("GetRandNonedge", "ergm", NULL);
return fun(tail,head,nwp);
}
void printedge(Edge e, TreeNode *edges){
static void (*fun)(Edge,TreeNode *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,TreeNode *)) R_FindSymbol("printedge", "ergm", NULL);
fun(e,edges);
}
void InOrderTreeWalk(TreeNode *edges, Edge x){
static void (*fun)(TreeNode *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(TreeNode *,Edge)) R_FindSymbol("InOrderTreeWalk", "ergm", NULL);
fun(edges,x);
}
void NetworkEdgeList(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("NetworkEdgeList", "ergm", NULL);
fun(nwp);
}
void ShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
static void (*fun)(Vertex *,Vertex *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex *,Vertex *,Edge)) R_FindSymbol("ShuffleEdges", "ergm", NULL);
fun(tails,heads,nedges);
}
void DetShuffleEdges(Vertex *tails, Vertex *heads, Edge nedges){
static void (*fun)(Vertex *,Vertex *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex *,Vertex *,Edge)) R_FindSymbol("DetShuffleEdges", "ergm", NULL);
fun(tails,heads,nedges);
}
Edge DesignMissing(Vertex a, Vertex b, Network *mnwp){
static Edge (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex,Vertex,Network *)) R_FindSymbol("DesignMissing", "ergm", NULL);
return fun(a,b,mnwp);
}
Edge EdgeTree2EdgeList(Vertex *tails, Vertex *heads, Network *nwp, Edge nmax){
static Edge (*fun)(Vertex *,Vertex *,Network *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex *,Vertex *,Network *,Edge)) R_FindSymbol("EdgeTree2EdgeList", "ergm", NULL);
return fun(tails,heads,nwp,nmax);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_etamap.h"
void ergm_eta(double *theta, SEXP etamap, double *eta){
static void (*fun)(double *,SEXP,double *) = NULL;
if(fun==NULL) fun = (void (*)(double *,SEXP,double *)) R_FindSymbol("ergm_eta", "ergm", NULL);
fun(theta,etamap,eta);
}
void ergm_etagrad(double *theta, SEXP etamap, double *eta){
static void (*fun)(double *,SEXP,double *) = NULL;
if(fun==NULL) fun = (void (*)(double *,SEXP,double *)) R_FindSymbol("ergm_etagrad", "ergm", NULL);
fun(theta,etamap,eta);
}
void ergm_etagradmult(double *theta, double *v, unsigned int nv, SEXP etamap, double *ans){
static void (*fun)(double *,double *,unsigned int,SEXP,double *) = NULL;
if(fun==NULL) fun = (void (*)(double *,double *,unsigned int,SEXP,double *)) R_FindSymbol("ergm_etagradmult", "ergm", NULL);
fun(theta,v,nv,etamap,ans);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_MHproposal.h"
DegreeBound* DegreeBoundInitialize(int *attribs, int *maxout, int *maxin,int *minout, int *minin, int condAllDegExact,int attriblength, Network *nwp){
static DegreeBound* (*fun)(int *,int *,int *,int *,int *,int,int,Network *) = NULL;
if(fun==NULL) fun = (DegreeBound* (*)(int *,int *,int *,int *,int *,int,int,Network *)) R_FindSymbol("DegreeBoundInitialize", "ergm", NULL);
return fun(attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,nwp);
}
void DegreeBoundDestroy(DegreeBound *bd){
static void (*fun)(DegreeBound *) = NULL;
if(fun==NULL) fun = (void (*)(DegreeBound *)) R_FindSymbol("DegreeBoundDestroy", "ergm", NULL);
fun(bd);
}
MHProposal * MHProposalInitialize(SEXP pR, Network *nwp, void **aux_storage){
static MHProposal * (*fun)(SEXP,Network *,void **) = NULL;
if(fun==NULL) fun = (MHProposal * (*)(SEXP,Network *,void **)) R_FindSymbol("MHProposalInitialize", "ergm", NULL);
return fun(pR,nwp,aux_storage);
}
void MHProposalDestroy(MHProposal *MHp, Network *nwp){
static void (*fun)(MHProposal *,Network *) = NULL;
if(fun==NULL) fun = (void (*)(MHProposal *,Network *)) R_FindSymbol("MHProposalDestroy", "ergm", NULL);
fun(MHp,nwp);
}
int CheckTogglesValid(MHProposal *MHp, Network *nwp){
static int (*fun)(MHProposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHProposal *,Network *)) R_FindSymbol("CheckTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}
int CheckConstrainedTogglesValid(MHProposal *MHp, Network *nwp){
static int (*fun)(MHProposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHProposal *,Network *)) R_FindSymbol("CheckConstrainedTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_model.h"
Model* ModelInitialize(SEXP mR, SEXP ext_stateR, Network *nwp, Rboolean noinit_s){
static Model* (*fun)(SEXP,SEXP,Network *,Rboolean) = NULL;
if(fun==NULL) fun = (Model* (*)(SEXP,SEXP,Network *,Rboolean)) R_FindSymbol("ModelInitialize", "ergm", NULL);
return fun(mR,ext_stateR,nwp,noinit_s);
}
void ModelDestroy(Network *nwp, Model *m){
static void (*fun)(Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(Network *,Model *)) R_FindSymbol("ModelDestroy", "ergm", NULL);
fun(nwp,m);
}
int GetIndexForAttrValue(int value){
static int (*fun)(int) = NULL;
if(fun==NULL) fun = (int (*)(int)) R_FindSymbol("GetIndexForAttrValue", "ergm", NULL);
return fun(value);
}
void ChangeStats(unsigned int ntoggles, Vertex *tails, Vertex *heads, Network *nwp, Model *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,Network *,Model *)) R_FindSymbol("ChangeStats", "ergm", NULL);
fun(ntoggles,tails,heads,nwp,m);
}
void ChangeStats1(Vertex tail, Vertex head, Network *nwp, Model *m, Rboolean edgeflag){
static void (*fun)(Vertex,Vertex,Network *,Model *,Rboolean) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,Network *,Model *,Rboolean)) R_FindSymbol("ChangeStats1", "ergm", NULL);
fun(tail,head,nwp,m,edgeflag);
}
void ZStats(Network *nwp, Model *m, Rboolean skip_s){
static void (*fun)(Network *,Model *,Rboolean) = NULL;
if(fun==NULL) fun = (void (*)(Network *,Model *,Rboolean)) R_FindSymbol("ZStats", "ergm", NULL);
fun(nwp,m,skip_s);
}
void EmptyNetworkStats(Model *m, Rboolean skip_s){
static void (*fun)(Model *,Rboolean) = NULL;
if(fun==NULL) fun = (void (*)(Model *,Rboolean)) R_FindSymbol("EmptyNetworkStats", "ergm", NULL);
fun(m,skip_s);
}
void SummStats(Edge n_edges, Vertex *tails, Vertex *heads, Network *nwp, Model *m){
static void (*fun)(Edge,Vertex *,Vertex *,Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,Vertex *,Vertex *,Network *,Model *)) R_FindSymbol("SummStats", "ergm", NULL);
fun(n_edges,tails,heads,nwp,m);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_rlebdm.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_Rutil.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_state.h"
ErgmState * ErgmStateInit(SEXP stateR,Rboolean empty, Rboolean noinit_s){
static ErgmState * (*fun)(SEXP,Rboolean,Rboolean) = NULL;
if(fun==NULL) fun = (ErgmState * (*)(SEXP,Rboolean,Rboolean)) R_FindSymbol("ErgmStateInit", "ergm", NULL);
return fun(stateR,empty,noinit_s);
}
SEXP ErgmStateRSave(SEXP startR, ErgmState *s){
static SEXP (*fun)(SEXP,ErgmState *) = NULL;
if(fun==NULL) fun = (SEXP (*)(SEXP,ErgmState *)) R_FindSymbol("ErgmStateRSave", "ergm", NULL);
return fun(startR,s);
}
void ErgmStateDestroy(ErgmState *s){
static void (*fun)(ErgmState *) = NULL;
if(fun==NULL) fun = (void (*)(ErgmState *)) R_FindSymbol("ErgmStateDestroy", "ergm", NULL);
fun(s);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_wtedgetree.h"
WtNetwork * WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static WtNetwork * (*fun)(Vertex *,Vertex *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (WtNetwork * (*)(Vertex *,Vertex *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("WtNetworkInitialize", "ergm", NULL);
return fun(tails,heads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
void WtNetworkDestroy(WtNetwork *nwp){
static void (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *)) R_FindSymbol("WtNetworkDestroy", "ergm", NULL);
fun(nwp);
}
WtNetwork * WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static WtNetwork * (*fun)(double *,double *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (WtNetwork * (*)(double *,double *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("WtNetworkInitializeD", "ergm", NULL);
return fun(tails,heads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
WtNetwork * WtNetworkCopy(WtNetwork *src){
static WtNetwork * (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (WtNetwork * (*)(WtNetwork *)) R_FindSymbol("WtNetworkCopy", "ergm", NULL);
return fun(src);
}
SEXP WtNetwork2Redgelist(WtNetwork *nwp){
static SEXP (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (SEXP (*)(WtNetwork *)) R_FindSymbol("WtNetwork2Redgelist", "ergm", NULL);
return fun(nwp);
}
WtNetwork * Redgelist2WtNetwork(SEXP elR, Rboolean empty){
static WtNetwork * (*fun)(SEXP,Rboolean) = NULL;
if(fun==NULL) fun = (WtNetwork * (*)(SEXP,Rboolean)) R_FindSymbol("Redgelist2WtNetwork", "ergm", NULL);
return fun(elR,empty);
}
void WtSetEdge(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static void (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtSetEdge", "ergm", NULL);
fun(tail,head,weight,nwp);
}
int WtToggleEdge(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static int (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtToggleEdge", "ergm", NULL);
return fun(tail,head,weight,nwp);
}
int WtAddEdgeToTrees(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static int (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtAddEdgeToTrees", "ergm", NULL);
return fun(tail,head,weight,nwp);
}
void WtAddHalfedgeToTree(Vertex a, Vertex b, double weight, WtTreeNode *edges, Edge *last_edge){
static void (*fun)(Vertex,Vertex,double,WtTreeNode *,Edge *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,double,WtTreeNode *,Edge *)) R_FindSymbol("WtAddHalfedgeToTree", "ergm", NULL);
fun(a,b,weight,edges,last_edge);
}
void WtCheckEdgetreeFull(WtNetwork *nwp){
static void (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *)) R_FindSymbol("WtCheckEdgetreeFull", "ergm", NULL);
fun(nwp);
}
int WtDeleteEdgeFromTrees(Vertex tail, Vertex head, WtNetwork *nwp){
static int (*fun)(Vertex,Vertex,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,WtNetwork *)) R_FindSymbol("WtDeleteEdgeFromTrees", "ergm", NULL);
return fun(tail,head,nwp);
}
int WtDeleteHalfedgeFromTree(Vertex a, Vertex b, WtTreeNode *edges,Edge *last_edge){
static int (*fun)(Vertex,Vertex,WtTreeNode *,Edge *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,WtTreeNode *,Edge *)) R_FindSymbol("WtDeleteHalfedgeFromTree", "ergm", NULL);
return fun(a,b,edges,last_edge);
}
void WtRelocateHalfedge(Edge from, Edge to, WtTreeNode *edges){
static void (*fun)(Edge,Edge,WtTreeNode *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,Edge,WtTreeNode *)) R_FindSymbol("WtRelocateHalfedge", "ergm", NULL);
fun(from,to,edges);
}
void AddOnWtNetworkToggle(WtNetwork *nwp, OnWtNetworkToggle callback, void *payload, unsigned int pos){
static void (*fun)(WtNetwork *,OnWtNetworkToggle,void *,unsigned int) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,OnWtNetworkToggle,void *,unsigned int)) R_FindSymbol("AddOnWtNetworkToggle", "ergm", NULL);
fun(nwp,callback,payload,pos);
}
void DeleteOnWtNetworkToggle(WtNetwork *nwp, OnWtNetworkToggle callback, void *payload){
static void (*fun)(WtNetwork *,OnWtNetworkToggle,void *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,OnWtNetworkToggle,void *)) R_FindSymbol("DeleteOnWtNetworkToggle", "ergm", NULL);
fun(nwp,callback,payload);
}
int WtFindithEdge(Vertex *tail, Vertex *head, double *weight, Edge i, WtNetwork *nwp){
static int (*fun)(Vertex *,Vertex *,double *,Edge,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,double *,Edge,WtNetwork *)) R_FindSymbol("WtFindithEdge", "ergm", NULL);
return fun(tail,head,weight,i,nwp);
}
int WtGetRandEdge(Vertex *tail, Vertex *head, double *weight, WtNetwork *nwp){
static int (*fun)(Vertex *,Vertex *,double *,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,double *,WtNetwork *)) R_FindSymbol("WtGetRandEdge", "ergm", NULL);
return fun(tail,head,weight,nwp);
}
int WtFindithNonedge(Vertex *tail, Vertex *head, Dyad i, WtNetwork *nwp){
static int (*fun)(Vertex *,Vertex *,Dyad,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,Dyad,WtNetwork *)) R_FindSymbol("WtFindithNonedge", "ergm", NULL);
return fun(tail,head,i,nwp);
}
int WtGetRandNonedge(Vertex *tail, Vertex *head, WtNetwork *nwp){
static int (*fun)(Vertex *,Vertex *,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex *,Vertex *,WtNetwork *)) R_FindSymbol("WtGetRandNonedge", "ergm", NULL);
return fun(tail,head,nwp);
}
void Wtprintedge(Edge e, WtTreeNode *edges){
static void (*fun)(Edge,WtTreeNode *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,WtTreeNode *)) R_FindSymbol("Wtprintedge", "ergm", NULL);
fun(e,edges);
}
void WtInOrderTreeWalk(WtTreeNode *edges, Edge x){
static void (*fun)(WtTreeNode *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(WtTreeNode *,Edge)) R_FindSymbol("WtInOrderTreeWalk", "ergm", NULL);
fun(edges,x);
}
void WtNetworkEdgeList(WtNetwork *nwp){
static void (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *)) R_FindSymbol("WtNetworkEdgeList", "ergm", NULL);
fun(nwp);
}
void WtShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
static void (*fun)(Vertex *,Vertex *,double *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex *,Vertex *,double *,Edge)) R_FindSymbol("WtShuffleEdges", "ergm", NULL);
fun(tails,heads,weights,nedges);
}
void WtDetShuffleEdges(Vertex *tails, Vertex *heads, double *weights, Edge nedges){
static void (*fun)(Vertex *,Vertex *,double *,Edge) = NULL;
if(fun==NULL) fun = (void (*)(Vertex *,Vertex *,double *,Edge)) R_FindSymbol("WtDetShuffleEdges", "ergm", NULL);
fun(tails,heads,weights,nedges);
}
Edge WtDesignMissing(Vertex a, Vertex b, WtNetwork *mnwp){
static Edge (*fun)(Vertex,Vertex,WtNetwork *) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex,Vertex,WtNetwork *)) R_FindSymbol("WtDesignMissing", "ergm", NULL);
return fun(a,b,mnwp);
}
Edge WtEdgeTree2EdgeList(Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, Edge nmax){
static Edge (*fun)(Vertex *,Vertex *,double *,WtNetwork *,Edge) = NULL;
if(fun==NULL) fun = (Edge (*)(Vertex *,Vertex *,double *,WtNetwork *,Edge)) R_FindSymbol("WtEdgeTree2EdgeList", "ergm", NULL);
return fun(tails,heads,weights,nwp,nmax);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_wtMHproposal.h"
WtMHProposal * WtMHProposalInitialize(SEXP pR, WtNetwork *nwp, void **aux_storage){
static WtMHProposal * (*fun)(SEXP,WtNetwork *,void **) = NULL;
if(fun==NULL) fun = (WtMHProposal * (*)(SEXP,WtNetwork *,void **)) R_FindSymbol("WtMHProposalInitialize", "ergm", NULL);
return fun(pR,nwp,aux_storage);
}
void WtMHProposalDestroy(WtMHProposal *MH, WtNetwork *nwp){
static void (*fun)(WtMHProposal *,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtMHProposal *,WtNetwork *)) R_FindSymbol("WtMHProposalDestroy", "ergm", NULL);
fun(MH,nwp);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_wtmodel.h"
WtModel* WtModelInitialize(SEXP mR, SEXP ext_stateR, WtNetwork *nwp, Rboolean noinit_s){
static WtModel* (*fun)(SEXP,SEXP,WtNetwork *,Rboolean) = NULL;
if(fun==NULL) fun = (WtModel* (*)(SEXP,SEXP,WtNetwork *,Rboolean)) R_FindSymbol("WtModelInitialize", "ergm", NULL);
return fun(mR,ext_stateR,nwp,noinit_s);
}
void WtModelDestroy(WtNetwork *nwp, WtModel *m){
static void (*fun)(WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,WtModel *)) R_FindSymbol("WtModelDestroy", "ergm", NULL);
fun(nwp,m);
}
void WtChangeStats(unsigned int ntoggles, Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, WtModel *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *)) R_FindSymbol("WtChangeStats", "ergm", NULL);
fun(ntoggles,tails,heads,weights,nwp,m);
}
void WtChangeStats1(Vertex tail, Vertex head, double weight, WtNetwork *nwp, WtModel *m, double edgeweight){
static void (*fun)(Vertex,Vertex,double,WtNetwork *,WtModel *,double) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,double,WtNetwork *,WtModel *,double)) R_FindSymbol("WtChangeStats1", "ergm", NULL);
fun(tail,head,weight,nwp,m,edgeweight);
}
void WtZStats(WtNetwork *nwp, WtModel *m, Rboolean skip_s){
static void (*fun)(WtNetwork *,WtModel *,Rboolean) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,WtModel *,Rboolean)) R_FindSymbol("WtZStats", "ergm", NULL);
fun(nwp,m,skip_s);
}
void WtEmptyNetworkStats(WtModel *m, Rboolean skip_s){
static void (*fun)(WtModel *,Rboolean) = NULL;
if(fun==NULL) fun = (void (*)(WtModel *,Rboolean)) R_FindSymbol("WtEmptyNetworkStats", "ergm", NULL);
fun(m,skip_s);
}
void WtSummStats(Edge n_edges, Vertex *tails, Vertex *heads, double *weights, WtNetwork *nwp, WtModel *m){
static void (*fun)(Edge,Vertex *,Vertex *,double *,WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(Edge,Vertex *,Vertex *,double *,WtNetwork *,WtModel *)) R_FindSymbol("WtSummStats", "ergm", NULL);
fun(n_edges,tails,heads,weights,nwp,m);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_wtstate.h"
WtErgmState * WtErgmStateInit(SEXP stateR,Rboolean empty, Rboolean noinit_s){
static WtErgmState * (*fun)(SEXP,Rboolean,Rboolean) = NULL;
if(fun==NULL) fun = (WtErgmState * (*)(SEXP,Rboolean,Rboolean)) R_FindSymbol("WtErgmStateInit", "ergm", NULL);
return fun(stateR,empty,noinit_s);
}
SEXP WtErgmStateRSave(SEXP startR, WtErgmState *s){
static SEXP (*fun)(SEXP,WtErgmState *) = NULL;
if(fun==NULL) fun = (SEXP (*)(SEXP,WtErgmState *)) R_FindSymbol("WtErgmStateRSave", "ergm", NULL);
return fun(startR,s);
}
void WtErgmStateDestroy(WtErgmState *s){
static void (*fun)(WtErgmState *) = NULL;
if(fun==NULL) fun = (void (*)(WtErgmState *)) R_FindSymbol("WtErgmStateDestroy", "ergm", NULL);
fun(s);
}
