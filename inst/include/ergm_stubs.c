
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
#include "ergm_changestat_operator.h"
Model * unpack_Model_as_double(double **x, Network *nwp){
static Model * (*fun)(double **,Network *) = NULL;
if(fun==NULL) fun = (Model * (*)(double **,Network *)) R_FindSymbol("unpack_Model_as_double", "ergm", NULL);
return fun(x,nwp);
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
void SetEdge(Vertex tail, Vertex head, unsigned int weight, Network *nwp){
static void (*fun)(Vertex,Vertex,unsigned int,Network *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,unsigned int,Network *)) R_FindSymbol("SetEdge", "ergm", NULL);
fun(tail,head,weight,nwp);
}
void SetEdgeWithTimestamp(Vertex tail, Vertex head, unsigned int weight, Network *nwp){
static void (*fun)(Vertex,Vertex,unsigned int,Network *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,unsigned int,Network *)) R_FindSymbol("SetEdgeWithTimestamp", "ergm", NULL);
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
int ToggleEdgeWithTimestamp(Vertex tail, Vertex head, Network *nwp){
static int (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,Network *)) R_FindSymbol("ToggleEdgeWithTimestamp", "ergm", NULL);
return fun(tail,head,nwp);
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
void TouchEdge(Vertex tail, Vertex head, Network *nwp){
static void (*fun)(Vertex,Vertex,Network *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,Network *)) R_FindSymbol("TouchEdge", "ergm", NULL);
fun(tail,head,nwp);
}
void ExpireTimestamps(unsigned int edges, unsigned int nonedges, Network *nwp){
static void (*fun)(unsigned int,unsigned int,Network *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,unsigned int,Network *)) R_FindSymbol("ExpireTimestamps", "ergm", NULL);
fun(edges,nonedges,nwp);
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
MHProposal * MHProposalInitialize(const char *MHProposaltype, const char *MHProposalpackage,double *inputs,Network *nwp,int *attribs, int *maxout, int *maxin,int *minout, int *minin, int condAllDegExact,int attriblength,void **aux_storage){
static MHProposal * (*fun)(const char *,const char *,double *,Network *,int *,int *,int *,int *,int *,int,int,void **) = NULL;
if(fun==NULL) fun = (MHProposal * (*)(const char *,const char *,double *,Network *,int *,int *,int *,int *,int *,int,int,void **)) R_FindSymbol("MHProposalInitialize", "ergm", NULL);
return fun(MHProposaltype,MHProposalpackage,inputs,nwp,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,aux_storage);
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
Model* ModelInitialize(const char *fnames, const char *sonames, double **inputs,int n_terms, Network *nwp, Rboolean noinit_s){
static Model* (*fun)(const char *,const char *,double **,int,Network *,Rboolean) = NULL;
if(fun==NULL) fun = (Model* (*)(const char *,const char *,double **,int,Network *,Rboolean)) R_FindSymbol("ModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms,nwp,noinit_s);
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
void ChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, Network *nwp, Model *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,Network *,Model *)) R_FindSymbol("ChangeStats", "ergm", NULL);
fun(ntoggles,toggletail,togglehead,nwp,m);
}

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_rlebdm.h"

#define STUBFILE
#include <stddef.h>
#include <R_ext/Rdynload.h>
#include "ergm_wtchangestat_operator.h"
WtModel * unpack_WtModel_as_double(double **x, WtNetwork *nwp){
static WtModel * (*fun)(double **,WtNetwork *) = NULL;
if(fun==NULL) fun = (WtModel * (*)(double **,WtNetwork *)) R_FindSymbol("unpack_WtModel_as_double", "ergm", NULL);
return fun(x,nwp);
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
void WtSetEdge(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static void (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtSetEdge", "ergm", NULL);
fun(tail,head,weight,nwp);
}
void WtSetEdgeWithTimestamp(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static void (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtSetEdgeWithTimestamp", "ergm", NULL);
fun(tail,head,weight,nwp);
}
int WtToggleEdge(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static int (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtToggleEdge", "ergm", NULL);
return fun(tail,head,weight,nwp);
}
int WtToggleEdgeWithTimestamp(Vertex tail, Vertex head, double weight, WtNetwork *nwp){
static int (*fun)(Vertex,Vertex,double,WtNetwork *) = NULL;
if(fun==NULL) fun = (int (*)(Vertex,Vertex,double,WtNetwork *)) R_FindSymbol("WtToggleEdgeWithTimestamp", "ergm", NULL);
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
void WtTouchEdge(Vertex tail, Vertex head, WtNetwork *nwp){
static void (*fun)(Vertex,Vertex,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(Vertex,Vertex,WtNetwork *)) R_FindSymbol("WtTouchEdge", "ergm", NULL);
fun(tail,head,nwp);
}
void WtExpireTimestamps(unsigned int edges, unsigned int nonedges, WtNetwork *nwp){
static void (*fun)(unsigned int,unsigned int,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,unsigned int,WtNetwork *)) R_FindSymbol("WtExpireTimestamps", "ergm", NULL);
fun(edges,nonedges,nwp);
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
WtMHProposal * WtMHProposalInitialize(const char *MHProposaltype, const char *MHProposalpackage,double *inputs,WtNetwork *nwp,void **aux_storage){
static WtMHProposal * (*fun)(const char *,const char *,double *,WtNetwork *,void **) = NULL;
if(fun==NULL) fun = (WtMHProposal * (*)(const char *,const char *,double *,WtNetwork *,void **)) R_FindSymbol("WtMHProposalInitialize", "ergm", NULL);
return fun(MHProposaltype,MHProposalpackage,inputs,nwp,aux_storage);
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
WtModel* WtModelInitialize(const char *fnames, const char *sonames, double **inputs,int n_terms, WtNetwork *nwp, Rboolean noinit_s){
static WtModel* (*fun)(const char *,const char *,double **,int,WtNetwork *,Rboolean) = NULL;
if(fun==NULL) fun = (WtModel* (*)(const char *,const char *,double **,int,WtNetwork *,Rboolean)) R_FindSymbol("WtModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms,nwp,noinit_s);
}
void WtModelDestroy(WtNetwork *nwp, WtModel *m){
static void (*fun)(WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,WtModel *)) R_FindSymbol("WtModelDestroy", "ergm", NULL);
fun(nwp,m);
}
void WtChangeStats(unsigned int ntoggles, Vertex *toggletail, Vertex *togglehead, double *toggleweight, WtNetwork *nwp, WtModel *m){
static void (*fun)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(unsigned int,Vertex *,Vertex *,double *,WtNetwork *,WtModel *)) R_FindSymbol("WtChangeStats", "ergm", NULL);
fun(ntoggles,toggletail,togglehead,toggleweight,nwp,m);
}
