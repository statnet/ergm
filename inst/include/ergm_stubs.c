#include "R_ext/Rdynload.h"
#include "ergm_changestat_common.do_not_include_directly.h"
double my_choose(double n, int r){
static double (*fun)(double,int) = NULL;
if(fun==NULL) fun = (double (*)(double,int)) R_FindSymbol("my_choose", "ergm", NULL);
return fun(n,r);
}
#undef CHOOSE
#undef EQUAL
#undef XOR
#undef XNOR
#undef OUTVAL
#undef INVAL
#undef N_NODES
#undef N_DYADS
#undef OUT_DEG
#undef IN_DEG
#undef DIRECTED
#undef N_EDGES
#undef BIPARTITE
#undef N_TAILS
#undef N_HEADS
#undef NEXT_INEDGE_NUM
#undef NEXT_OUTEDGE_NUM
#undef CHANGE_STAT
#undef N_CHANGE_STATS
#undef INPUT_PARAM
#undef N_INPUT_PARAMS
#undef ZERO_ALL_CHANGESTATS
#undef INPUT_ATTRIB
#include "R_ext/Rdynload.h"
#include "ergm_changestat_operator.h"
Model * unpack_Model_as_double(double **x){
static Model * (*fun)(double **) = NULL;
if(fun==NULL) fun = (Model * (*)(double **)) R_FindSymbol("unpack_Model_as_double", "ergm", NULL);
return fun(x);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef CHOOSE
#undef EQUAL
#undef XOR
#undef XNOR
#undef OUTVAL
#undef INVAL
#undef N_NODES
#undef N_DYADS
#undef OUT_DEG
#undef IN_DEG
#undef DIRECTED
#undef N_EDGES
#undef BIPARTITE
#undef N_TAILS
#undef N_HEADS
#undef NEXT_INEDGE_NUM
#undef NEXT_OUTEDGE_NUM
#undef CHANGE_STAT
#undef N_CHANGE_STATS
#undef INPUT_PARAM
#undef N_INPUT_PARAMS
#undef ZERO_ALL_CHANGESTATS
#undef INPUT_ATTRIB
#undef TAIL
#undef HEAD
#undef IS_OUTEDGE
#undef IS_INEDGE
#undef IS_UNDIRECTED_EDGE
#undef MIN_OUTEDGE
#undef MIN_INEDGE
#undef NEXT_OUTEDGE
#undef NEXT_INEDGE
#undef NEXT_OUTEDGE_PRE
#undef NEXT_INEDGE_PRE
#undef STEP_THROUGH_OUTEDGES
#undef STEP_THROUGH_INEDGES
#undef STEP_THROUGH_OUTEDGES_PRE
#undef STEP_THROUGH_INEDGES_PRE
#undef STEP_THROUGH_OUTEDGES_DECL
#undef STEP_THROUGH_INEDGES_DECL
#undef STEP_THROUGH_OUTEDGES_PRE_DECL
#undef STEP_THROUGH_INEDGES_PRE_DECL
#undef EXEC_THROUGH_OUTEDGES
#undef EXEC_THROUGH_INEDGES
#undef EXEC_THROUGH_EDGES
#undef EXEC_THROUGH_OUTEDGES_PRE
#undef EXEC_THROUGH_INEDGES_PRE
#undef EXEC_THROUGH_EDGES_PRE
#undef EXEC_THROUGH_FOUTEDGES
#undef EXEC_THROUGH_FINEDGES
#undef EXEC_THROUGH_FOUTEDGES_PRE
#undef EXEC_THROUGH_FINEDGES_PRE
#undef EXEC_THROUGH_NET_EDGES
#undef EXEC_THROUGH_NET_EDGES_PRE
#undef TOGGLE
#undef TOGGLE_DISCORD
#undef GETWT
#undef SETWT
#undef FOR_EACH_TOGGLE
#undef IF_MORE_TO_COME
#undef TOGGLE_IF_MORE_TO_COME
#undef TOGGLE_DISCORD_IF_MORE_TO_COME
#undef UNDO_PREVIOUS
#undef UNDO_PREVIOUS_TOGGLES
#undef UNDO_PREVIOUS_DISCORD_TOGGLES
#undef CHANGESTAT_FN
#undef C_CHANGESTAT_FN
#undef D_CHANGESTAT_FN
#undef I_CHANGESTAT_FN
#undef U_CHANGESTAT_FN
#undef F_CHANGESTAT_FN
#undef S_CHANGESTAT_FN
#undef D_FROM_S
#undef D_FROM_S_FN
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef BD_LOOP
#undef BD_COND_LOOP
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef MH_I_FN
#undef MH_U_FN
#undef MH_P_FN
#undef MH_F_FN
#undef FOR_EACH_TERM
#undef EXEC_THROUGH_TERMS
#undef FOR_EACH_TERM_INREVERSE
#undef EXEC_THROUGH_TERMS_INREVERSE
#undef EXEC_THROUGH_TERMS_INTO
#undef IFDEBUG_BACKUP_DSTATS
#undef IFDEBUG_RESTORE_DSTATS
#undef UPDATE_STORAGE_COND
#undef UPDATE_STORAGE
#undef STORAGE
#undef ALLOC_STORAGE
#undef GET_STORAGE
#undef AUX_STORAGE
#undef ALLOC_AUX_STORAGE
#undef GET_AUX_STORAGE
#undef AUX_STORAGE_NUM
#undef GET_AUX_STORAGE_NUM
#undef ALLOC_AUX_SOCIOMATRIX
#undef FREE_AUX_SOCIOMATRIX
#include "R_ext/Rdynload.h"
#include "ergm_dyad_hashmap.h"
void PrintDyadMapUInt(StoreDyadMapUInt *h){
static void (*fun)(StoreDyadMapUInt *) = NULL;
if(fun==NULL) fun = (void (*)(StoreDyadMapUInt *)) R_FindSymbol("PrintDyadMapUInt", "ergm", NULL);
fun(h);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef kh_inline
#undef klib_unused
#undef __ac_isempty
#undef __ac_isdel
#undef __ac_iseither
#undef __ac_set_isdel_false
#undef __ac_set_isempty_false
#undef __ac_set_isboth_false
#undef __ac_set_isdel_true
#undef __ac_fsize
#undef kroundup32
#undef kcalloc
#undef kmalloc
#undef krealloc
#undef kfree
#undef __KHASH_TYPE
#undef __KHASH_PROTOTYPES
#undef __KHASH_IMPL
#undef KHASH_DECLARE
#undef KHASH_INIT2
#undef KHASH_INIT
#undef kh_int_hash_func
#undef kh_int_hash_equal
#undef kh_int64_hash_func
#undef kh_int64_hash_equal
#undef kh_str_hash_func
#undef kh_str_hash_equal
#undef kh_int_hash_func2
#undef khash_t
#undef kh_init
#undef kh_destroy
#undef kh_clear
#undef kh_resize
#undef kh_put
#undef kh_get
#undef kh_getval
#undef kh_del
#undef kh_exist
#undef kh_key
#undef kh_val
#undef kh_value
#undef kh_begin
#undef kh_end
#undef kh_size
#undef kh_n_buckets
#undef kh_foreach
#undef kh_foreach_value
#undef KHASH_SET_INIT_INT
#undef KHASH_MAP_INIT_INT
#undef KHASH_SET_INIT_INT64
#undef KHASH_MAP_INIT_INT64
#undef KHASH_SET_INIT_STR
#undef KHASH_MAP_INIT_STR
#undef kh_set
#undef kh_unset
#undef TH
#undef THD
#undef THU
#undef ROT_INT
#undef kh_vertexvertex_hash_func
#undef kh_vertexvertex_hash_equal
#include "R_ext/Rdynload.h"
#include "ergm_edgelist.h"
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#include "R_ext/Rdynload.h"
#include "ergm_edgetree.h"
Network NetworkInitialize(Vertex *tails, Vertex *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network (*fun)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network (*)(Vertex *,Vertex *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitialize", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
void NetworkDestroy(Network *nwp){
static void (*fun)(Network *) = NULL;
if(fun==NULL) fun = (void (*)(Network *)) R_FindSymbol("NetworkDestroy", "ergm", NULL);
fun(nwp);
}
Network NetworkInitializeD(double *tails, double *heads, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static Network (*fun)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (Network (*)(double *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("NetworkInitializeD", "ergm", NULL);
return fun(tails,heads,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
Network * NetworkCopy(Network *dest, Network *src){
static Network * (*fun)(Network *,Network *) = NULL;
if(fun==NULL) fun = (Network * (*)(Network *,Network *)) R_FindSymbol("NetworkCopy", "ergm", NULL);
return fun(dest,src);
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
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#include "R_ext/Rdynload.h"
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
void MH_init(MHproposal *MHp,char *MHproposaltype, char *MHproposalpackage,double *inputs,int fVerbose,Network *nwp,int *attribs, int *maxout, int *maxin,int *minout, int *minin, int condAllDegExact,int attriblength,void **aux_storage){
static void (*fun)(MHproposal *,char *,char *,double *,int,Network *,int *,int *,int *,int *,int *,int,int,void **) = NULL;
if(fun==NULL) fun = (void (*)(MHproposal *,char *,char *,double *,int,Network *,int *,int *,int *,int *,int *,int,int,void **)) R_FindSymbol("MH_init", "ergm", NULL);
fun(MHp,MHproposaltype,MHproposalpackage,inputs,fVerbose,nwp,attribs,maxout,maxin,minout,minin,condAllDegExact,attriblength,aux_storage);
}
void MH_free(MHproposal *MHp, Network *nwp){
static void (*fun)(MHproposal *,Network *) = NULL;
if(fun==NULL) fun = (void (*)(MHproposal *,Network *)) R_FindSymbol("MH_free", "ergm", NULL);
fun(MHp,nwp);
}
int CheckTogglesValid(MHproposal *MHp, Network *nwp){
static int (*fun)(MHproposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHproposal *,Network *)) R_FindSymbol("CheckTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}
int CheckConstrainedTogglesValid(MHproposal *MHp, Network *nwp){
static int (*fun)(MHproposal *,Network *) = NULL;
if(fun==NULL) fun = (int (*)(MHproposal *,Network *)) R_FindSymbol("CheckConstrainedTogglesValid", "ergm", NULL);
return fun(MHp,nwp);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef BD_LOOP
#undef BD_COND_LOOP
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef MH_I_FN
#undef MH_U_FN
#undef MH_P_FN
#undef MH_F_FN
#include "R_ext/Rdynload.h"
#include "ergm_model.h"
Model* ModelInitialize(char *fnames, char *sonames, double **inputs,int n_terms){
static Model* (*fun)(char *,char *,double **,int) = NULL;
if(fun==NULL) fun = (Model* (*)(char *,char *,double **,int)) R_FindSymbol("ModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms);
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
void InitStats(Network *nwp, Model *m){
static void (*fun)(Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(Network *,Model *)) R_FindSymbol("InitStats", "ergm", NULL);
fun(nwp,m);
}
void DestroyStats(Network *nwp, Model *m){
static void (*fun)(Network *,Model *) = NULL;
if(fun==NULL) fun = (void (*)(Network *,Model *)) R_FindSymbol("DestroyStats", "ergm", NULL);
fun(nwp,m);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef CHOOSE
#undef EQUAL
#undef XOR
#undef XNOR
#undef OUTVAL
#undef INVAL
#undef N_NODES
#undef N_DYADS
#undef OUT_DEG
#undef IN_DEG
#undef DIRECTED
#undef N_EDGES
#undef BIPARTITE
#undef N_TAILS
#undef N_HEADS
#undef NEXT_INEDGE_NUM
#undef NEXT_OUTEDGE_NUM
#undef CHANGE_STAT
#undef N_CHANGE_STATS
#undef INPUT_PARAM
#undef N_INPUT_PARAMS
#undef ZERO_ALL_CHANGESTATS
#undef INPUT_ATTRIB
#undef TAIL
#undef HEAD
#undef IS_OUTEDGE
#undef IS_INEDGE
#undef IS_UNDIRECTED_EDGE
#undef MIN_OUTEDGE
#undef MIN_INEDGE
#undef NEXT_OUTEDGE
#undef NEXT_INEDGE
#undef NEXT_OUTEDGE_PRE
#undef NEXT_INEDGE_PRE
#undef STEP_THROUGH_OUTEDGES
#undef STEP_THROUGH_INEDGES
#undef STEP_THROUGH_OUTEDGES_PRE
#undef STEP_THROUGH_INEDGES_PRE
#undef STEP_THROUGH_OUTEDGES_DECL
#undef STEP_THROUGH_INEDGES_DECL
#undef STEP_THROUGH_OUTEDGES_PRE_DECL
#undef STEP_THROUGH_INEDGES_PRE_DECL
#undef EXEC_THROUGH_OUTEDGES
#undef EXEC_THROUGH_INEDGES
#undef EXEC_THROUGH_EDGES
#undef EXEC_THROUGH_OUTEDGES_PRE
#undef EXEC_THROUGH_INEDGES_PRE
#undef EXEC_THROUGH_EDGES_PRE
#undef EXEC_THROUGH_FOUTEDGES
#undef EXEC_THROUGH_FINEDGES
#undef EXEC_THROUGH_FOUTEDGES_PRE
#undef EXEC_THROUGH_FINEDGES_PRE
#undef EXEC_THROUGH_NET_EDGES
#undef EXEC_THROUGH_NET_EDGES_PRE
#undef TOGGLE
#undef TOGGLE_DISCORD
#undef GETWT
#undef SETWT
#undef FOR_EACH_TOGGLE
#undef IF_MORE_TO_COME
#undef TOGGLE_IF_MORE_TO_COME
#undef TOGGLE_DISCORD_IF_MORE_TO_COME
#undef UNDO_PREVIOUS
#undef UNDO_PREVIOUS_TOGGLES
#undef UNDO_PREVIOUS_DISCORD_TOGGLES
#undef CHANGESTAT_FN
#undef C_CHANGESTAT_FN
#undef D_CHANGESTAT_FN
#undef I_CHANGESTAT_FN
#undef U_CHANGESTAT_FN
#undef F_CHANGESTAT_FN
#undef S_CHANGESTAT_FN
#undef D_FROM_S
#undef D_FROM_S_FN
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef BD_LOOP
#undef BD_COND_LOOP
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef MH_I_FN
#undef MH_U_FN
#undef MH_P_FN
#undef MH_F_FN
#undef FOR_EACH_TERM
#undef EXEC_THROUGH_TERMS
#undef FOR_EACH_TERM_INREVERSE
#undef EXEC_THROUGH_TERMS_INREVERSE
#undef EXEC_THROUGH_TERMS_INTO
#undef IFDEBUG_BACKUP_DSTATS
#undef IFDEBUG_RESTORE_DSTATS
#undef UPDATE_STORAGE_COND
#undef UPDATE_STORAGE
#include "R_ext/Rdynload.h"
#include "ergm_rlebdm.h"
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef TH2Dyad
#undef Dyad2T
#undef Dyad2H
#include "R_ext/Rdynload.h"
#include "ergm_wtchangestat_operator.h"
WtModel * unpack_WtModel_as_double(double **x){
static WtModel * (*fun)(double **) = NULL;
if(fun==NULL) fun = (WtModel * (*)(double **)) R_FindSymbol("unpack_WtModel_as_double", "ergm", NULL);
return fun(x);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef ENSURE_TH_ORDER
#undef CHOOSE
#undef EQUAL
#undef XOR
#undef XNOR
#undef OUTVAL
#undef INVAL
#undef N_NODES
#undef N_DYADS
#undef OUT_DEG
#undef IN_DEG
#undef DIRECTED
#undef N_EDGES
#undef BIPARTITE
#undef N_TAILS
#undef N_HEADS
#undef NEXT_INEDGE_NUM
#undef NEXT_OUTEDGE_NUM
#undef CHANGE_STAT
#undef N_CHANGE_STATS
#undef INPUT_PARAM
#undef N_INPUT_PARAMS
#undef ZERO_ALL_CHANGESTATS
#undef INPUT_ATTRIB
#undef WtIS_OUTEDGE
#undef WtIS_INEDGE
#undef WtIS_UNDIRECTED_EDGE
#undef WtMIN_OUTEDGE
#undef WtMIN_INEDGE
#undef WtNEXT_OUTEDGE
#undef WtNEXT_INEDGE
#undef WtNEXT_OUTEDGE_PRE
#undef WtNEXT_INEDGE_PRE
#undef OUTWT
#undef INWT
#undef WtSTEP_THROUGH_OUTEDGES
#undef WtSTEP_THROUGH_INEDGES
#undef WtSTEP_THROUGH_OUTEDGES_PRE
#undef WtSTEP_THROUGH_INEDGES_PRE
#undef WtSTEP_THROUGH_OUTEDGES_DECL
#undef WtSTEP_THROUGH_INEDGES_DECL
#undef WtSTEP_THROUGH_OUTEDGES_PRE_DECL
#undef WtSTEP_THROUGH_INEDGES_PRE_DECL
#undef WtEXEC_THROUGH_OUTEDGES
#undef WtEXEC_THROUGH_INEDGES
#undef WtEXEC_THROUGH_EDGES
#undef WtEXEC_THROUGH_OUTEDGES_PRE
#undef WtEXEC_THROUGH_INEDGES_PRE
#undef WtEXEC_THROUGH_EDGES_PRE
#undef WtEXEC_THROUGH_FOUTEDGES
#undef WtEXEC_THROUGH_FINEDGES
#undef WtEXEC_THROUGH_FOUTEDGES_PRE
#undef WtEXEC_THROUGH_FINEDGES_PRE
#undef WtEXEC_THROUGH_NET_EDGES
#undef WtEXEC_THROUGH_NET_EDGES_PRE
#undef WtGETWT
#undef WtSETWT
#undef WtFOR_EACH_TOGGLE
#undef WtTAIL
#undef WtHEAD
#undef WtNEWWT
#undef WtOLDWT
#undef WtGETOLDTOGGLEINFO
#undef WtGETTOGGLEINFO
#undef WtGETNEWTOGGLEINFO
#undef WtSETWT_WITH_BACKUP
#undef WtUNDO_SETWT
#undef WtIF_MORE_TO_COME
#undef WtSETWT_IF_MORE_TO_COME
#undef WtUNDO_PREVIOUS
#undef WtUNDO_PREVIOUS_SETWTS
#undef WtEXEC_THROUGH_TOGGLES
#undef SAMEDYAD
#undef WtGETOLDWT
#undef WtGETNEWWT
#undef WtGETNEWWTOLD
#undef WtC_CHANGESTAT_FN
#undef WtD_CHANGESTAT_FN
#undef WtI_CHANGESTAT_FN
#undef WtU_CHANGESTAT_FN
#undef WtF_CHANGESTAT_FN
#undef WtS_CHANGESTAT_FN
#undef WtD_FROM_S
#undef WtD_FROM_S_FN
#undef IS_OUTEDGE
#undef IS_INEDGE
#undef IS_UNDIRECTED_EDGE
#undef MIN_OUTEDGE
#undef MIN_INEDGE
#undef NEXT_OUTEDGE
#undef NEXT_INEDGE
#undef NEXT_OUTEDGE_PRE
#undef NEXT_INEDGE_PRE
#undef STEP_THROUGH_OUTEDGES
#undef STEP_THROUGH_INEDGES
#undef STEP_THROUGH_OUTEDGES_PRE
#undef STEP_THROUGH_INEDGES_PRE
#undef STEP_THROUGH_OUTEDGES_DECL
#undef STEP_THROUGH_INEDGES_DECL
#undef STEP_THROUGH_OUTEDGES_PRE_DECL
#undef STEP_THROUGH_INEDGES_PRE_DECL
#undef EXEC_THROUGH_OUTEDGES
#undef EXEC_THROUGH_INEDGES
#undef EXEC_THROUGH_EDGES
#undef EXEC_THROUGH_OUTEDGES_PRE
#undef EXEC_THROUGH_INEDGES_PRE
#undef EXEC_THROUGH_EDGES_PRE
#undef EXEC_THROUGH_FOUTEDGES
#undef EXEC_THROUGH_FINEDGES
#undef EXEC_THROUGH_FOUTEDGES_PRE
#undef EXEC_THROUGH_FINEDGES_PRE
#undef EXEC_THROUGH_NET_EDGES
#undef EXEC_THROUGH_NET_EDGES_PRE
#undef GETWT
#undef SETWT
#undef FOR_EACH_TOGGLE
#undef TAIL
#undef HEAD
#undef NEWWT
#undef OLDWT
#undef TOGGLEIND
#undef GETOLDTOGGLEINFO
#undef GETTOGGLEINFO
#undef GETNEWTOGGLEINFO
#undef SETWT_WITH_BACKUP
#undef UNDO_SETWT
#undef IF_MORE_TO_COME
#undef SETWT_IF_MORE_TO_COME
#undef UNDO_PREVIOUS
#undef UNDO_PREVIOUS_SETWTS
#undef EXEC_THROUGH_TOGGLES
#undef GETOLDWT
#undef GETNEWWT
#undef GETNEWWTOLD
#undef D_FROM_S
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef Mweight
#undef WtMH_I_FN
#undef WtMH_U_FN
#undef WtMH_P_FN
#undef WtMH_F_FN
#undef WtFOR_EACH_TERM
#undef WtEXEC_THROUGH_TERMS
#undef WtFOR_EACH_TERM_INREVERSE
#undef WtEXEC_THROUGH_TERMS_INREVERSE
#undef WtEXEC_THROUGH_TERMS_INTO
#undef IFDEBUG_BACKUP_DSTATS
#undef IFDEBUG_RESTORE_DSTATS
#undef WtUPDATE_STORAGE_COND
#undef WtUPDATE_STORAGE
#undef STORAGE
#undef ALLOC_STORAGE
#undef GET_STORAGE
#undef AUX_STORAGE
#undef ALLOC_AUX_STORAGE
#undef GET_AUX_STORAGE
#undef AUX_STORAGE_NUM
#undef GET_AUX_STORAGE_NUM
#undef ALLOC_AUX_SOCIOMATRIX
#undef FREE_AUX_SOCIOMATRIX
#include "R_ext/Rdynload.h"
#include "ergm_wtedgetree.h"
WtNetwork WtNetworkInitialize(Vertex *tails, Vertex *heads, double *weights, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static WtNetwork (*fun)(Vertex *,Vertex *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (WtNetwork (*)(Vertex *,Vertex *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("WtNetworkInitialize", "ergm", NULL);
return fun(tails,heads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
void WtNetworkDestroy(WtNetwork *nwp){
static void (*fun)(WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *)) R_FindSymbol("WtNetworkDestroy", "ergm", NULL);
fun(nwp);
}
WtNetwork WtNetworkInitializeD(double *tails, double *heads, double *weights, Edge nedges,Vertex nnodes, int directed_flag, Vertex bipartite,int lasttoggle_flag, int time, int *lasttoggle){
static WtNetwork (*fun)(double *,double *,double *,Edge,Vertex,int,Vertex,int,int,int *) = NULL;
if(fun==NULL) fun = (WtNetwork (*)(double *,double *,double *,Edge,Vertex,int,Vertex,int,int,int *)) R_FindSymbol("WtNetworkInitializeD", "ergm", NULL);
return fun(tails,heads,weights,nedges,nnodes,directed_flag,bipartite,lasttoggle_flag,time,lasttoggle);
}
WtNetwork * WtNetworkCopy(WtNetwork *dest, WtNetwork *src){
static WtNetwork * (*fun)(WtNetwork *,WtNetwork *) = NULL;
if(fun==NULL) fun = (WtNetwork * (*)(WtNetwork *,WtNetwork *)) R_FindSymbol("WtNetworkCopy", "ergm", NULL);
return fun(dest,src);
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
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef ENSURE_TH_ORDER
#include "R_ext/Rdynload.h"
#include "ergm_wtMHproposal.h"
void WtMH_init(WtMHproposal *MH,char *MHproposaltype, char *MHproposalpackage,double *inputs,int fVerbose,WtNetwork *nwp,void **aux_storage){
static void (*fun)(WtMHproposal *,char *,char *,double *,int,WtNetwork *,void **) = NULL;
if(fun==NULL) fun = (void (*)(WtMHproposal *,char *,char *,double *,int,WtNetwork *,void **)) R_FindSymbol("WtMH_init", "ergm", NULL);
fun(MH,MHproposaltype,MHproposalpackage,inputs,fVerbose,nwp,aux_storage);
}
void WtMH_free(WtMHproposal *MH, WtNetwork *nwp){
static void (*fun)(WtMHproposal *,WtNetwork *) = NULL;
if(fun==NULL) fun = (void (*)(WtMHproposal *,WtNetwork *)) R_FindSymbol("WtMH_free", "ergm", NULL);
fun(MH,nwp);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef ENSURE_TH_ORDER
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef Mweight
#undef WtMH_I_FN
#undef WtMH_U_FN
#undef WtMH_P_FN
#undef WtMH_F_FN
#include "R_ext/Rdynload.h"
#include "ergm_wtmodel.h"
WtModel* WtModelInitialize(char *fnames, char *sonames, double **inputs,int n_terms){
static WtModel* (*fun)(char *,char *,double **,int) = NULL;
if(fun==NULL) fun = (WtModel* (*)(char *,char *,double **,int)) R_FindSymbol("WtModelInitialize", "ergm", NULL);
return fun(fnames,sonames,inputs,n_terms);
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
void WtInitStats(WtNetwork *nwp, WtModel *m){
static void (*fun)(WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,WtModel *)) R_FindSymbol("WtInitStats", "ergm", NULL);
fun(nwp,m);
}
void WtDestroyStats(WtNetwork *nwp, WtModel *m){
static void (*fun)(WtNetwork *,WtModel *) = NULL;
if(fun==NULL) fun = (void (*)(WtNetwork *,WtModel *)) R_FindSymbol("WtDestroyStats", "ergm", NULL);
fun(nwp,m);
}
#undef MIN
#undef MAX
#undef DYADCOUNT
#undef ENSURE_TH_ORDER
#undef GetRandDyad
#undef ENSURE_TH_ORDER
#undef CHOOSE
#undef EQUAL
#undef XOR
#undef XNOR
#undef OUTVAL
#undef INVAL
#undef N_NODES
#undef N_DYADS
#undef OUT_DEG
#undef IN_DEG
#undef DIRECTED
#undef N_EDGES
#undef BIPARTITE
#undef N_TAILS
#undef N_HEADS
#undef NEXT_INEDGE_NUM
#undef NEXT_OUTEDGE_NUM
#undef CHANGE_STAT
#undef N_CHANGE_STATS
#undef INPUT_PARAM
#undef N_INPUT_PARAMS
#undef ZERO_ALL_CHANGESTATS
#undef INPUT_ATTRIB
#undef WtIS_OUTEDGE
#undef WtIS_INEDGE
#undef WtIS_UNDIRECTED_EDGE
#undef WtMIN_OUTEDGE
#undef WtMIN_INEDGE
#undef WtNEXT_OUTEDGE
#undef WtNEXT_INEDGE
#undef WtNEXT_OUTEDGE_PRE
#undef WtNEXT_INEDGE_PRE
#undef OUTWT
#undef INWT
#undef WtSTEP_THROUGH_OUTEDGES
#undef WtSTEP_THROUGH_INEDGES
#undef WtSTEP_THROUGH_OUTEDGES_PRE
#undef WtSTEP_THROUGH_INEDGES_PRE
#undef WtSTEP_THROUGH_OUTEDGES_DECL
#undef WtSTEP_THROUGH_INEDGES_DECL
#undef WtSTEP_THROUGH_OUTEDGES_PRE_DECL
#undef WtSTEP_THROUGH_INEDGES_PRE_DECL
#undef WtEXEC_THROUGH_OUTEDGES
#undef WtEXEC_THROUGH_INEDGES
#undef WtEXEC_THROUGH_EDGES
#undef WtEXEC_THROUGH_OUTEDGES_PRE
#undef WtEXEC_THROUGH_INEDGES_PRE
#undef WtEXEC_THROUGH_EDGES_PRE
#undef WtEXEC_THROUGH_FOUTEDGES
#undef WtEXEC_THROUGH_FINEDGES
#undef WtEXEC_THROUGH_FOUTEDGES_PRE
#undef WtEXEC_THROUGH_FINEDGES_PRE
#undef WtEXEC_THROUGH_NET_EDGES
#undef WtEXEC_THROUGH_NET_EDGES_PRE
#undef WtGETWT
#undef WtSETWT
#undef WtFOR_EACH_TOGGLE
#undef WtTAIL
#undef WtHEAD
#undef WtNEWWT
#undef WtOLDWT
#undef WtGETOLDTOGGLEINFO
#undef WtGETTOGGLEINFO
#undef WtGETNEWTOGGLEINFO
#undef WtSETWT_WITH_BACKUP
#undef WtUNDO_SETWT
#undef WtIF_MORE_TO_COME
#undef WtSETWT_IF_MORE_TO_COME
#undef WtUNDO_PREVIOUS
#undef WtUNDO_PREVIOUS_SETWTS
#undef WtEXEC_THROUGH_TOGGLES
#undef SAMEDYAD
#undef WtGETOLDWT
#undef WtGETNEWWT
#undef WtGETNEWWTOLD
#undef WtC_CHANGESTAT_FN
#undef WtD_CHANGESTAT_FN
#undef WtI_CHANGESTAT_FN
#undef WtU_CHANGESTAT_FN
#undef WtF_CHANGESTAT_FN
#undef WtS_CHANGESTAT_FN
#undef WtD_FROM_S
#undef WtD_FROM_S_FN
#undef IS_OUTEDGE
#undef IS_INEDGE
#undef IS_UNDIRECTED_EDGE
#undef MIN_OUTEDGE
#undef MIN_INEDGE
#undef NEXT_OUTEDGE
#undef NEXT_INEDGE
#undef NEXT_OUTEDGE_PRE
#undef NEXT_INEDGE_PRE
#undef STEP_THROUGH_OUTEDGES
#undef STEP_THROUGH_INEDGES
#undef STEP_THROUGH_OUTEDGES_PRE
#undef STEP_THROUGH_INEDGES_PRE
#undef STEP_THROUGH_OUTEDGES_DECL
#undef STEP_THROUGH_INEDGES_DECL
#undef STEP_THROUGH_OUTEDGES_PRE_DECL
#undef STEP_THROUGH_INEDGES_PRE_DECL
#undef EXEC_THROUGH_OUTEDGES
#undef EXEC_THROUGH_INEDGES
#undef EXEC_THROUGH_EDGES
#undef EXEC_THROUGH_OUTEDGES_PRE
#undef EXEC_THROUGH_INEDGES_PRE
#undef EXEC_THROUGH_EDGES_PRE
#undef EXEC_THROUGH_FOUTEDGES
#undef EXEC_THROUGH_FINEDGES
#undef EXEC_THROUGH_FOUTEDGES_PRE
#undef EXEC_THROUGH_FINEDGES_PRE
#undef EXEC_THROUGH_NET_EDGES
#undef EXEC_THROUGH_NET_EDGES_PRE
#undef GETWT
#undef SETWT
#undef FOR_EACH_TOGGLE
#undef TAIL
#undef HEAD
#undef NEWWT
#undef OLDWT
#undef TOGGLEIND
#undef GETOLDTOGGLEINFO
#undef GETTOGGLEINFO
#undef GETNEWTOGGLEINFO
#undef SETWT_WITH_BACKUP
#undef UNDO_SETWT
#undef IF_MORE_TO_COME
#undef SETWT_IF_MORE_TO_COME
#undef UNDO_PREVIOUS
#undef UNDO_PREVIOUS_SETWTS
#undef EXEC_THROUGH_TOGGLES
#undef GETOLDWT
#undef GETNEWWT
#undef GETNEWWTOLD
#undef D_FROM_S
#undef NO_EDGE
#undef OLD_EDGE
#undef NEW_EDGE
#undef CAN_IGNORE
#undef MAX_TRIES
#undef MH_FAILED
#undef MH_UNRECOVERABLE
#undef MH_IMPOSSIBLE
#undef MH_UNSUCCESSFUL
#undef MH_CONSTRAINT
#undef MH_QUIT_UNSUCCESSFUL
#undef XOR
#undef XNOR
#undef MH_INPUTS
#undef Mtail
#undef Mhead
#undef Mweight
#undef WtMH_I_FN
#undef WtMH_U_FN
#undef WtMH_P_FN
#undef WtMH_F_FN
#undef WtFOR_EACH_TERM
#undef WtEXEC_THROUGH_TERMS
#undef WtFOR_EACH_TERM_INREVERSE
#undef WtEXEC_THROUGH_TERMS_INREVERSE
#undef WtEXEC_THROUGH_TERMS_INTO
#undef IFDEBUG_BACKUP_DSTATS
#undef IFDEBUG_RESTORE_DSTATS
#undef WtUPDATE_STORAGE_COND
#undef WtUPDATE_STORAGE
