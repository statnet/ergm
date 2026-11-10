#include <cpp/ergm_network.h>
#include <cpp/ergm_changestat.h>

struct DegcorStorage {
  double degree2;
  double degree3;
  double nedges;
  double cross;
};

double degcor_statistic(double degree2, double degree3, double cross, double nedges) {
  if(nedges <= 0) return 0;
  double mean = degree2 / (2.0 * nedges);
  double sigma2 = degree3 / (2.0 * nedges) - mean * mean;
  return sigma2 > 0 ? (cross / (2.0 * nedges) - mean * mean) / sigma2 : 0;
}

double degcor_cross_change(ergm::ErgmCppNetwork& nw, Vertex tail, Vertex head, Rboolean edgestate) {
  double change = 0;
  double taildeg = nw.degree(tail);
  double headdeg = nw.degree(head);
  if(edgestate) {
    change -= taildeg * headdeg;
    for(auto node: nw.neighbors(head))
      if(node != tail) change -= nw.degree(node);
    for(auto node: nw.neighbors(tail))
      if(node != head) change -= nw.degree(node);
  }else{
    change += (taildeg + 1.0) * (headdeg + 1.0);
    for(auto node: nw.neighbors(head)) change += nw.degree(node);
    for(auto node: nw.neighbors(tail)) change += nw.degree(node);
  }
  return 2.0 * change;
}

I_CHANGESTAT_CPP(degcor, DegcorStorage, {
    mt.storage = R_Calloc(1, DegcorStorage);
    mt.storage->nedges = nw.ptr->nedges;
    for(auto node: nw.nodes()) {
      double degree = nw.degree(node);
      mt.storage->degree2 += degree * degree;
      mt.storage->degree3 += degree * degree * degree;
    }
    for(auto [edge_tail, edge_head]: nw.edges())
      mt.storage->cross += 2.0 * nw.degree(edge_tail) * nw.degree(edge_head);
  })

C_CHANGESTAT_CPP(degcor, DegcorStorage, {
    double echange = edgestate ? -1.0 : 1.0;
    double taildeg = nw.degree(tail);
    double headdeg = nw.degree(head);
    double degree2 = mt.storage->degree2 + 2.0 * echange * (taildeg + headdeg) + 2.0;
    double degree3 = mt.storage->degree3 + 3.0 * echange * (taildeg * taildeg + headdeg * headdeg)
      + 3.0 * (taildeg + headdeg) + 2.0 * echange;
    double nedges = mt.storage->nedges + echange;
    double cross = mt.storage->cross + degcor_cross_change(nw, tail, head, edgestate);
    mt.stat[0] += degcor_statistic(degree2, degree3, cross, nedges)
      - degcor_statistic(mt.storage->degree2, mt.storage->degree3, mt.storage->cross, mt.storage->nedges);
  })

U_CHANGESTAT_CPP(degcor, DegcorStorage, {
    double change = edgestate ? -1.0 : 1.0;
    double taildeg = nw.degree(tail);
    double headdeg = nw.degree(head);
    mt.storage->degree2 += 2.0 * change * (taildeg + headdeg) + 2.0;
    mt.storage->degree3 += 3.0 * change * (taildeg * taildeg + headdeg * headdeg)
      + 3.0 * (taildeg + headdeg) + 2.0 * change;
    mt.storage->nedges += change;
    mt.storage->cross += degcor_cross_change(nw, tail, head, edgestate);
  })

C_CHANGESTAT_CPP(triangle, {
    int change = 0;
    if(mt.iinput.size()){
      bool diff = mt.iinput[0] != 0;
      int tailattr = mt.iinput[tail];
      if(tailattr && tailattr == mt.iinput[head]){
        for(auto k: nw.neighbors(head)) {
          if(tailattr == mt.iinput[k]){
            if (nw.dir) change += nw(k, tail) + nw(tail, k);
            else change += nw(k,tail);
          }
        }
        mt.stat[diff ? tailattr-1 : 0] += edgestate ? -change : change;
      }
    }else{ /* no attribute matching */
      for(auto k: nw.neighbors(head)) {
        if (nw.dir) change += nw(k, tail) + nw(tail, k);
        else change += nw(k,tail);
      }
      mt.stat[0] += edgestate ? -change : change;
    }
  })


/*****************
 changestat: d_cycle
*****************/
void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode,
                           Vertex *visited, Vertex curlen, Vertex *countv, Vertex maxlen, bool semi);
void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head,
                           Vertex *countv, Vertex maxlen, bool semi);

I_CHANGESTAT_CPP(cycle, Vertex, {
    mt.storage = R_Calloc(mt.iinput[1] * 2, Vertex);
  })

C_CHANGESTAT_CPP(cycle, Vertex, {
    int emult;

  /*Perform initial setup*/
  bool semi = (bool) mt.iinput[0];             /*Are we using semicycles?*/
  Vertex maxlen = (Vertex) mt.iinput[1];      /*Get max cycle length*/

  /* *** don't forget tail -> head */
  /*Clear out the count vector*/
  memset(mt.storage, 0, sizeof(*mt.storage)*(maxlen-1));
    /*In semi-cycle case, this toggle can't matter if there is a*/
    /*head->tail edge in the graph; not counting saves much time.*/
    if(!(semi&&(IS_OUTEDGE(head,tail)))){
      /*Count the cycles associated with this edge*/
      edgewise_cycle_census(nwp, tail, head, mt.storage, maxlen, semi);

      /*Make the change, as needed*/
      emult = edgestate ? -1 : 1;
      for(unsigned int j=0, k=0; j<maxlen-1;j++)
        if(mt.iinput[2+j] > 0)
          mt.stat[k++] += emult * (int) mt.storage[j];
    }
  })

/*****************
 edgewise_path_recurse:  Called by d_cycle
*****************/
void edgewise_path_recurse(Network *nwp, Vertex dest, Vertex curnode,
                           Vertex *visited, Vertex curlen, Vertex *countv, Vertex maxlen, bool semi) {
  Vertex v;
  Edge e;

  /*If we've found a path to the destination, increment the census vector*/
  if(DIRECTED){  /*Use outedges, or both if counting semi-paths*/
    if(!semi)
      countv[curlen] += IS_OUTEDGE(curnode, dest);
    else
      countv[curlen] += (IS_OUTEDGE(curnode, dest) || IS_INEDGE(curnode, dest));
  }else{   /*For undirected graphs, edges go from low to high*/
    if(curnode<dest)
      countv[curlen] += IS_OUTEDGE(curnode, dest);
    else
      countv[curlen] += IS_INEDGE(curnode, dest);
  }

  /*If possible, keep searching for novel paths*/
  if(curlen<maxlen-2){
    visited[curlen+1]=curnode; /*Add current node to visited list*/

    /*Recurse on all unvisited neighbors of curnode*/
    STEP_THROUGH_OUTEDGES(curnode,e,v){
      bool rflag = true;
      for(Vertex i=0;(i<=curlen)&&(rflag);i++)  /*Check earlier nodes in path*/
        rflag=(v!=visited[i]);
      if(rflag)
        edgewise_path_recurse(nwp,dest,v,visited,curlen+1,countv,maxlen, semi);
    }
    if(semi||(!DIRECTED)){ /*If semi or !directed, need in-neighbors too*/
      STEP_THROUGH_INEDGES(curnode,e,v){
        bool rflag = ((!DIRECTED)||(!(IS_OUTEDGE(curnode,v))));
        for(Vertex i=0;(i<=curlen)&&(rflag);i++)  /*Check earlier nodes in path*/
          rflag=(v!=visited[i]);
        if(rflag)
          edgewise_path_recurse(nwp,dest,v,visited,curlen+1,countv,maxlen, semi);
      }
    }
  }
}

/*****************
 edgewise_cycle_census:  Called by d_cycle
*****************/
void edgewise_cycle_census(Network *nwp, Vertex tail, Vertex head,
                           Vertex *countv, Vertex maxlen, bool semi) {
  /* *** don't forget tail -> head */
  Vertex *visited;
  Vertex v;
  Edge e;

  /*First, check for a 2-cycle (but only if directed and !semi)*/
  if(DIRECTED && (!semi) && IS_OUTEDGE(head,tail))
    countv[0]++;
  if(N_NODES == 2)
    return;                 /*Failsafe for graphs of order 2*/

  /*Perform the recursive path count*/
  visited=countv + maxlen; /*Locate the list of visited nodes*/
  memset(visited, 0, sizeof(*visited)*maxlen);
  visited[0]=tail;
  visited[1]=head;

  /*Recurse on each neighbor of head*/
  STEP_THROUGH_OUTEDGES(head,e,v){
    if(v!=tail)
      edgewise_path_recurse(nwp,tail,v,visited,1,countv,maxlen,semi);
  }
  if(semi||(!DIRECTED)){ /*If semi or !directed, need in-neighbors too*/
    STEP_THROUGH_INEDGES(head,e,v){
      if((v!=tail)&&((!DIRECTED)||(!(IS_OUTEDGE(head,v)))))
        edgewise_path_recurse(nwp,tail,v,visited,1,countv,maxlen, semi);
    }
  }
}
