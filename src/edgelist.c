#include "edgelist.h"
/*********************
 unsigned int dEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented in edgelist.h using an array of double.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

unsigned int dEdgeListSearch(Vertex head, Vertex tail, double *el){
  unsigned int nedges=el[0];
  unsigned int u=nedges,l=1,m;

  do{
    m = l + (u-l)/2;
    if(tail>el[nedges+m] || (tail==el[nedges+m] && head>el[m])) l = m+1;
    else u = m-1;
  }while(!(tail==el[nedges+m] && head==el[m]) && l<=u);

  if(l>u) return(0);
  else return(m);
}

/*********************
 unsigned int iEdgeListSearch

 A function to check whether a given edge is in an edgelist formatted
 as documented in edgelist.h using an array of int.

 Returns the index of the edge (counting from 1) if found, 0 if not.
*********************/

unsigned int iEdgeListSearch(Vertex head, Vertex tail, int *el){
  unsigned int nedges=el[0];
  unsigned int u=nedges,l=1,m;

  do{
    m = l + (u-l)/2;
    if(tail>el[nedges+m] || (tail==el[nedges+m] && head>el[m])) l = m+1;
    else u = m-1;
  }while(!(tail==el[nedges+m] && head==el[m]) && l<=u);

  if(l>u) return(0);
  else return(m);
}
