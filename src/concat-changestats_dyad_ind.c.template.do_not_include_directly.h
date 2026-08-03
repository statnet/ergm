/*  File src/changestats_dyad_ind.c.template.do_not_include_directly.h in
 *  package ergm, part of the Statnet suite of packages for network analysis,
 *  https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_Rutil.h"
#include "ergm_storage.h"
#include "ergm_edgelist.h"

/********************  changestats:  A    ***********/
/*****************
 changestat: c_absdiff
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_absdiff)) {
  double p = INPUT_ATTRIB[0];
  if(p==1.0){
    CHANGE_STAT[0] = ECHANGE(fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]));
  } else {
    CHANGE_STAT[0] = ECHANGE(pow(fabs(INPUT_ATTRIB[tail] - INPUT_ATTRIB[head]), p));
  }
}

/*****************
 changestat: d_absdiffcat
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_absdiffcat)) {
  /* *** don't forget tail -> head */
  double tailval = INPUT_ATTRIB[tail-1],
    headval = INPUT_ATTRIB[head-1],
    absdiff = fabs(tailval - headval);
    if (absdiff>0) {
      for (unsigned int j=0; j<N_CHANGE_STATS; j++) {
        CHANGE_STAT[j] += (absdiff==INPUT_PARAM[j]) ? ECHANGE1 : 0.0;
      }
    }
}

/*****************
 changestat: attrcov
*****************/

typedef struct {
  int *nodecov;
  double **mat;
} SVARIANT(attrcov_storage);

ETYPE(I_CHANGESTAT_FN)(SVARIANT(i_attrcov)) {
  ALLOC_STORAGE(1, SVARIANT(attrcov_storage), sto);
  sto->nodecov = INTEGER(getListElement(mtp->R, "nodecov"));

  int nr = asInteger(getListElement(mtp->R, "nr"));
  int nc = asInteger(getListElement(mtp->R, "nc"));

  // rows vary faster than columns because we did not transpose a$mat in the InitErgmTerm function
  sto->mat = R_Calloc(nc, double *);
  sto->mat[0] = REAL(getListElement(mtp->R, "mat"));
  for(int i = 1; i < nc; i++) {
    sto->mat[i] = sto->mat[i - 1] + nr;
  }
}

ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_attrcov)) {
  GET_STORAGE(SVARIANT(attrcov_storage), sto);
  // head comes before tail here because we did not transpose a$mat in the InitErgmTerm function
  CHANGE_STAT[0] += ECHANGE(sto->mat[sto->nodecov[head]][sto->nodecov[tail]]);
}

ETYPE(F_CHANGESTAT_FN)(SVARIANT(f_attrcov)) {
  GET_STORAGE(SVARIANT(attrcov_storage), sto);
  R_Free(sto->mat);
}

/*****************
 changestat: d_b2cov
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_b2cov)) {
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */
  Vertex nb1 = BIPARTITE;
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[head-nb1+o-1];
	CHANGE_STAT[j] += ECHANGE(sum);
    }
}

/*****************
 changestat: d_b2factor
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_b2factor)) {
  /* *** don't forget tail -> head */
    int headpos = INPUT_ATTRIB[head-1-BIPARTITE];
    if (headpos!=-1) CHANGE_STAT[headpos] += ECHANGE1;
}

/********************  changestats:  D    ***********/

/*****************
 changestat: d_density
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_density)) {
  Dyad ndyads = N_DYADS;

  /* *** don't forget tail -> head */
  CHANGE_STAT[0] += ECHANGE(1.0 / ndyads);
}

/*****************
 changestat: d_diff
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_diff)) {
  double p = INPUT_PARAM[0]; // Conveniently, nodal covariate starts at 1.
  int mul = IINPUT_PARAM[0], sign_code = IINPUT_PARAM[1];

  /* *** don't forget tail -> head */
    double change = (INPUT_PARAM[tail] - INPUT_PARAM[head])*mul;
    switch(sign_code){
    case 1: // identity
      break;
    case 2: // abs
      change = fabs(change);
      break;
    case 3: // positive only
      change = change<0 ? 0 : change;
      break;
    case 4: // negative only
      change = change>0 ? 0 : change;
      break;
    default:
      error("Invalid sign action code passed to d_diff.");
      break;
    }

    if(p==0.0){ // Special case: take the sign of the difference instead.
      change = sign(change);
    }else if(p!=1.0){
      change = pow(change, p);
    }

    CHANGE_STAT[0] += ECHANGE(change);
}

/********************  changestats:  E    ***********/
/*****************
 changestat: d_edgecov
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_edgecov)) {
  double val;
  int nrow, noffset;

  noffset = BIPARTITE;
  if(noffset > 0){
    /*   nrow = (N_NODES)-(long int)(INPUT_PARAM[0]); */
    nrow = noffset;
  }else{
    nrow = (long int)(INPUT_PARAM[0]);
  }

  /* *** don't forget tail -> head */
    /*Get the initial edge state*/
    /*Get the covariate value*/
    val = INPUT_ATTRIB[(head-1-noffset)*nrow+(tail-1)];
    /*  Rprintf("tail %d head %d nrow %d val %f\n", tail, head, nrow, val); */
    /*Update the change statistic, based on the toggle type*/
    CHANGE_STAT[0] += ECHANGE(val);
}

/*****************
 changestat: c_edges
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_edges)) {
  CHANGE_STAT[0] = ECHANGE1;
}

/********************  changestats:  L    ***********/

/*****************
 changestat: d_meandeg
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_meandeg)) {
  // Effectively, change is 2/n if undirected and 1/n if directed.
  if(DIRECTED) CHANGE_STAT[0] = ECHANGE(1.0/N_NODES);
  else CHANGE_STAT[0] = ECHANGE(2.0/N_NODES);
}

/*****************
 changestat: d_mixmat
 General mixing matrix (mm) implementation.
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_mixmat)){
  unsigned int symm = IINPUT_PARAM[0] & 1;
  unsigned int marg = IINPUT_PARAM[0] & 2;
  int *tx = IINPUT_PARAM;
  int *hx = BIPARTITE? IINPUT_PARAM : IINPUT_PARAM + N_NODES;
  int *cells = BIPARTITE? IINPUT_PARAM + N_NODES + 1: IINPUT_PARAM + N_NODES*2 + 1;

  unsigned int diag = tx[tail]==tx[head] && hx[tail]==hx[head];
  for(unsigned int j=0; j<N_CHANGE_STATS; j++){
    unsigned int thmatch = tx[tail]==cells[j*2] && hx[head]==cells[j*2+1];
    unsigned int htmatch = tx[head]==cells[j*2] && hx[tail]==cells[j*2+1];

    int w = DIRECTED || BIPARTITE? thmatch :
      (symm ? thmatch||htmatch : thmatch+htmatch)*(symm && marg && diag?2:1);
    if(w) CHANGE_STAT[j] += ECHANGE(w);
  }
}

/*****************
 changestat: d_nodecov
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodecov)) {
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[tail+o-1] + INPUT_ATTRIB[head+o-1];
	CHANGE_STAT[j] += ECHANGE(sum);
    }
}

/*****************
 changestat: d_nodefactor
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodefactor)) {
  int tailpos = IINPUT_ATTRIB[tail-1];
  int headpos = IINPUT_ATTRIB[head-1];
  if (tailpos!=-1) CHANGE_STAT[tailpos] += ECHANGE1;
  if (headpos!=-1) CHANGE_STAT[headpos] += ECHANGE1;
}

/*****************
 changestat: d_nodeicov
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodeicov)) {
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[head+o-1];
	CHANGE_STAT[j] += ECHANGE(sum);
      }
}

/*****************
 changestat: d_nodeifactor
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodeifactor)) {
  int headpos = INPUT_ATTRIB[head-1];
  if (headpos!=-1) CHANGE_STAT[headpos] += ECHANGE1;
}

/*****************
 changestat: d_nodematch
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodematch)) {
  double matchval;
  Vertex ninputs;
  int j;

  ninputs = N_INPUT_PARAMS - N_NODES;

  /* *** don't forget tail -> head */
    matchval = INPUT_PARAM[tail+ninputs-1];
    if (matchval == INPUT_PARAM[head+ninputs-1]) { /* We have a match! */
      if (ninputs==0) {/* diff=F in network statistic specification */
        CHANGE_STAT[0] += ECHANGE1;
      } else { /* diff=T */
        for (j=0; j<ninputs; j++) {
          if (matchval == INPUT_PARAM[j])
            CHANGE_STAT[j] += ECHANGE1;
        }
      }
    }
}

/*****************
 changestat: nodemix
*****************/

typedef struct {
  int *nodecov;
  int **indmat;
} SVARIANT(nodemix_storage);

ETYPE(I_CHANGESTAT_FN)(SVARIANT(i_nodemix)) {
  ALLOC_STORAGE(1, SVARIANT(nodemix_storage), sto);
  sto->nodecov = INTEGER(getListElement(mtp->R, "nodecov"));

  int nr = asInteger(getListElement(mtp->R, "nr"));
  int nc = asInteger(getListElement(mtp->R, "nc"));

  sto->indmat = R_Calloc(nr, int *);
  sto->indmat[0] = INTEGER(getListElement(mtp->R, "indmat"));
  for(int i = 1; i < nr; i++) {
    sto->indmat[i] = sto->indmat[i - 1] + nc;
  }
}

ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodemix)) {
  GET_STORAGE(SVARIANT(nodemix_storage), sto);
  int index = sto->indmat[sto->nodecov[tail]][sto->nodecov[head]];
  if(index >= 0) {
    CHANGE_STAT[index] += ECHANGE1;
  }
}

ETYPE(F_CHANGESTAT_FN)(SVARIANT(f_nodemix)) {
  GET_STORAGE(SVARIANT(nodemix_storage), sto);
  R_Free(sto->indmat);
}

/*****************
 changestat: d_nodeocov
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodeocov)) {
  unsigned int oshift = N_INPUT_PARAMS / N_CHANGE_STATS;

  /* *** don't forget tail -> head */
      for(unsigned int j=0, o=0; j<N_CHANGE_STATS; j++, o+=oshift){
	double sum = INPUT_ATTRIB[tail+o-1];
	CHANGE_STAT[j] += ECHANGE(sum);
      }
}

/*****************
 changestat: d_nodeofactor
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_nodeofactor)) {
  int tailpos = INPUT_ATTRIB[tail-1];
  if (tailpos!=-1) CHANGE_STAT[tailpos] += ECHANGE1;
}

/********************  changestats:  O    ***********/

/*****************
 changestat: d_receiver
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_receiver)) {
  int j;
  Vertex deg;

  /* *** don't forget tail -> head */
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while((deg != head) && (j < (N_CHANGE_STATS-1))){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(deg==head){CHANGE_STAT[j] += ECHANGE1;}
}

/********************  changestats:  S    ***********/
/*****************
 changestat: d_sender
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_sender)) {
  int j;
  Vertex deg;

  /* *** don't forget tail -> head */
    j=0;
    deg = (Vertex)INPUT_PARAM[j];
    while((deg != tail) && (j < (N_CHANGE_STATS-1))){
      j++;
      deg = (Vertex)INPUT_PARAM[j];
    }
    if(deg==tail){CHANGE_STAT[j] += ECHANGE1;}
}

/*****************
 changestat: d_smalldiff
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_smalldiff)) {

  /* *** don't forget tail -> head */
    CHANGE_STAT[0] += (fabs(INPUT_ATTRIB[tail-1] - INPUT_ATTRIB[head-1])
    > INPUT_PARAM[0]) ? 0.0 : ECHANGE1;
}

/*****************
 changestat: d_sociality
*****************/
ETYPE(C_CHANGESTAT_FN)(SVARIANT(c_sociality)) {
  int j;
  Vertex deg;
  int ninputs, nstats;
  double tailattr;

  ninputs = (int)N_INPUT_PARAMS;
  nstats  = (int)N_CHANGE_STATS;

  /* *** don't forget tail -> head */
  if(ninputs>nstats+1){
    /* match on attributes */
      tailattr = INPUT_ATTRIB[tail-1+nstats+1]; // +1 for the "guard" value between vertex IDs and attribute vector
      if(tailattr == INPUT_ATTRIB[head-1+nstats+1]){
	j=0;
	deg = (Vertex)INPUT_PARAM[j];
	while(deg != tail && j < nstats){
	  j++;
	  deg = (Vertex)INPUT_PARAM[j];
	}
	if(j < nstats){CHANGE_STAT[j] += ECHANGE1;}
	j=0;
	deg = (Vertex)INPUT_PARAM[j];
	while(deg != head && j < nstats){
	  j++;
	  deg = (Vertex)INPUT_PARAM[j];
	}
	if(j < nstats){CHANGE_STAT[j] += ECHANGE1;}
      }

}else{
    /* *** don't forget tail -> head */
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != tail && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats){CHANGE_STAT[j] += ECHANGE1;}
      j=0;
      deg = (Vertex)INPUT_PARAM[j];
      while(deg != head && j < nstats){
	j++;
	deg = (Vertex)INPUT_PARAM[j];
      }
      if(j < nstats){CHANGE_STAT[j] += ECHANGE1;}

}

}

/*****************
 changestat: c_distance
*****************/
C_CHANGESTAT_FN(c_distance) {
  Vertex t,h;
  int j;
  int nv,dim,logd,sphd;
  double dis,dexp,*coord,sphr,logmind,logdoff,scale,mdis,pw;

  /*Set things up*/
  nv = (int)INPUT_PARAM[0];     /*Number of vertices*/
  dim = (int)INPUT_PARAM[1];    /*Dimension of coordinate space*/
  dexp = INPUT_PARAM[2];        /*Minkowski exponent*/
  logd = (int)INPUT_PARAM[3];   /*Should we use the log distance?*/
  sphd = (int)INPUT_PARAM[4];   /*Should we use lat/lon distances on the geosphere?*/
  sphr = INPUT_PARAM[5];        /*Radius to use for spherical coordinates*/
  logmind = INPUT_PARAM[6];     /*Distance minimum for the log case*/
  logdoff = INPUT_PARAM[7];     /*Distance offset for the log case*/
  scale = INPUT_PARAM[8];       /*Scaling adjustment for raw distance*/
  pw = INPUT_PARAM[9];          /*Power to which distance should be raised*/
  coord = INPUT_PARAM+10;       /*Pointer to the coordinate matrix*/

  /*Compute the changescore*/
  t=tail-1;  /*0-indexed values, for C convenience*/
  h=head-1;
  if(sphd){  /*Compute great circle distances on the sphere*/
    dis=spheredist(coord[t],coord[t+nv],coord[h],coord[h+nv],sphr);
  }else{
    dis=0.0;
  }
  /*Compute Minkowski distances in free space (adding extra dimensions to
    great circle distance if sphd and we were given dim>2)*/
  mdis=0.0;
  for(j=0+2*sphd;j<dim;j++){  /*Calculate the t,h distance*/
    if(dexp==1.0){
      mdis+=fabs(coord[t+nv*j]-coord[h+nv*j]);
    }else if(dexp==2.0){
      mdis+=(coord[t+nv*j]-coord[h+nv*j]) * (coord[t+nv*j]-coord[h+nv*j]);
    }else{
      mdis+=pow(fabs(coord[t+nv*j]-coord[h+nv*j]),dexp);
    }
  }
  if(dexp==2.0){
    mdis=sqrt(mdis);
  }else if(dexp!=1.0){
    mdis=pow(mdis,1.0/dexp);
  }
  dis+=mdis;   /*Add the two distance types*/
  dis*=scale;  /*Apply the scaling coefficient*/
  if(pw==0.5)  /*Apply the power transformation*/
    dis=sqrt(dis);
  else if(pw==2.0)
    dis*=dis;
  else if(pw!=1.0)
    dis=pow(dis,pw);
  /*If using log distances, log 'em.  But to prevent divergence, threshold
    from below by logmind (after adding the offset).*/
  if(logd)
    dis=log(MAX(logmind,dis+logdoff));
  /*Check for an edge, and update accordingly*/
  if(DIRECTED)
    CHANGE_STAT[0] += IS_OUTEDGE(tail,head) ? -dis : dis;
  else
    CHANGE_STAT[0] += IS_UNDIRECTED_EDGE(tail,head) ? -dis : dis;
}
