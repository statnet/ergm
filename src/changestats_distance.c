
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
