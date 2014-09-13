/*  File src/conddegmple.c in package ergm, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "conddegmple.h"

/*****************
 changestat: d_conddegmple
*****************/
D_CHANGESTAT_FN(d_conddegmple) 
{
  
//Rprintf("ntoggles %d\n", ntoggles);
//Rprintf("tails[0] %d tails[1] %d tails[2] %d tails[3] %d\n", tails[0],tails[1],tails[2],tails[3]);
//Rprintf("heads[0] %d heads[1] %d heads[2] %d heads[3] %d\n", heads[0],heads[1],heads[2],heads[3]);

  if(ntoggles==4){
   CHANGE_STAT[0] = (tails[0] < tails[3]) ? 2.0 : 1.0;
  }else{
   CHANGE_STAT[0] = 0.0;
  }
//Rprintf("CHANGE_STAT[0] %f\n", CHANGE_STAT[0]);
}
