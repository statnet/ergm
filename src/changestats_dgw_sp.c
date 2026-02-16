/*  File src/changestats_dgw_sp.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "changestats_dgw_sp.h"
#include "ergm_storage.h"
#include "changestats.h"

#define all_calcs(term)                         \
  dvec_calc(term)                               \
       dist_calc(term)                          \
       gw_calc(term)

#define all_calcs2(term)                        \
  dvec_calc2(term)                              \
       dist_calc2(term)                         \
       gw_calc2(term)


#define sp_args tail,head,mtp,nwp,edgestate,spcache,N_CHANGE_STATS,dvec,CHANGE_STAT

#define dvec_calc(term)                                                 \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, int nd, Vertex *dvec, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change({                                                   \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += ((L2+echange == deg) - (L2 == deg));                 \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = dvec[j];                                         \
          cs[j] += (echange)*(L2 == deg);                               \
        }                                                               \
      });                                                               \
  }

#define dvec_calc2(term)                                                \
  static inline void term ## _calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, int nd, Vertex *dvec, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change({                                                   \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += ((L2+echange == deg) - (L2 == deg))*2;               \
        }                                                               \
      },{                                                               \
        for(unsigned int j = 0; j < nd; j++){                           \
          Vertex deg = (Vertex)dvec[j];                                 \
          cs[j] += (echange)*(L2 == deg);                               \
        }                                                               \
      });                                                               \
  }


#define spd_args tail,head,mtp,nwp,edgestate,spcache,N_CHANGE_STATS,CHANGE_STAT

#define dist_calc(term)                                                 \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, int nd, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change({                                                   \
        int nL2 = L2 + echange;                                         \
        if(nL2 > nd) cutoff_error(mtp);                                 \
        if(L2) cs[L2-1]--;                                              \
        if(nL2) cs[nL2-1]++;                                            \
      },{                                                               \
        if(L2 > nd) cutoff_error(mtp);                                  \
        if(L2) cs[L2-1] += echange;                                     \
      });                                                               \
  }

#define dist_calc2(term)                                                \
  static inline void term ## _dist_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, int nd, double *cs) { \
    int echange = edgestate ? -1 : 1;                                   \
    term ## _change({                                                   \
        int nL2 = L2 + echange;                                         \
        if(nL2 > nd) cutoff_error(mtp);                                 \
        if(L2) cs[L2-1]-=2;                                             \
        if(nL2) cs[nL2-1]+=2;                                           \
      },{                                                               \
        if(L2 > nd) cutoff_error(mtp);                                  \
        if(L2) cs[L2-1] += echange;                                     \
      });                                                               \
  }


#define gwsp_args tail,head,mtp,nwp,edgestate,spcache,alpha,loneexpa

#define gw_calc(term)                                                   \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, double alpha, double loneexpa) { \
    double cumchange = 0;                                               \
    term ## _change({                                                   \
        cumchange += alpha ? exp(loneexpa*(L2-edgestate)) : L2-edgestate == 0; \
      },{                                                               \
        cumchange += alpha ? exp(alpha + log1mexp(-loneexpa*L2)) : L2 != 0; \
      });                                                               \
    return cumchange;                                                   \
  }


#define gw_calc2(term)                                                  \
  static inline double term ## _gw_calc(Vertex tail, Vertex head, ModelTerm *mtp, Network *nwp, Rboolean edgestate, StoreStrictDyadMapUInt *spcache, double alpha, double loneexpa) { \
    double cumchange = 0;                                               \
    term ## _change({                                                   \
        cumchange += (alpha ? exp(loneexpa*(L2-edgestate)) : L2-edgestate == 0) * 2; \
      },{                                                               \
        cumchange += alpha ? exp(alpha + log1mexp(-loneexpa*L2)) : L2 != 0;    \
      });                                                               \
    return cumchange;                                                   \
  }


all_calcs(dspUTP)
all_calcs(dspOTP)
all_calcs(dspITP)
all_calcs2(dspOSP)
all_calcs2(dspISP)
all_calcs2(dspRTP)

/*****************
 changestat: d_dsp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
C_CHANGESTAT_FN(c_ddsp) { 
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the L2 stats list*/

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP: dspUTP_calc(sp_args); break;
  case L2OTP: dspOTP_calc(sp_args); break;
  case L2ITP: dspITP_calc(sp_args); break;
  case L2RTP: dspRTP_calc(sp_args); break;
  case L2OSP: dspOSP_calc(sp_args); break;
  case L2ISP: dspISP_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_ddspdist) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP: dspUTP_dist_calc(spd_args); break;
  case L2OTP: dspOTP_dist_calc(spd_args); break;
  case L2ITP: dspITP_dist_calc(spd_args); break;
  case L2RTP: dspRTP_dist_calc(spd_args); break;
  case L2OSP: dspOSP_dist_calc(spd_args); break;
  case L2ISP: dspISP_dist_calc(spd_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwdsp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted DSP stats; the specific type of DSP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwdsp) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute log(1-exp(-alpha))*/
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  double cumchange = 0;

  /*Obtain the DSP changescores (by type)*/
  switch(type){
  case L2UTP: cumchange = dspUTP_gw_calc(gwsp_args); break;
  case L2OTP: cumchange = dspOTP_gw_calc(gwsp_args); break;
  case L2ITP: cumchange = dspITP_gw_calc(gwsp_args); break;
  case L2RTP: cumchange = dspRTP_gw_calc(gwsp_args); break;
  case L2OSP: cumchange = dspOSP_gw_calc(gwsp_args); break;
  case L2ISP: cumchange = dspISP_gw_calc(gwsp_args); break;
  }
  
  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
}


all_calcs(espUTP)
     all_calcs(espOTP)
     all_calcs(espITP)
     all_calcs(espOSP)
     all_calcs(espISP)
     all_calcs(espRTP)


/*****************
 changestat: d_esp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
     C_CHANGESTAT_FN(c_desp) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the ESP stats list*/

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP: espUTP_calc(sp_args); break;
  case L2OTP: espOTP_calc(sp_args); break;
  case L2ITP: espITP_calc(sp_args); break;
  case L2RTP: espRTP_calc(sp_args); break;
  case L2OSP: espOSP_calc(sp_args); break;
  case L2ISP: espISP_calc(sp_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_despdist) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP: espUTP_dist_calc(spd_args); break;
  case L2OTP: espOTP_dist_calc(spd_args); break;
  case L2ITP: espITP_dist_calc(spd_args); break;
  case L2RTP: espRTP_dist_calc(spd_args); break;
  case L2OSP: espOSP_dist_calc(spd_args); break;
  case L2ISP: espISP_dist_calc(spd_args); break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwesp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted ESP stats; the specific type of ESP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwesp) { 
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute (1-exp(-alpha))*/
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  double cumchange = 0;

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP: cumchange = espUTP_gw_calc(gwsp_args); break;
  case L2OTP: cumchange = espOTP_gw_calc(gwsp_args); break;
  case L2ITP: cumchange = espITP_gw_calc(gwsp_args); break;
  case L2RTP: cumchange = espRTP_gw_calc(gwsp_args); break;
  case L2OSP: cumchange = espOSP_gw_calc(gwsp_args); break;
  case L2ISP: cumchange = espISP_gw_calc(gwsp_args); break;
  }
  
  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
}


/*****************
 changestat: d_nsp
*****************/
/*
  Note that d_esp is a meta-function, dispatching actual changescore
  calculation to one of the esp*_calc routines, based on the selected shared
  partner type code.

  Type codes are as follows (where (i,j) is the focal edge):

  UTP - Undirected two-path (undirected graphs only)
  OTP - Outgoing two-path (i->k->j)
  ITP - Incoming two-path (i<-k<-j)
  RTP - Reciprocated two-path (i<->k<->j)
  OSP - Outgoing shared partner (i->k<-j)
  ISP - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  UTP should always be used for undirected graphs; OTP is the traditional directed default.
*/
#define NEGATE_CHANGE_STATS for(unsigned int i = 0; i < N_CHANGE_STATS; i++) CHANGE_STAT[i] *= -1;
C_CHANGESTAT_FN(c_dnsp) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  Vertex *dvec = (Vertex*) IINPUT_PARAM+1;           /*Get the pointer to the NSP stats list*/
  
  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP:
    espUTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspUTP_calc(sp_args);
    break;
  case L2OTP:
    espOTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOTP_calc(sp_args);
    break;
  case L2ITP:
    espITP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspITP_calc(sp_args);
    break;
  case L2RTP:
    espRTP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspRTP_calc(sp_args);
    break;
  case L2OSP:
    espOSP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspOSP_calc(sp_args);
    break;
  case L2ISP:
    espISP_calc(sp_args);
    NEGATE_CHANGE_STATS;
    dspISP_calc(sp_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/  
}


C_CHANGESTAT_FN(c_dnspdist) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP:
    espUTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspUTP_dist_calc(spd_args);
    break;
  case L2OTP:
    espOTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOTP_dist_calc(spd_args);
    break;
  case L2ITP:
    espITP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspITP_dist_calc(spd_args);
    break;
  case L2RTP:
    espRTP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspRTP_dist_calc(spd_args);
    break;
  case L2OSP:
    espOSP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspOSP_dist_calc(spd_args);
    break;
  case L2ISP:
    espISP_dist_calc(spd_args);
    NEGATE_CHANGE_STATS;
    dspISP_dist_calc(spd_args);
    break;
  }
  /*We're done!  (Changestats were written in by the calc routine.)*/
}


/*****************
 changestat: d_gwnsp
*****************/

/*
  Note that d_gwesp is a meta-function for all geometrically weighted NSP stats; the specific type of NSP to be employed is determined by the type argument (INPUT_PARAM[1]).  Type codes are as follows (where (i,j) is the focal edge):

  OTP (0) - Outgoing two-path (i->k->j)
  ITP (1) - Incoming two-path (i<-k<-j)
  RTP (2) - Reciprocated two-path (i<->k<->j)
  OSP (3) - Outgoing shared partner (i->k<-j)
  ISP (4) - Incoming shared partner (i<-k->j)

  Only one type may be specified per esp term.  The default, OTP, retains the original behavior of esp/gwesp.  In the case of undirected graphs, OTP should be used (the others assume a directed network memory structure, and are not safe in the undirected case).
*/
C_CHANGESTAT_FN(c_dgwnsp) {
  /*Set things up*/
  StoreStrictDyadMapUInt *spcache = N_AUX ? AUX_STORAGE : NULL;
  double alpha = INPUT_PARAM[0];       /*Get alpha*/
  double loneexpa = log1mexp(alpha);    /*Precompute (1-exp(-alpha))*/
  L2Type type = (L2Type) IINPUT_PARAM[0];     /*Get the L2 type code to be used*/
  double cumchange = 0;

  /*Obtain the changescores (by type)*/
  switch(type){
  case L2UTP:
    cumchange = dspUTP_gw_calc(gwsp_args) - espUTP_gw_calc(gwsp_args);
    break;
  case L2OTP:
    cumchange = dspOTP_gw_calc(gwsp_args) - espOTP_gw_calc(gwsp_args);
    break;
  case L2ITP:
    cumchange = dspITP_gw_calc(gwsp_args) - espITP_gw_calc(gwsp_args);
    break;
  case L2RTP:
    cumchange = dspRTP_gw_calc(gwsp_args) - espRTP_gw_calc(gwsp_args);
    break;
  case L2OSP:
    cumchange = dspOSP_gw_calc(gwsp_args) - espOSP_gw_calc(gwsp_args);
    break;
  case L2ISP:
    cumchange = dspISP_gw_calc(gwsp_args) - espISP_gw_calc(gwsp_args);
    break;
  }

  CHANGE_STAT[0] = edgestate ? -cumchange : cumchange;
}


/*****************
 changestat: c_ddspbwrap
*****************/
C_CHANGESTAT_FN(c_ddspbwrap) {
  c_ddsp(tail, head, mtp, nwp, edgestate);

  // correct for double counting of directed vs. undirected dyads
  for(int ind = 0; ind < N_CHANGE_STATS; ind++) CHANGE_STAT[ind] /= 2.0;
}


C_CHANGESTAT_FN(c_ddspdistbwrap) {
  c_ddspdist(tail, head, mtp, nwp, edgestate);

  // correct for double counting of directed vs. undirected dyads
  for(int ind = 0; ind < N_CHANGE_STATS; ind++) CHANGE_STAT[ind] /= 2.0;
}


/*****************
 changestat: c_dgwdspbwrap
*****************/
C_CHANGESTAT_FN(c_dgwdspbwrap) {
  c_dgwdsp(tail, head, mtp, nwp, edgestate);
  
  // correct for double counting of directed vs. undirected dyads
  CHANGE_STAT[0] /= 2.0;
}

