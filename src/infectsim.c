/* C speedup of infectsim function */

#include "infectsim.h"



/**********************************************************/      

/* In the variables names below,
   ND_xxxx is for network.data
   R_xxxx is for result */
void InfectSimLoop (int *time, int *N, 
                    int *Nedges, int *discordantedges,
                    int *ND_part1, int *ND_part2,
                    int *ND_startdate, int *ND_enddate,
                    int *R_serostatus, int *R_time_infect,
                    int *R_infector_ID, int *R_infector_recency,
                    int *R_ego_numpartners, int *R_partner_numpartners,
                    int *R_seed, int *R_race,
                    int *R_infector_race,
                    int *lrecvec, double *recvec) {
  int ndiscord;
  int i, j, day, edge, R1, R2, P1, P2;
  int id_infect=0, id_suscept=0;
  int rel_duration, day_of_rel, recency;
  double prob_trans;
  
  GetRNGstate();  /* R function enabling uniform RNG */
  for(day=0; day<*time; day++) {
    ndiscord = 0;
    for(j=0; j<*Nedges; j++) {
      if (day <= ND_enddate[j] && day >= ND_startdate[j]) {
        P1 = ND_part1[j];
        P2 = ND_part2[j];
        if (R_serostatus[P1-1] != R_serostatus[P2-1]) {
          discordantedges[ndiscord++] = j+1;
        }
      }
    }
    
    for(i=0; i<ndiscord; i++) {
      edge = discordantedges[i];
      P1 = ND_part1[edge-1];
      P2 = ND_part2[edge-1];
      R1 = R_serostatus[P1-1];
      R2 = R_serostatus[P2-1];
      if (R2==0) {
        id_infect = P1;
        id_suscept = P2;
      } else if (R1==0) {
        id_infect = P2;
        id_suscept = P1;
      } else if (R1>0 && R2>0) {
        if (R_time_infect[P1-1] == day) {
          id_infect = P2;
          id_suscept = P1;
        } else {
          id_infect = P1;
          id_suscept = P2;
        }
      }
      rel_duration = ND_enddate[edge-1] - ND_startdate[edge-1];
      day_of_rel = day - ND_startdate[edge-1];
      recency = day - R_time_infect[id_infect-1];
      prob_trans = recency == 0 ? 0.0 :
                   (recency > *lrecvec? recvec[*lrecvec-1]:
                   recvec[recency-1]);
      if (unif_rand() < prob_trans) {
        R_time_infect[id_suscept-1] = day;
        R_serostatus[id_suscept-1] = 1;
        R_infector_ID[id_suscept-1] = id_infect;
        R_infector_recency[id_suscept-1] = recency;
        R_infector_race[id_suscept-1] = R_race[id_infect-1];
        R_seed[id_suscept-1] = R_seed[id_infect-1];
        R_ego_numpartners[id_suscept-1] = 0;
        R_partner_numpartners[id_suscept-1] = 0;
        for(j=0; j<*Nedges; j++) {
          if (day <= ND_enddate[j] && day >= ND_startdate[j]) {
            P1 = ND_part1[j];
            P2 = ND_part2[j];
            if (P1 == id_suscept || P2 == id_suscept) {
              ++R_ego_numpartners[id_suscept-1];
            }
            if (P1 == id_infect || P2 == id_infect) {
              ++R_partner_numpartners[id_suscept-1];
            }
          }
        }
      }
    }
    
  }
  PutRNGstate();  /* Disable RNG before returning */
}

