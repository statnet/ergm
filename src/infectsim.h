#ifndef INFECTSIM_H 
#define INFECTSIM_H

#include <R.h>
#include "wtedgetree.h"
#include "MCMC.h"

/* Function prototypes */
void InfectSimLoop (int *time, int *N, 
                    int *Nedges, int *discordantedges,
                    int *ND_part1, int *ND_part2,
                    int *ND_startdate, int *ND_enddate,
                    int *R_serostatus, int *R_time_infect,
                    int *R_infector_ID, int *R_infector_recency,
                    int *R_ego_numpartners, int *R_partner_numpartners,
                    int *R_seed, int *R_race,
                    int *R_infector_race,
                    int *lrecvec, double *recvec);

      
#endif

