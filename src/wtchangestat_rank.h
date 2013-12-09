#ifndef WTCHANGESTAT_RANK_H
#define WTCHANGESTAT_RANK_H

#include "wtedgetree.h"
#include "wtchangestat.h"


// Case 1: One ego swaps the values of two alters.
// Case 2: The alters are adjacent?
// Otherwise, run the toggle loop with code wherein one ego changes the value of one alter.

#define GETOLDWT2(a,b) (SAMEDYAD(TAIL,HEAD1,a,b) ? OLDWT1 : (SAMEDYAD(TAIL,HEAD2,a,b) ? OLDWT2 : GETWT(a,b)))
#define GETNEWWT2(a,b) (SAMEDYAD(TAIL,HEAD1,a,b) ? NEWWT1 : (SAMEDYAD(TAIL,HEAD2,a,b) ? NEWWT2 : GETWT(a,b)))
#define GETNEWWT2OLD(a,b,old) (SAMEDYAD(TAIL,HEAD1,a,b) ? NEWWT1 : (SAMEDYAD(TAIL,HEAD2,a,b) ? NEWWT2 : (old)))

#define OPTIMAL_RANK_D(case1sub,defaultsub){				\
    if(ntoggles==2 && tails[0]==tails[1]){				\
      Vertex TAIL=tails[0], HEAD1=heads[0], HEAD2=heads[1];		\
      double OLDWT1=GETWT(TAIL, HEAD1), OLDWT2=GETWT(TAIL, HEAD2);	\
      double NEWWT1=weights[0], NEWWT2=weights[1];			\
      if(NEWWT1==OLDWT2 && NEWWT2==OLDWT1){				\
	ZERO_ALL_CHANGESTATS();						\
	{case1sub};							\
      }									\
    }else{								\
      EXEC_THROUGH_TOGGLES({defaultsub});				\
    }									\
  };									

#endif
