#include <cpp/ergm_network.h>
#include <cpp/ergm_changestat.h>

C_CHANGESTAT_CPP(triangle, {
    int change = 0;
    if(mt.iinput.size()){
      bool diff = mt.iinput[0] != 0;
      unsigned int tailattr = mt.iinput[tail];
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
