#include <cpp/ergm_wtnetwork.h>
#include <cpp/ergm_wtchangestat.h>

/*****************
 stat: cyclicalweights
*****************/

WtC_CHANGESTAT_CPP(cyclicalweights, {
    unsigned int path = mt.dattrib[0];
    unsigned int combine = mt.dattrib[1];
    unsigned int compare =  mt.dattrib[2];

    // (tail,head) as the focus dyad
    // This means that the strongest 2-path doesn't change.
    double two_paths = 0;
    for(auto [k, yhk]: nw.out_neighbors(head)) {
      double two_path;

      switch(path){
      case 1: two_path = fmin(nw(k,tail), yhk); break; // min
      case 2: two_path = sqrt(nw(k,tail) * yhk); break; // geomean
      default: // never reached, but prevents a warning
        two_path = 0;
      }

      switch(combine){
      case 1: two_paths = fmax(two_paths, two_path); break; // max
      case 2: two_paths += two_path; break; // sum
      }
    }

    switch(compare){
    case 1: mt.stat[0] += fmin(two_paths, weight) - fmin(two_paths, edgestate); break; // min
    case 2: mt.stat[0] += sqrt(two_paths * weight) - sqrt(two_paths * edgestate); break; // geomean
    }

    // (tail,head) as the first link in the 2-path
    // This means that only the strongest 2-path may change.
    for(auto [j, yjt]: nw.in_neighbors(tail)) {
      if(j==head) continue;

      double old_two_paths = 0;
      double new_two_paths = 0;

      for(auto [k, ykj]: nw.in_neighbors(j)) {
        double old_ytk = (k==head) ? edgestate : nw(tail, k);
        double new_ytk = (k==head) ? weight : old_ytk;
        double old_two_path;
        double new_two_path;

        switch(path){
        case 1: // min
          old_two_path = fmin(old_ytk, ykj);
          new_two_path = fmin(new_ytk, ykj);
          break;
        case 2: // geomean
          old_two_path = sqrt(old_ytk * ykj);
          new_two_path = sqrt(new_ytk * ykj);
          break;
        default: // never reached, but prevents a warning
          old_two_path = 0;
          new_two_path = 0;
        }

        switch(combine){
        case 1:
          old_two_paths = fmax(old_two_paths, old_two_path);// max
          new_two_paths = fmax(new_two_paths, new_two_path);
          break; // max
        case 2:
          old_two_paths += old_two_path;
          new_two_paths += new_two_path;
          break; // sum
        }

      }

      switch(compare){
      case 1: mt.stat[0] += fmin(new_two_paths, yjt) - fmin(old_two_paths, yjt); break; // min
      case 2: mt.stat[0] += sqrt(new_two_paths * yjt) - sqrt(old_two_paths * yjt); break; // geomean
      }
    }

    // (tail,head) as the second link of the 2-path
    // This means that only the strongest 2-path may change.
    for(auto [i, yhi]: nw.out_neighbors(head)) {
      if(i==tail) continue;

      double old_two_paths = 0;
      double new_two_paths = 0;

      for(auto [k, yik]: nw.out_neighbors(i)) {
        double old_ykh = (k==tail) ? edgestate : nw(k,head);
        double new_ykh = (k==tail) ? weight : old_ykh;
        double old_two_path;
        double new_two_path;

        switch(path){
        case 1: // min
          old_two_path = fmin(old_ykh, yik);
          new_two_path = fmin(new_ykh, yik);
          break;
        case 2: // geomean
          old_two_path = sqrt(old_ykh * yik);
          new_two_path = sqrt(new_ykh * yik);
          break;
        default: // never reached, but prevents a warning
          old_two_path = 0;
          new_two_path = 0;
        }

        switch(combine){
        case 1:
          old_two_paths = fmax(old_two_paths, old_two_path);// max
          new_two_paths = fmax(new_two_paths, new_two_path);
          break; // max
        case 2:
          old_two_paths += old_two_path;
          new_two_paths += new_two_path;
          break; // sum
        }
      }

      switch(compare){
      case 1: mt.stat[0] += fmin(new_two_paths, yhi) - fmin(old_two_paths, yhi); break; // min
      case 2: mt.stat[0] += sqrt(new_two_paths * yhi) - sqrt(old_two_paths * yhi); break; // geomean
      }
    }
  });


WtS_CHANGESTAT_CPP(cyclicalweights, {
    unsigned int path = mt.dattrib[0];
    unsigned int combine = mt.dattrib[1];
    unsigned int compare =  mt.dattrib[2];

    /* unsigned int threshold = mt.dattrib[3]; */

    for (auto [tail, head, yth]: nw.edges()) {
      double two_paths = 0;
      for(auto [node3, yh3]: nw.out_neighbors(head)) {
        double two_path;

        switch(path){
        case 1: two_path = fmin(nw(node3,tail), yh3); break; // min
        case 2: two_path = sqrt(nw(node3,tail) * yh3); break; // geomean
        default: // never reached, but prevents a warning
          two_path = 0;
        }

        switch(combine){
        case 1: two_paths = fmax(two_paths, two_path); break; // max
        case 2: two_paths += two_path; break; // sum
        }
      }

      switch(compare){
      case 1: mt.stat[0] += fmin(two_paths, yth); break; // min
      case 2: mt.stat[0] += sqrt(two_paths * yth); break; // geomean
      }

    }
  });
