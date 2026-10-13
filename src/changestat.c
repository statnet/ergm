/*  File src/changestat.c in package ergm, part of the Statnet suite of
 *  packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free, open
 *  source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2003-2026 Statnet Commons
 */
#include "ergm_changestat.h"

/*****************
 my_choose:
 Simple routine to return simple binomial coefficients quickly, 
 avoiding costly call to choose() function.  Note:  my_choose is
 usually not called directly; use CHOOSE macro instead.
*****************/
double my_choose(double n, int r) {
  const double recip_factorial[21] = {1.0, 1.0, 0.5,
	 1.66666666666667e-01, 4.16666666666667e-02, 8.33333333333333e-03,
         1.38888888888889e-03, 1.98412698412698e-04, 2.48015873015873e-05,
         2.75573192239859e-06, 2.75573192239859e-07, 2.50521083854417e-08,
         2.08767569878681e-09, 1.60590438368216e-10, 1.14707455977297e-11,
         7.64716373181982e-13, 4.77947733238739e-14, 2.81145725434552e-15,
         1.56192069685862e-16, 8.22063524662433e-18, 4.11031762331216e-19};
  double ans;

  if (r>20)
    return choose (n, (double)r); /* Use complicated function for large r */
  for(ans=recip_factorial[r]; r>0; r--)
    ans*=(n--);
  return ans;
}

/*
Utility function to calculate distances on a sphere of radius r, between two points
in spherical coordinates.  Our notation betrays our assumption that these are
lat/lon points on the geosphere, but technically you could use any sphere you
wanted, by passing an alternative radius (so long as your coordinates were given
in angular units).
*/
double spheredist(double lat0, double lon0, double lat1, double lon1, double r) {
  double rlat0,rlat1,rlon0,rlon1,cosrl0,cosrl1,cosdl,sinrl0,sinrl1,sindl;

  rlat0=lat0/180.0*M_PI;
  rlat1=lat1/180.0*M_PI;
  rlon0=lon0/180.0*M_PI;
  rlon1=lon1/180.0*M_PI;
  cosrl0=cos(rlat0);
  cosrl1=cos(rlat1);
  cosdl=cos(rlon0-rlon1);
  sinrl0=sin(rlat0);
  sinrl1=sin(rlat1);
  sindl=sin(rlon0-rlon1);
  return r*atan2(sqrt((cosrl0*sindl) * (cosrl0*sindl) + (cosrl1*sinrl0-sinrl1*cosrl0*cosdl) * (cosrl1*sinrl0-sinrl1*cosrl0*cosdl)), sinrl0*sinrl1+cosrl0*cosrl1*cosdl);
}
