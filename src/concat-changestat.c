
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
