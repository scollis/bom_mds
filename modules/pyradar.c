#include "rsl.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int volume_exists(Radar* myradar, int volnum){
   int exists_state;
   Volume *myvolume;
   exists_state=1;
   myvolume=myradar->v[volnum];
   if (myvolume==NULL) exists_state=0;
   return exists_state;
}

Radar *get_radar(char* str) {
  Radar *myradar;
  printf("We are loading %s \n", str);
  myradar = RSL_anyformat_to_radar(str, NULL);
  return myradar;
}

int get_nvolumes(Radar* myradar){
  int nvols;
  nvols=myradar->h.nvolumes;
   return nvols;
}

int free_radar(Radar* myradar){
   int dummy;
   dummy=1;
   RSL_free_radar(myradar);
   return dummy;
}

Radar_header get_radar_header(Radar* myradar){
   Radar_header my_header;
   my_header= myradar->h;
   return my_header;
}

Volume *get_volume(Radar* myradar, int volnum){
   Volume *myvolume;
   myvolume=myradar->v[volnum];
   return myvolume;
}

int get_nsweeps(Volume* myvolume){
   int nsweeps;
   nsweeps=myvolume->h.nsweeps;
   return nsweeps;
}

float get_lat(Radar* myradar){
  float latnum;
  latnum= myradar->h.latd + myradar->h.latm/60.0 + myradar->h.lats/(60.0*60.0);
  return latnum;
}

float get_lon(Radar* myradar){
  float lonnum;
  lonnum= myradar->h.lond + myradar->h.lonm/60.0 + myradar->h.lons/(60.0*60.0);
  return lonnum;
}


Sweep *get_sweep(Volume* myvolume, int sweepnum){
   Sweep *mysweep;
   mysweep=myvolume->sweep[sweepnum];
   return mysweep;
}

Ray *get_ray(Sweep* mysweep, int raynum){
   Ray *myray;
   myray=mysweep->ray[raynum];
   return myray;
}

float get_value(Ray* myray, int gatenum){
   float myvalue;
   myvalue=myray->h.f(myray->range[gatenum]);
   return myvalue;
}

float get_value_debug(Ray* myray, int gatenum){
   float myvalue;
   myvalue=myray->h.f(myray->range[gatenum]);
   printf("The value we are getting is %f \n", myvalue);
   return myvalue;
}

float get_storage(Ray* myray, int gatenum){
   float myvalue;
   myvalue=myray->range[gatenum];
   /*printf("The value we are getting is %f \n", myvalue);*/
   return myvalue;
}

int get_ngates(Ray* myray){
   int nbins;
   nbins=myray->h.nbins;
   return nbins;
}

int get_first_gate(Ray* myray){
   int bin1;
   bin1=myray->h.range_bin1;
   return bin1;
}

int get_gate_size(Ray* myray){
   int gate_size;
   gate_size=myray->h.gate_size;
   return gate_size;
}

int get_nrays(Sweep* mysweep){
   int nrays;
   nrays=mysweep->h.nrays;
   return nrays;
}

float get_azimuth(Ray* myray){
   float azimuth;
   azimuth=myray->h.azimuth;
   return azimuth;
}

float get_elev(Ray* myray){
   float elev;
   elev=myray->h.elev;
   return elev;
}

int get_yr(Ray* myray){
   int yr;
   yr=myray->h.year;
   return yr;
}
int get_month(Ray* myray){
   int mo;
   mo=myray->h.month;
   return mo;

}
int get_day(Ray* myray){
   int day;
   day=myray->h.day;
   return day;

}
int get_hour(Ray* myray){
   int hr;
   hr=myray->h.hour;
   return hr;

}
int get_min(Ray* myray){
   int min;
   min=myray->h.minute;
   return min;

}
float get_sec(Ray* myray){
   float sec;
   sec=myray->h.sec;
   return sec;
}

char* get_radarname(Radar* myradar){
   char *nm =(char *) malloc(64);
   strncpy(nm, myradar->h.radar_name, 64);
   return nm;
}



