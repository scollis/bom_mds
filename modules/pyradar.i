%module pyradar
%import include/rsl.h
%{
#include "rsl.h"
extern int volume_exists(Radar* myradar, int volnum);
extern Radar *get_radar(char* str);
extern int get_nvolumes(Radar* myradar);
extern int free_radar(Radar* myradar);
extern Radar_header get_radar_header(Radar* myradar);
extern Volume *get_volume(Radar* myradar, int volnum);
extern Sweep *get_sweep(Volume* myvolume, int sweepnum);
extern Ray *get_ray(Sweep* mysweep, int raynum);
extern float get_value(Ray* myray, int gatenum);
extern float get_value_debug(Ray* myray, int gatenum);
extern float get_storage(Ray* myray, int gatenum);
extern int get_ngates(Ray* myray);
extern int get_first_gate(Ray* myray);
extern int get_gate_size(Ray* myray);
extern int get_nrays(Sweep* mysweep);
extern float get_azimuth(Ray* myray);
extern float get_elev(Ray* myray);
extern int get_yr(Ray* myray);
extern int get_month(Ray* myray);
extern int get_day(Ray* myray);
extern int get_hour(Ray* myray);
extern int get_min(Ray* myray);
extern float get_sec(Ray* myray);
extern char* get_radarname(Radar* myradar);
extern int get_nsweeps(Volume* myvolume);
extern float get_lon(Radar* myradar);
extern float get_lat(Radar* myradar);
%}
%include pyradar.c
