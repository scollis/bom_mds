#make stats on a wind cube
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
#import simul_winds
import propigation	
#import gracon_vel2d
#import gracon_vel2d_3d
import read_rays
import radar_math
import radar_to_cart
import mathematics
import netcdf_utis
#import grad_conj_solver_plus
import grad_conj_solver_plus_plus
import grad_conj_solver_3d
import read_sounding
import met
#from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random
#import pickle_zip
import dealias
import parse_ini
import read_radar
import matplotlib.numerix.ma as M

def std_datestr(date_obj, datetype):
	if datetype=="lassen":
		mydatestr=std_lass_datestr(date_obj)
	elif datetype=="uf":
		mydatestr=std_uf_datestr(date_obj)
	return mydatestr

def std_uf_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %date_dict

def std_lass_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d%(HH)02d%(MM)02d" %date_dict

def cube_stats(date_str):
	#date_str='20060123 1600'
	tim_date=num2date(datestr2num(date_str))
	radar1, radar2=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/'+std_datestr(tim_date, 'uf')+'_winds.nc')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar1['radar_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar1['radar_loc'][1]+360.0*radar1['xar']/(rad_at_radar*2.0*pi)
	lats=radar1['radar_loc'][0] + 360.0*radar1['yar']/(Re*2.0*pi)	
	angs=array(propigation.make_lobe_grid(radar2['radar_loc'], radar1['radar_loc'], lats,lons))
	mywts=met.make_mask_bad1(radar2, radar1, angs, 1.0, 80.0)
	submask=met.make_submask(mywts)
	X=[radar1['u_array'], radar1['v_array'], radar1['w_array']]
	req=[ 'alt(m)',  'wspd(m/s)',  'wdir(degs)', 'tdry(degs)','press(hPa)' ]
	first_sonde,second_sonde = read_sounding.get_two_best_conc_sondes(date_str, req_vars=req)
	interp_sonde=read_sounding.interp_sonde_time(first_sonde, second_sonde, tim_date, radar1['levs'])
	costs=grad_conj_solver_3d.return_cost(X, radar2, radar1, mywts, interp_sonde, submask,loud=True)
	disag=costs[1]/(mywts.sum()*2.0)
	print disag

def cube_stats_f(date_str):
	#date_str='20060123 1600'
	tim_date=num2date(datestr2num(date_str))
	radar1, radar2=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/'+std_datestr(tim_date, 'uf')+'_winds.nc')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar1['radar_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar1['radar_loc'][1]+360.0*radar1['xar']/(rad_at_radar*2.0*pi)
	lats=radar1['radar_loc'][0] + 360.0*radar1['yar']/(Re*2.0*pi)	
	angs=array(propigation.make_lobe_grid(radar2['radar_loc'], radar1['radar_loc'], lats,lons))
	mywts=met.make_mask_bad1(radar2, radar1, angs, 1.0, 80.0)
	submask=met.make_submask(mywts)
	X=[radar1['u_array'], radar1['v_array'], radar1['w_array']]
	req=[ 'alt(m)',  'wspd(m/s)',  'wdir(degs)', 'tdry(degs)','press(hPa)' ]
	first_sonde,second_sonde = read_sounding.get_two_best_conc_sondes(date_str, req_vars=req)
	interp_sonde=read_sounding.interp_sonde_time(first_sonde, second_sonde, tim_date, radar1['levs'])
	costs=grad_conj_solver_3d.return_cost(X, radar2, radar1, mywts, interp_sonde, submask,loud=True)
	disag=costs[1]/(mywts.sum()*2.0)
	print disag


if __name__ == "__main__":
	cube_stats(sys.argv[1])#, sys.argv[2],sys.argv[3])


