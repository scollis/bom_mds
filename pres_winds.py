import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
import simul_winds
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
from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random
import pickle_zip
import parse_ini

def std_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %date_dict




def pres_winds(datestr,  **kwargs): #latstr, lonstr,
	ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
	ini_dict=parse_ini.parse_ini(ini_fname)
	parm=kwargs.get('parm','CZ')
	dateobj=num2date(datestr2num(datestr))
	tim_date=num2date(datestr2num(datestr))
	radar1, radar2=netcdf_utis.load_cube('/data/cube_data/'+std_datestr(tim_date)+'_winds.nc')
	print radar1['radar_name']
	print radar2['radar_name']
	lvl=kwargs.get('lvl',3000.0)
	lvl_num=argsort(abs(radar1['levs']-lvl))[0]
	print "Level_num=", lvl_num
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar1['radar_loc'][0]*pi/180.0))
	lons=radar1['radar_loc'][1]+360.0*radar1['xar']/(rad_at_radar*2.0*pi)
	lats=radar1['radar_loc'][0] + 360.0*radar1['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=[-12.2492,  131.0444]
	angs=array(propigation.make_lobe_grid(radar2['radar_loc'], radar1['radar_loc'], lats,lons))
	mywts=met.make_mask_bad1(radar2, radar1, angs, 1.0, 80.0)
	if ini_dict['cross'][0]=='max' or ini_dict['cross'][1]=='max':
		 #maskedcz=radar1['CZ'][:,:,lvl_num]*mywts[:,:,lvl_num]
		 #i,j=mathematics.where_closest_2d(maskedcz.max(), maskedcz)
		 maskedw=radar1['w_array'][:,:,7]*mywts[:,:,7]
		 i,j=mathematics.where_closest_2d(maskedw.max(), maskedw)
		 print i,j
		 lat=lats[j[0]]
		 lon=lons[i[0]]
		 print "Max w at ", lat, lon 
	else:
		 lon=ini_dict['cross'][1]#float(lonstr)
		 lat=ini_dict['cross'][0]#float(latstr)
	f=figure()
	alat, alon, alvl=pres.plot_slices(lat, lon, lvl, radar1, lats, lons, radar1['levs'], radar1['u_array'], radar1['v_array'], radar1['w_array'], angs, mywts, par=parm, w_mag=1.0,box=ini_dict['pres_box'], bquiver=[0.05, 0.75], ksp=0.05,qscale=ini_dict['qscale'])
	t1='Gunn Point reflectivity (dBZ) and reconstructed winds (m/s)\n sliced at %(alat)2.2fS and %(alon)3.2fE and %(alvl)d Metres on  %(day)02d/%(mon)02d/%(yr)04d at ' %{'day':tim_date.day, 'mon':tim_date.month, 'yr':tim_date.year,'alat':abs(alat), 'alon':alon, 'alvl':alvl}
	t2=" %(HH)02d%(MM)02dZ" %{'HH':tim_date.hour, 'MM':tim_date.minute}
	f.text( .1, .92, t1+t2) 
	inte_part=1000*(float(int(lat))-lat)
	print  {'alat':abs(alat), 'alon':alon, 'alvl':alvl}
	#ff=os.getenv('HOME')+'/bom_mds/output/recons_'+std_datestr(tim_date)[0:-5]+'/slicer3_%(alat)2.02f_%(alon)3.02f_%(alvl)05d_' %{'alat':abs(alat), 'alon':alon, 'alvl':alvl}
	ff=os.getenv('HOME')+'/bom_mds/output/tests/slicer_'+std_datestr(tim_date)+'_%(alat)2.02f_%(alon)3.02f_%(alvl)05d_' %{'alat':abs(alat), 'alon':alon, 'alvl':alvl}
	print ff
	savefig(ff+t2+'.png', dpi=200)
	close(f)
	#f=figure()
	#mo=pres.reconstruction_plot_pcolor(lats, lons, radar1, lvl_num, parm, radar1['u_array'][:,:,lvl_num], radar1['v_array'][:,:,lvl_num],radar1['w_array'][:,:,lvl_num], angs, mywts[:,:,lvl_num],w_mag=1.0,box=ini_dict['pres_box'], bquiver=[0.05, 0.75], ksp=0.05, qscale=ini_dict['qscale'])
	#mydir=os.getenv('HOME')+'/bom_mds/output/tests/pcolor'
	#savefig(mydir+t2+'.png')
	#reconstruction_plot_pcolor(lats, lons, data_cube, lvl_num, parm, u, v,w, angs, mask, **kwargs):

if __name__ == "__main__":
	pres_winds(sys.argv[1], parm=sys.argv[2], lvl=float(sys.argv[3]))#, sys.argv[2],sys.argv[3])


