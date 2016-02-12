import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
import simul_winds
import propigation	
import gracon_vel2d
import gracon_vel2d_3d
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

def speed_test():
	ber_loc=[-12.457, 130.925]
	gp_loc=[-12.2492,  131.0444]
	use_guess='sonde'
	sonde_file='/bm/gdata/scollis/twpice/darwin.txt'
	tim='1350'
	ber, gp=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/20060122_1350_winds_ver1.nc')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	u=gp['u_array']
	v=gp['v_array']
	dx=dy=gp['xar'][1]-gp['xar'][1]
	tim_date=num2date(datestr2num('20060122 '+tim))
	sonde_list=read_sounding.read_sounding_within_a_day(sonde_file, tim_date)
	launch_dates=[sonde['date_list'][0] for sonde in sonde_list]
	launch_date_offset=[date2num(sonde['date_list'][0])- date2num(tim_date)  for sonde in sonde_list]
	best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[0]]
	print 'Time of radar: ', tim_date, ' Time of sonde_launch: ', best_sonde['date_list'][0], ' Time of sonde_termination: ', best_sonde['date_list'][-1]
	interp_sonde=read_sounding.interp_sounding(best_sonde, gp['levs'])
	tdry=interp_sonde['tdry(degs)']
	pressure=interp_sonde['press(hPa)']
	w=u*0.0
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	mywts=met.make_mask(ber, gp, angs, 5.0, 80.0)
	submask=met.make_submask(mywts)
	t0=systime()
	gv_w, poo=grad_conj_solver_3d.continuity_cost2(u*0.0, u, v,w, dx, dx, gp['levs'], pressure, tdry, mywts, submask)
	print (poo*submask).mean()
	print "total_time=", systime()-t0

if __name__ == "__main__":
	speed_test()
