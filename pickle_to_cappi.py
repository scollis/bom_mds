###################################################################
#Pickle_to_cappi.py: Take a python pickle of a radar scan and generate a CAPPI and save it as a netcdf file
###################################################################
#                                                                 
###################################################################
# Scott Collis, CAWCR, May 2008                                   #
###################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from matplotlib import use as use_interface
# chose a non-GUI backend
use_interface( 'Agg' )
#import simul_winds
import propigation	
#import gracon_vel2d
import radar_to_cart
import mathematics
import netcdf_utis
from time import time as systime
from numpy import linspace, array, arctan, pi, random
import pickle_zip
from pylab import num2date, datestr2num

def std_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %date_dict

def pickle_to_cappi(gp_pickle, ber_pickle, **kwargs):
	path=kwargs.get('path', '/bm/gkeep/scollis/deal_ber/')
	debug=kwargs.get('debug', False)
	if debug: print "Loading Pickles"
	gp=pickle_zip.load(path+gp_pickle)
	ber=pickle_zip.load(path+ber_pickle)	
	ber_loc=[-12.457, 130.925]
	gp_loc=[-12.2492,  131.0444]
	displace=mathematics.corner_to_point(gp_loc, ber_loc)
	ldict={'lat_0':gp_loc[0], 'lon_0':gp_loc[1],'llcrnrlat':-13.0, 'llcrnrlon':130.2, 'urcrnrlat':-12.0 , 'urcrnrlon':131.2, 'lat_ts':gp_loc[0]}
	levs=linspace(500,10000, 30)
	xar=linspace(-50.,50., 100)*1000.0
	yar=linspace(-50.,50., 100)*1000.0
	gp_cube=radar_to_cart.make_cube(gp, xar, yar, levs)
	ber_cube=radar_to_cart.make_cube(ber, xar-displace[0], yar-displace[1], levs)
	print gp_cube['CZ'].shape
	netcdf_utis.save_data_cube(ber_cube, gp_cube, '/bm/gdata/scollis/cube_data/'+std_datestr(gp[0]['date'])+'_deal.nc', gp_loc)

def pick_a_pickle(datestr, **kwargs):
	path=kwargs.get('path', '/bm/gkeep/scollis/deal_ber/')
	target_date=num2date(datestr2num(datestr))
	target_cpol='C-POL_deal_'+std_datestr(target_date)+'.pickle.gz'
	target_ber='Berrimah_deal_'+std_datestr(target_date)+'.pickle.gz'
	
	if not(target_cpol in os.listdir(path)):
		raise IOError, 'No time found for Cpol, need '+ target_cpol
	if not(target_ber in os.listdir(path)):
		raise IOError, 'No time found for Berrimah, need '+ target_ber
	pickle_to_cappi(target_cpol, target_ber, debug=True)



if __name__ == "__main__":
	t0=systime()
	print "the uber cool test"
	print sys.argv
	pick_a_pickle(sys.argv[1])
	#save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	#simple_reconstruction()
	#test_pert_winds()
	print "Finished running runtime=",systime()-t0, "Seconds"

