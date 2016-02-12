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
import read_rays
import radar_to_cart
import mathematics
import netcdf_utis
from time import time as systime
from numpy import linspace, array, arctan, pi, random



def save_cube_test(date,gpnum, bernum):
	gp_0740=read_rays.construct_lassen_scan(path='/bm/gscratch/scollis/gunn_pt/'+date+gpnum+'/')
	ber_0740=read_rays.construct_uf_scan(path='/bm/gscratch/scollis/berrimah/'+date+'_'+bernum+'/')	
	print gp_0740[0].keys
	print gp_0740[0]['CZ'][0,:]
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	displace=mathematics.corner_to_point(gp_loc, ber_loc)
	ldict={'lat_0':gp_loc[0], 'lon_0':gp_loc[1],'llcrnrlat':-13.0, 'llcrnrlon':130.2, 'urcrnrlat':-12.0 , 'urcrnrlon':131.2, 'lat_ts':gp_loc[0]}
	levs=linspace(500,10000, 30)
	xar=linspace(-50.,50., 100)*1000.0
	yar=linspace(-50.,50., 100)*1000.0
	gp_cube=radar_to_cart.make_cube(gp_0740, xar, yar, levs)
	ber_cube=radar_to_cart.make_cube(ber_0740, xar-displace[0], yar-displace[1], levs)
	print gp_cube['CZ'].shape
	netcdf_utis.save_data_cube(ber_cube, gp_cube, '/bm/gdata/scollis/cube_data/'+date+'_'+bernum[0:4]+'_ver2.nc', gp_loc)


if __name__ == "__main__":
	t0=systime()
	print "the uber cool test"
	print sys.argv
	save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	#simple_reconstruction()
	#test_pert_winds()
	print "Finished running runtime=",systime()-t0, "Seconds"

