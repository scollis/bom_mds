###################################################################
#Pickle_to_cappi.py: Take a python pickle of a radar scan and generate a CAPPI and save it as a netcdf file
###################################################################
#                                                                 
###################################################################
# Scott Collis, CAWCR, May 2008                                   #
###################################################################
__author__ = "Scott Collis, s.collis@bom.gov.au"
__version__ = "1.0"

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
from pylab import num2date, datestr2num
import parse_ini
import read_radar

def std_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %date_dict

def radar_to_cappi(radar1_filename, radar2_filename,**kwargs):
	ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
	loud=kwargs.get('loud', False)
	ini_dict=parse_ini.parse_ini(ini_fname)
	if loud: print "Loading radar file 1"
	if 'radar1_path' in ini_dict.keys():
		radar1=read_radar.load_radar(ini_dict['radar1_path']+radar1_filename)
	else:
		radar1=pyradar.load_radar(radar1_filename)
	if loud: print "Loading radar file 2"
	if 'radar2_path' in ini_dict.keys():
		radar2=read_radar.load_radar(ini_dict['radar2_path']+radar2_filename)
	else:
		radar2=pyradar.load_radar(radar2_filename)
	cappi_z_bounds=ini_dict.get('cappi_z_bounds', [500,15000])
	cappi_xy_bounds=ini_dict.get('cappi_xy_bounds', [-50000, 50000])
	cappi_resolution=ini_dict.get('cappi_resolution', [100, 40])
	levs=linspace(cappi_z_bounds[0], cappi_z_bounds[1], cappi_resolution[1])
	xar=linspace(cappi_xy_bounds[0], cappi_xy_bounds[1], cappi_resolution[0])
	yar=linspace(cappi_xy_bounds[0], cappi_xy_bounds[1], cappi_resolution[0])
	displace=mathematics.corner_to_point(radar1[0]['radar_loc'], radar2[0]['radar_loc'])
	if loud: print "Cappi-ing radar 1"
	radar1_cube=radar_to_cart.make_cube(radar1, xar, yar, levs)
	if loud: print "Cappi-ing radar 2"
	radar2_cube=radar_to_cart.make_cube(radar2, xar, yar, levs, displacement=displace)
	cube_fname=ini_dict['cube_path']+'cappi_'+std_datestr(radar1_cube['date'])+'.nc'
	netcdf_utis.save_data_cube(radar1_cube, radar2_cube, cube_fname)


