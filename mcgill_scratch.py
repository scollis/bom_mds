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
import cappi_v2


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



kwargs={}
datestr='200601220700'

#check to see if we have the radar files
#check to see if there are deailased files
#kwargs={}
loud=kwargs.get('loud', False)
ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
dateobj=num2date(datestr2num(datestr))
ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
ini_dict=parse_ini.parse_ini(ini_fname)
radar1_deal_list=os.listdir(ini_dict['radar1_path'])
radar2_deal_list=os.listdir(ini_dict['radar2_path'])
radar1_raw_list=os.listdir(ini_dict['radar1_raw_path'])
radar2_raw_list=os.listdir(ini_dict['radar2_raw_path'])

radar1_target=ini_dict['radar1_prefix']+std_datestr(dateobj, ini_dict['radar1_type'])
radar2_target=ini_dict['radar2_prefix']+std_datestr(dateobj, ini_dict['radar2_type'])

poss_deal_files1=[]
for item in radar1_deal_list:
    if radar1_target in item: poss_deal_files1.append(item)

if len(poss_deal_files1)==0:
    poss_raw_files1=[]
    print "No dealiased files found... Dealiasing"
    for item in radar1_raw_list:
        if radar1_target in item: poss_raw_files1.append(item)
    if len(poss_raw_files1)==0:
        #print "no files found"
        raise IOError, 'Radar 2 File not there'
        #return
    else:
        print "Dealiasing "+poss_raw_files1[0]
        radar1_filename=dealias.dealias_arb(poss_raw_files1[0], ini_dict['radar1_type'], ini_dict['radar1_raw_path'], ini_dict['radar1_path'], ini_dict['radar1_prefix'])
else:
    radar1_filename=poss_deal_files1[0]


    poss_deal_files2=[]
    for item in radar2_deal_list:
        if radar2_target in item: poss_deal_files2.append(item)

    if len(poss_deal_files2)==0:
        poss_raw_files2=[]
        print "No dealiased files found... Dealiasing"
        for item in radar2_raw_list:
            if radar2_target in item: poss_raw_files2.append(item)
        if len(poss_raw_files2)==0:
            #print "no files found"
            raise IOError, 'Radar 2 File not there'
        else:
            radar2_filename=dealias.dealias_arb(poss_raw_files2[0], ini_dict['radar2_type'], ini_dict['radar2_raw_path'], ini_dict['radar2_path'], ini_dict['radar2_prefix'])
    else:
        radar2_filename=poss_deal_files2[0]



if loud: print "Loading radar file 1"

	
if 'radar1_path' in ini_dict.keys():
    radar1=read_radar.load_radar(ini_dict['radar1_path']+radar1_filename)
else:
    radar1=read_radar.load_radar(radar1_filename)


if loud: print "Loading radar file 2"

	
if 'radar2_path' in ini_dict.keys():
    radar2=read_radar.load_radar(ini_dict['radar2_path']+radar2_filename)
else:
    radar2=read_radar.load_radar(radar2_filename)



	pres.plot_ppi(radar2[2],'VE', fig_path='/scratch/bom_mds_dumps/', fig_name='radar2_ve.png')
	pres.plot_ppi(radar2[2],'CZ', fig_path='/scratch/bom_mds_dumps/', fig_name='radar2_cz.png')
	cappi_z_bounds=ini_dict.get('cappi_z_bounds', [500,15000])
	cappi_xy_bounds=ini_dict.get('cappi_xy_bounds', [-50000, 50000])
	cappi_resolution=ini_dict.get('cappi_resolution', [100, 40])
	levs=linspace(cappi_z_bounds[0], cappi_z_bounds[1], cappi_resolution[1])
	xar=linspace(cappi_xy_bounds[0], cappi_xy_bounds[1], cappi_resolution[0])
	yar=linspace(cappi_xy_bounds[0], cappi_xy_bounds[1], cappi_resolution[0])
	displace=mathematics.corner_to_point(radar1[0]['radar_loc'], radar2[0]['radar_loc'])
	if loud: print "Cappi-ing radar 1"
	
	#radar1_cube_=radar_to_cart.make_cube(radar1, xar, yar, levs)
	radar1_cube=cappi_v2.make_cube_all(radar1,xar, yar,levs)
	#max_el=array([scan['Elev'][0] for scan in radar1]).max()
	#radar1_cube=cappi_v2.blend(radar1_cube_v,radar1_cube_h, max_el,loud=True)
	
	if loud: print "Cappi-ing radar 2"
	
	#radar2_cube_v=radar_to_cart.make_cube(radar2, xar, yar, levs, displacement=displace)
radar2_cube=cappi_v2.make_cube_all(radar2,xar, yar,levs, displacement=displace)
	#max_el=array([scan['Elev'][0] for scan in radar2]).max()
	#radar2_cube=cappi_v2.blend(radar2_cube_v,radar2_cube_h, max_el,loud=True)
	#radar2_cube_v=radar_to_cart.make_cube(radar2, xar, yar, levs, displacement=displace)
cube_fname=ini_dict['cube_path']+'cappi_'+std_datestr(radar1_cube['date'], "uf")+'.nc'
	#netcdf_utis.save_data_cube(radar1_cube, radar2_cube, cube_fname)
	
	#Initial Guess



cube_fname=ini_dict['cube_path']+'cappi_'+std_datestr(radar1_cube['date'], "uf")+'.nc'
req=[ 'alt(m)',  'wspd(m/s)',  'wdir(degs)', 'tdry(degs)','press(hPa)' ]
first_sonde,second_sonde = read_sounding.get_two_best_conc_sondes(datestr, req_vars=req)
interp_sonde=read_sounding.interp_sonde_time(first_sonde, second_sonde, dateobj, levs)

if ini_dict['initial_guess']=='sonde':
	#using a sonde for out initial gues
	u_ig=ones(radar1_cube['CZ'].shape, dtype=float)
	v_ig=ones(radar1_cube['CZ'].shape, dtype=float)
	w_ig=zeros(radar1_cube['CZ'].shape, dtype=float)
	for k in range(len(levs)):
		u_ig[:,:,k]=1.0*u_ig[:,:,k]*interp_sonde['wspd(m/s)'][k]*sin(pi*interp_sonde['wdir(degs)'][k]/180.0)
		v_ig[:,:,k]=1.0*v_ig[:,:,k]*interp_sonde['wspd(m/s)'][k]*cos(pi*interp_sonde['wdir(degs)'][k]/180.0)
else:
	u_ig=zeros(radar1_cube['CZ'].shape, dtype=float)
	v_ig=zeros(radar1_cube['CZ'].shape, dtype=float)
	w_ig=zeros(radar1_cube['CZ'].shape, dtype=float)

Re=6371.0*1000.0
rad_at_radar=Re*sin(pi/2.0 -abs(radar1[0]['radar_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
lons=radar1[0]['radar_loc'][1]+360.0*xar/(rad_at_radar*2.0*pi)
lats=radar1[0]['radar_loc'][0] + 360.0*yar/(Re*2.0*pi)	
#Masking
angs=array(propigation.make_lobe_grid(radar2[0]['radar_loc'], radar1[0]['radar_loc'], lats,lons))
mywts=met.make_mask_bad1(radar2_cube, radar1_cube, angs, 1.0, 80.0)

print "Mean gp masked Velocity ", (radar1_cube['VE']*mywts).mean()
print "min gp masked Velocity ", (radar1_cube['VE']*mywts).min()
print "max gp masked Velocity ", (radar1_cube['VE']*mywts).max()
print "Mean Berrimah masked Velocity ", (radar2_cube['VE']*mywts).mean()
print "min Berrimah masked Velocity ", (radar2_cube['VE']*mywts).min()
print "max Berrimah masked Velocity ", (radar2_cube['VE']*mywts).max()
print "Mean gp masked CZ ", (radar1_cube['CZ']*mywts).mean()
print "min gp masked CZ ", (radar1_cube['CZ']*mywts).min()
print "max gp masked CZ ", (radar1_cube['CZ']*mywts).max()
print "Mean Berrimah masked CZ ", (radar2_cube['CZ']*mywts).mean()
print "min Berrimah masked CZ ", (radar2_cube['CZ']*mywts).min()
print "max Berrimah masked CZ ", (radar2_cube['CZ']*mywts).max()
print "Number of masked points", (mywts.shape[0]*mywts.shape[1]*mywts.shape[2])-mywts.sum()
print "Number of unmasked points ", mywts.sum()
print "**********************FALLSPEED INFO****************************"
#def terminal_velocity(refl, temps, levs, display=False):
tdry=interp_sonde['tdry(degs)']
pressure=interp_sonde['press(hPa)']
dummy=met.terminal_velocity(radar1_cube['CZ']*mywts, tdry, radar1_cube['levs'], display=True)
print "**********************FALLSPEED INFO****************************"
f=0.0
X=[u_ig,v_ig,w_ig]
G,F,X=grad_conj_solver_3d.gracon_3d_packaged(X ,radar2_cube, radar1_cube, mywts, interp_sonde)
u_array,v_array,w_array=X 
radar1_cube.update({'u_array':u_array, 'v_array':v_array, 'w_array':w_array})
netcdf_utis.save_data_cube(radar1_cube, radar2_cube,  '/data/cube_data/'+std_datestr(dateobj, "uf") +'_winds.nc')
	











